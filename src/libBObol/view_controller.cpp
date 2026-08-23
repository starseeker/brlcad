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
#include "BObol/BExportAction.h"
#include "BObol/BHUDLabelOverlay.h"
#include "BObol/BLineLayerOverlay.h"
#include "BObol/BLodService.h"
#include "BObol/BLodUpdateAction.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BMeshLodSubmitAction.h"
#include "BObol/BMeshResidencyAction.h"
#include "BObol/BMeshShape.h"
#include "BObol/BPickDetail.h"
#include "BObol/BSceneController.h"
#include "BObol/BSnapAction.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "BObol/BViewLod.h"
#include "BObol/BViewStore.h"
#include "cad_assembly_private.h"
#include "lod_coordinator_private.h"
#include "retained_allocation_private.h"
#include "view_controller_private.h"
#include "raytrace.h"
#include "rt/view.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <queue>
#include <string.h>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Inventor/SbName.h>
#include <Inventor/SbPlane.h>
#include <Inventor/SoDB.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>
#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/gl.h>
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

static const char *controller_database_id(const struct db_i *dbip);
std::vector<SoBRLDatabaseSource *> controller_render_database_sources(
    const BObolViewController *controller);
std::vector<SoBRLDatabaseSource *> controller_render_database_source_roots(
    const BObolViewController *controller);
static void controller_synchronize_compact_cad_presentations(
    BObolViewController *controller);
static std::vector<SoBRLMeshShape *> controller_render_mesh_shapes(
    const BObolViewController *controller);

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

/* Kept out of the public controller layout so callback teardown is ABI-safe. */
struct ControllerFrameRequestState {
    ControllerFrameRequestState(BObolFrameRequestCallback callback_in,
	void *user_data_in) :
	callback(callback_in),
	userData(user_data_in),
	dispatches(0),
	closing(false)
    {
    }

    bool beginDispatch(void)
    {
	std::lock_guard<std::mutex> lock(this->mutex);
	if (this->closing)
	    return false;
	this->dispatches++;
	return true;
    }

    void finishDispatch(void)
    {
	std::lock_guard<std::mutex> lock(this->mutex);
	this->dispatches--;
	if (!this->dispatches)
	    this->cv.notify_all();
    }

    void close(void)
    {
	std::unique_lock<std::mutex> lock(this->mutex);
	this->closing = true;
	this->cv.wait(lock, [this]() {
	    return this->dispatches == 0;
	});
    }

    BObolFrameRequestCallback callback;
    void *userData;
    std::mutex mutex;
    std::condition_variable cv;
    unsigned int dispatches;
    bool closing;
};

static void
controller_initialize_render_action(SoRenderManager *manager)
{
    SoGLRenderAction *action = manager ? manager->getGLRenderAction() : NULL;
    if (!action)
	return;
    /* Rendering providers are host/controller policy, never process-global
     * fallback.  A host binds its concrete manager before it can render. */
    action->setContextManager(NULL);
    action->setCacheContext(SoGLCacheContextElement::getUniqueCacheContext());
}

struct ControllerGLFunctions {
    void (*clearColor)(GLclampf, GLclampf, GLclampf, GLclampf) = NULL;
    void (*clear)(GLbitfield) = NULL;
    void (*getIntegerv)(GLenum, GLint *) = NULL;
    void (*pushAttrib)(GLbitfield) = NULL;
    void (*popAttrib)(void) = NULL;
    void (*enable)(GLenum) = NULL;
    void (*disable)(GLenum) = NULL;
    void (*depthMask)(GLboolean) = NULL;
    void (*matrixMode)(GLenum) = NULL;
    void (*pushMatrix)(void) = NULL;
    void (*popMatrix)(void) = NULL;
    void (*loadIdentity)(void) = NULL;
    void (*begin)(GLenum) = NULL;
    void (*end)(void) = NULL;
    void (*color3f)(GLfloat, GLfloat, GLfloat) = NULL;
    void (*vertex2f)(GLfloat, GLfloat) = NULL;

    void load(SoDB::ContextManager *m)
    {
#define CONTROLLER_GL_LOAD(member, name) \
	member = reinterpret_cast<decltype(member)>(m->getProcAddress(name))
	CONTROLLER_GL_LOAD(clearColor, "glClearColor");
	CONTROLLER_GL_LOAD(clear, "glClear");
	CONTROLLER_GL_LOAD(getIntegerv, "glGetIntegerv");
	CONTROLLER_GL_LOAD(pushAttrib, "glPushAttrib");
	CONTROLLER_GL_LOAD(popAttrib, "glPopAttrib");
	CONTROLLER_GL_LOAD(enable, "glEnable");
	CONTROLLER_GL_LOAD(disable, "glDisable");
	CONTROLLER_GL_LOAD(depthMask, "glDepthMask");
	CONTROLLER_GL_LOAD(matrixMode, "glMatrixMode");
	CONTROLLER_GL_LOAD(pushMatrix, "glPushMatrix");
	CONTROLLER_GL_LOAD(popMatrix, "glPopMatrix");
	CONTROLLER_GL_LOAD(loadIdentity, "glLoadIdentity");
	CONTROLLER_GL_LOAD(begin, "glBegin");
	CONTROLLER_GL_LOAD(end, "glEnd");
	CONTROLLER_GL_LOAD(color3f, "glColor3f");
	CONTROLLER_GL_LOAD(vertex2f, "glVertex2f");
#undef CONTROLLER_GL_LOAD
    }

    bool complete(void) const
    {
	return clearColor && clear && getIntegerv && pushAttrib && popAttrib &&
	    enable && disable && depthMask && matrixMode && pushMatrix &&
	    popMatrix && loadIdentity && begin && end && color3f && vertex2f;
    }
};

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
    this->coordinatorPhase = BOBOL_LOD_COORDINATOR_FALLBACK;
    this->coordinatorEvent = BOBOL_LOD_EVENT_INITIALIZE;
    this->coordinatorTransitionSerial = 0;
    this->coordinatorProgressSequence = 0;
    this->coordinatorDispatchSerial = 0;
    this->coordinatorStagnantDispatchCount = 0;
    this->coordinatorInvariantViolationCount = 0;
    this->coordinatorInvariantMask = 0;
    this->coordinatorInvariantHistoryMask = 0;
    this->viewQualityHistoryEntryCount = 0;
    this->viewQualityHistoryRememberCount = 0;
    this->viewQualityHistoryRecallCount = 0;
    this->viewRevision = 0;
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
    this->residentMeshBytes = 0;
    this->stableResidentMeshBytes = 0;
    this->reservedResidentMeshGrowthBytes = 0;
    this->residentMeshLimitBytes = 0;
    this->memoryLimitedPayloadCount = 0;
    this->activeWorkingSetBytes = 0;
    this->peakWorkingSetBytes = 0;
    this->residentCompactionCount = 0;
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
    this->viewReady = FALSE;
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

SbMatrix
bobol_sbmatrix_from_brl_mat(const mat_t mat)
{
    if (!mat)
	return SbMatrix::identity();

    return SbMatrix(
	       static_cast<float>(mat[0]),  static_cast<float>(mat[4]),
	       static_cast<float>(mat[8]),  static_cast<float>(mat[12]),
	       static_cast<float>(mat[1]),  static_cast<float>(mat[5]),
	       static_cast<float>(mat[9]),  static_cast<float>(mat[13]),
	       static_cast<float>(mat[2]),  static_cast<float>(mat[6]),
	       static_cast<float>(mat[10]), static_cast<float>(mat[14]),
	       static_cast<float>(mat[3]),  static_cast<float>(mat[7]),
	       static_cast<float>(mat[11]), static_cast<float>(mat[15]));
}

SbRotation
bobol_camera_orientation_from_brl_rotation(const mat_t rotation)
{
    if (!rotation)
	return SbRotation::identity();

    SbMatrix cameraAxes = SbMatrix::identity();
    cameraAxes[0][0] = static_cast<float>(rotation[0]);
    cameraAxes[0][1] = static_cast<float>(rotation[1]);
    cameraAxes[0][2] = static_cast<float>(rotation[2]);
    cameraAxes[1][0] = static_cast<float>(rotation[4]);
    cameraAxes[1][1] = static_cast<float>(rotation[5]);
    cameraAxes[1][2] = static_cast<float>(rotation[6]);
    cameraAxes[2][0] = static_cast<float>(rotation[8]);
    cameraAxes[2][1] = static_cast<float>(rotation[9]);
    cameraAxes[2][2] = static_cast<float>(rotation[10]);
    return SbRotation(cameraAxes);
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

static void
controller_collect_mesh_shapes(SoNode *node,
			       std::vector<SoBRLMeshShape *> &shapes)
{
    if (!node)
	return;

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	if (std::find(shapes.begin(), shapes.end(), shape) == shapes.end())
	    shapes.push_back(shape);
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	controller_collect_mesh_shapes(group->getChild(i), shapes);
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

/*
 * Apply resident compact-occurrence presentation deltas before Coin begins a
 * render traversal.  SoBRLDatabaseSource retains a defensive lazy update in
 * GLRender/getBoundingBox for non-controller users, but mutating a CAD batch
 * from inside its own traversal is one frame too late: the current traversal
 * has already selected the old instance set.  Edit hide/restore consequently
 * produced a blank or stale frame even though no geometry had been evicted.
 *
 * This is not realization or LoD planning.  The source journal makes the
 * normal path proportional to the changed occurrence set (often one entry),
 * and the immutable resident part payload remains attached to the same CAD
 * assembly.  Source roots are sufficient because a compact registry owns the
 * descendants it represents; stopping there also avoids a hierarchy walk per
 * occurrence in large scenes.
 */
static void
controller_synchronize_compact_cad_presentations(
    BObolViewController *controller)
{
    if (!controller)
	return;
    BObolViewLodState *viewState = controller->getViewLodState();
    if (!viewState)
	return;

    const std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_source_roots(controller);
    for (SoBRLDatabaseSource *source : sources) {
	if (!source || !source->isCompactOccurrenceRegistry())
	    continue;
	if (source->currentCompactViewLodAssembly(viewState))
	    continue;
	const std::vector<const BObolViewLodState::CadPayload *> noPayloads;
	(void)source->compactViewLodAssembly(noPayloads, viewState);
    }
}

static std::vector<SoBRLMeshShape *>
controller_render_mesh_shapes(const BObolViewController *controller)
{
    std::vector<SoBRLMeshShape *> shapes;
    if (!controller)
	return shapes;

    SoNode *root = controller->getRenderSceneRoot();
    if (!root)
	root = controller->getSceneRoot();
    controller_collect_mesh_shapes(root, shapes);
    return shapes;
}

static SbString
rt_pick_result_path(const BObolRtPickResult &pick)
{
    return pick.detail.getPath();
}

static SbBool
rt_pick_result_path_recorded(
    const std::vector<BObolRtPickResult> &results,
    const BObolRtPickResult &candidate)
{
    const SbString candidatePath = rt_pick_result_path(candidate);
    if (candidatePath.getLength() == 0)
	return FALSE;

    for (size_t i = 0; i < results.size(); i++) {
	if (bu_strcmp(rt_pick_result_path(results[i]).getString(),
		   candidatePath.getString()) == 0)
	    return TRUE;
    }

    return FALSE;
}

static void
insert_rt_pick_result(std::vector<BObolRtPickResult> &results,
		      const BObolRtPickResult &pick,
		      SbBool pickAll)
{
    if (pickAll) {
	std::vector<BObolRtPickResult>::iterator it = results.begin();
	while (it != results.end() && it->distance <= pick.distance)
	    ++it;
	results.insert(it, pick);
	return;
    }

    if (results.empty() || pick.distance < results[0].distance)
	results.assign(1, pick);
}

static const char *
controller_path_skip_slashes(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static const char *
controller_path_leaf(const char *path)
{
    const char *clean = controller_path_skip_slashes(path);
    const char *slash = strrchr(clean, '/');
    return slash && slash[1] ? slash + 1 : clean;
}

static SbBool
controller_path_contains_request(const char *sourcePath,
				 const char *requestPath)
{
    const char *source = controller_path_skip_slashes(sourcePath);
    const char *request = controller_path_skip_slashes(requestPath);
    if (!source[0] || !request[0])
	return FALSE;

    size_t sourceLen = strlen(source);
    if (bu_strncmp(request, source, sourceLen) != 0)
	return FALSE;

    return request[sourceLen] == '\0' || request[sourceLen] == '/' ?
	   TRUE : FALSE;
}

static SbBool
controller_source_matches_request(SoBRLDatabaseSource *source,
				  const BObolSourceMeshRequest &request)
{
    if (!source)
	return FALSE;

    const char *sourcePath = source->path.getValue().getString();
    if (controller_path_contains_request(sourcePath,
					 request.path.getString()))
	return TRUE;

    if (request.path.getLength() > 0)
	return FALSE;

    const char *sourceName = request.sourceName.getString();
    if (!sourceName || !sourceName[0])
	return FALSE;

    const char *cleanSourcePath = controller_path_skip_slashes(sourcePath);
    return bu_strcmp(sourceName, cleanSourcePath) == 0 ||
	   bu_strcmp(sourceName, controller_path_leaf(sourcePath)) == 0 ? TRUE :
	   FALSE;
}

static SoBRLDatabaseSource *
controller_database_source_for_request(BObolViewController *controller,
				       const BObolSourceMeshRequest &request)
{
    SoBRLDatabaseSource *singleSource = NULL;
    SoBRLDatabaseSource *matchedSource = NULL;
    int sourceCountWithDatabase = 0;
    int matchedSourceCount = 0;

    if (!controller)
	return NULL;

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_sources(controller);
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (!source || !source->getDatabase())
	    continue;

	sourceCountWithDatabase++;
	if (!singleSource)
	    singleSource = source;

	if (controller_source_matches_request(source, request)) {
	    matchedSource = source;
	    matchedSourceCount++;
	}
    }

    if (matchedSourceCount == 1)
	return matchedSource;
    if (matchedSourceCount > 1)
	return NULL;

    return sourceCountWithDatabase == 1 ? singleSource : NULL;
}

static int
controller_database_source_count(BObolViewController *controller)
{
    int count = 0;

    if (!controller)
	return 0;

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_sources(controller);
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (source && source->getDatabase())
	    count++;
    }

    return count;
}

static void
controller_source_request_template(BObolLodRequest &requestTemplate,
				   SoBRLDatabaseSource *source)
{
    requestTemplate.clear();
    if (!source)
	return;

    requestTemplate.databaseId =
	controller_database_id(source->getDatabase());
    requestTemplate.sourceRevision = source->sourceRevision.getValue();
}

static SbBool
controller_consume_one_source_full_detail_result(
    SoBRLExportAction &action,
    const BObolSourceMeshRequest &sourceRequest,
    const BObolLodResult &result)
{
    return action.appendSourceBackedFullDetailResult(sourceRequest, result);
}

static SbBool
controller_consume_one_source_full_detail_result(
    SoBRLMeasureAction &action,
    const BObolSourceMeshRequest &sourceRequest,
    const BObolLodResult &result)
{
    return action.consumeSourceBackedFullDetailResult(sourceRequest, result);
}

static SbBool
controller_consume_one_source_full_detail_result(
    SoBRLSnapAction &action,
    const BObolSourceMeshRequest &sourceRequest,
    const BObolLodResult &result)
{
    return action.consumeSourceBackedFullDetailResult(sourceRequest, result);
}

template <typename Action>
static int
controller_consume_source_full_detail(BObolViewController *controller,
				      Action &action,
				      uint64_t generation,
				      int *submittedRequestCount)
{
    if (submittedRequestCount)
	*submittedRequestCount = 0;

    if (!controller)
	return 0;

    const int requestCount =
	action.getSourceBackedFullDetailRequestCount();
    if (requestCount <= 0)
	return 0;

    BObolLodService *service = controller->getLodService();
    if (!service || !service->isRunning())
	return 0;

    std::vector<BObolLodRequest> expectedRequests;
    std::vector<BObolSourceMeshRequest> expectedSourceRequests;
    std::vector<int> expectedRequestIndices;
    std::vector<BObolSourceMeshRequest> submitSourceRequests;
    std::vector<SoBRLDatabaseSource *> requestSources;
    std::vector<int> submitRequestIndices;
    const int databaseSourceCount = controller_database_source_count(controller);
    for (int i = 0; i < requestCount; i++) {
	const BObolSourceMeshRequest &sourceRequest =
	    action.getSourceBackedFullDetailRequest(i);
	SoBRLDatabaseSource *source =
	    controller_database_source_for_request(controller, sourceRequest);
	if (source) {
	    BObolLodRequest requestTemplate;
	    controller_source_request_template(requestTemplate, source);

	    BObolLodRequest request;
	    if (action.makeSourceBackedFullDetailLodRequest(i, request,
		    &requestTemplate)) {
		expectedRequests.push_back(request);
		expectedSourceRequests.push_back(sourceRequest);
		expectedRequestIndices.push_back(i);
		submitSourceRequests.push_back(sourceRequest);
		requestSources.push_back(source);
		submitRequestIndices.push_back(i);
	    }
	}

	if (databaseSourceCount <= 1) {
	    BObolLodRequest request;
	    if (action.makeSourceBackedFullDetailLodRequest(i, request)) {
		expectedRequests.push_back(request);
		expectedSourceRequests.push_back(sourceRequest);
		expectedRequestIndices.push_back(i);
	    }
	}
    }
    if (expectedRequests.empty() && submitSourceRequests.empty())
	return 0;

    std::vector<SbBool> requestMatched(
	static_cast<size_t>(requestCount), FALSE);
    int consumed = 0;
    if (!expectedRequests.empty()) {
	std::vector<BObolLodResult> sourceResults;
	service->drainMatchingResults(sourceResults, expectedRequests);
	if (!sourceResults.empty()) {
	    std::vector<SbBool> used(sourceResults.size(), FALSE);
	    std::vector<SbBool> requestConsumed(
		static_cast<size_t>(requestCount), FALSE);

	    for (size_t i = 0; i < expectedRequests.size(); i++) {
		const int requestIndex = expectedRequestIndices[i];
		if (requestIndex < 0 ||
		    static_cast<size_t>(requestIndex) >=
		    requestConsumed.size() ||
		    requestConsumed[static_cast<size_t>(requestIndex)])
		    continue;

		for (size_t j = 0; j < sourceResults.size(); j++) {
		    if (used[j] ||
			!bobol_lod_result_matches_request(
			    sourceResults[j], expectedRequests[i]))
			continue;

		    used[j] = TRUE;
		    requestMatched[static_cast<size_t>(requestIndex)] =
			TRUE;
		    if (controller_consume_one_source_full_detail_result(action,
			    expectedSourceRequests[i], sourceResults[j])) {
			requestConsumed[static_cast<size_t>(requestIndex)] =
			    TRUE;
			consumed++;
		    }
		    break;
		}
	    }
	}
    }

    int submitted = 0;
    for (size_t i = 0; i < submitSourceRequests.size(); i++) {
	const int requestIndex = submitRequestIndices[i];
	if (requestIndex >= 0 &&
	    static_cast<size_t>(requestIndex) < requestMatched.size() &&
	    requestMatched[static_cast<size_t>(requestIndex)])
	    continue;

	BObolLodRequest requestTemplate;
	controller_source_request_template(requestTemplate, requestSources[i]);
	if (bobol_lod_submit_rt_source_full_detail_request(service,
		generation, submitSourceRequests[i],
		requestSources[i]->getDatabase(), &requestTemplate,
		controller->getMaxExactFullDetailFaceCount(),
		controller->getMaxExactFullDetailPointCount()) != 0)
	    submitted++;
    }
    if (submittedRequestCount)
	*submittedRequestCount = submitted;

    return consumed;
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

/*
 * All normal views in a process share resident mesh assets, worker threads,
 * cache writes, and memory governors.  Each controller still owns an isolated
 * generation and resident-demand consumer.  The weak broker lets the service
 * shut down naturally after the last view releases it.
 */
static std::shared_ptr<BObolLodService>
controller_acquire_managed_lod_service(size_t workerCount)
{
    struct ManagedServiceBroker {
	std::mutex mutex;
	std::weak_ptr<BObolLodService> service;
    };
    static ManagedServiceBroker broker;

    std::lock_guard<std::mutex> lock(broker.mutex);
    std::shared_ptr<BObolLodService> service = broker.service.lock();
    if (!service) {
	service = std::make_shared<BObolLodService>();
	if (!service->start(workerCount, TRUE))
	    return std::shared_ptr<BObolLodService>();
	broker.service = service;
    } else if (!service->isRunning()) {
	if (!service->start(workerCount, TRUE))
	    return std::shared_ptr<BObolLodService>();
    } else if (!service->ensureWorkerCount(workerCount)) {
	return std::shared_ptr<BObolLodService>();
    }
    return service;
}

BObolViewController::BObolViewController(void) :
    d(new Impl(this))
{
    this->d->viewAttachment->ref();
    SoBRLCadRenderBatch *batch = new SoBRLCadRenderBatch;
    batch->ref();
    batch->addChild(this->d->viewport->getRoot());
    this->d->renderBatchRoot = batch;
    SoSeparator *presentation = new SoSeparator;
    presentation->ref();
    this->d->framebufferUnderlayRoot = new SoGroup;
    /* Unlike underlay/overlay this root is parented by GED's per-view render
     * composition, so the controller keeps a reference across rebinds. */
    this->d->framebufferInterlayRoot = new SoGroup;
    this->d->framebufferInterlayRoot->ref();
    this->d->framebufferOverlayRoot = new SoGroup;
    presentation->addChild(this->d->framebufferUnderlayRoot);
    presentation->addChild(batch);
    presentation->addChild(this->d->framebufferOverlayRoot);
    this->d->renderPresentationRoot = presentation;
    controller_initialize_render_action(this->d->renderManager);
    controller_configure_render_environment(this->d->viewport);
}

BObolViewController::BObolViewController(SoNode *root, SoCamera *camera) :
    d(new Impl(this))
{
    this->d->viewAttachment->ref();
    SoBRLCadRenderBatch *batch = new SoBRLCadRenderBatch;
    batch->ref();
    batch->addChild(this->d->viewport->getRoot());
    this->d->renderBatchRoot = batch;
    SoSeparator *presentation = new SoSeparator;
    presentation->ref();
    this->d->framebufferUnderlayRoot = new SoGroup;
    /* Unlike underlay/overlay this root is parented by GED's per-view render
     * composition, so the controller keeps a reference across rebinds. */
    this->d->framebufferInterlayRoot = new SoGroup;
    this->d->framebufferInterlayRoot->ref();
    this->d->framebufferOverlayRoot = new SoGroup;
    presentation->addChild(this->d->framebufferUnderlayRoot);
    presentation->addChild(batch);
    presentation->addChild(this->d->framebufferOverlayRoot);
    this->d->renderPresentationRoot = presentation;
    controller_initialize_render_action(this->d->renderManager);
    controller_configure_render_environment(this->d->viewport);
    this->setSceneRoot(root);
    this->setCamera(camera);
}

BObolViewController::~BObolViewController(void)
{
	this->setPresentationSyncCallback(NULL, NULL);
    ControllerFrameRequestState *frameRequestState = NULL;
    {
	std::lock_guard<std::mutex> lock(this->d->frameRequestMutex);
	frameRequestState = static_cast<ControllerFrameRequestState *>(
	    this->d->frameRequestUserData);
	this->d->frameRequestCallback = NULL;
	this->d->frameRequestUserData = NULL;
    }
    if (frameRequestState) {
	frameRequestState->close();
	delete frameRequestState;
    }
    delete this->d->imageRenderer;
    this->d->imageRenderer = NULL;
    this->d->imageRendererManager = NULL;
    this->clearProgressiveProviders();
    this->setLodService(NULL);
    this->d->managedLodService.reset();
    this->d->managedLodWorkerCount = 0;
    this->clearRtPickCaches();
    delete this->d->featureStore;
    this->d->featureStore = NULL;
    delete this->d->polygonStore;
    this->d->polygonStore = NULL;
    delete this->d->selectionStore;
    this->d->selectionStore = NULL;
    this->setCamera(NULL);
    this->setSceneRoot(NULL);
    this->d->viewAttachment->unref();
    this->d->viewAttachment = NULL;
    this->d->renderManager->setSceneGraph(NULL);
    delete this->d->renderManager;
    this->d->renderManager = NULL;
    if (this->d->renderPresentationRoot) {
	this->d->renderPresentationRoot->unref();
	this->d->renderPresentationRoot = NULL;
    }
    this->d->framebufferUnderlayRoot = NULL;
    if (this->d->framebufferInterlayRoot) {
	this->d->framebufferInterlayRoot->unref();
	this->d->framebufferInterlayRoot = NULL;
    }
    this->d->framebufferOverlayRoot = NULL;
    if (this->d->renderBatchRoot) {
	this->d->renderBatchRoot->unref();
	this->d->renderBatchRoot = NULL;
    }
    delete this->d->viewport;
    this->d->viewport = NULL;
}

void
BObolViewController::setViewportSceneGraphWithLod(SoNode *root)
{
    SoBRLCadRenderBatch *batch =
	dynamic_cast<SoBRLCadRenderBatch *>(this->d->renderBatchRoot);
    if (batch)
	batch->setBatchSourceRoot(NULL);
    if (batch)
	batch->setSoftwareWireMode(this->d->softwareWireMode);
    if (this->d->renderLodRoot) {
	this->d->viewport->setSceneGraph(NULL);
	this->d->renderLodRoot->unref();
	this->d->renderLodRoot = NULL;
    }

    if (!root) {
	this->d->viewport->setSceneGraph(NULL);
	return;
    }

    SoBRLViewLodGroup *wrapper = new SoBRLViewLodGroup;
    wrapper->ref();
    wrapper->setViewLodState(this->d->viewAttachment->getViewLodState());
    wrapper->setSoftwareWireMode(this->d->softwareWireMode);
    wrapper->addChild(root);
    this->d->renderLodRoot = wrapper;
    this->d->viewport->setSceneGraph(wrapper);
    if (batch)
	batch->setBatchSourceRoot(root);
}

void
BObolViewController::setSceneRoot(SoNode *root)
{
    this->cancelActiveLodGeneration();
    this->clearRtPickCaches();
    this->d->viewAttachment->setSceneRoot(root);
    this->d->sceneController.setSceneRoot(root);
    if (root && this->d->framebufferInterlayRoot) {
	/* A controller used without GED has no separate local feature root.
	 * Keep interlay visible after its scene; hosted GED replaces this with
	 * the more precise shared/interlay/local composition. */
	SoSeparator *renderRoot = new SoSeparator;
	renderRoot->addChild(root);
	renderRoot->addChild(this->d->framebufferInterlayRoot);
	this->setViewportSceneGraphWithLod(renderRoot);
    } else {
	this->setViewportSceneGraphWithLod(root);
    }
    this->syncRenderManager();
    this->requestRender("scene-root");
}

SoNode *
BObolViewController::getSceneRoot(void) const
{
    return this->d->viewAttachment->getSceneRoot();
}

void
BObolViewController::setRenderSceneRoot(SoNode *root)
{
    this->cancelActiveLodGeneration();
    this->clearRtPickCaches();
    this->d->viewAttachment->clearViewLodState();
    this->setViewportSceneGraphWithLod(root);
    this->syncRenderManager();
    this->requestRender("render-scene-root");
}

SoNode *
BObolViewController::getRenderSceneRoot(void) const
{
    return this->d->viewport->getSceneGraph();
}

SoNode *
BObolViewController::getRenderRoot(void) const
{
	if (!this->d->endpointGraphicalRenderingEnabled.load(
		std::memory_order_acquire))
	    return NULL;
	return this->d->renderPresentationRoot ? this->d->renderPresentationRoot :
	this->d->renderBatchRoot ? this->d->renderBatchRoot :
	this->d->viewport->getRoot();
}

SoGroup *
BObolViewController::getFramebufferUnderlayRoot(void) const
{
    return this->d->framebufferUnderlayRoot;
}

SoGroup *
BObolViewController::getFramebufferInterlayRoot(void) const
{
    return this->d->framebufferInterlayRoot;
}

SoGroup *
BObolViewController::getFramebufferOverlayRoot(void) const
{
    return this->d->framebufferOverlayRoot;
}

void
BObolViewController::setViewAttachment(BObolViewAttachment *attachment)
{
    if (!attachment || attachment == this->d->viewAttachment)
	return;

    this->cancelActiveLodGeneration();

    SoNode *root = this->getSceneRoot();
    if (root)
	root->ref();

    SoNode *renderScene = NULL;
    if (this->d->renderLodRoot && this->d->renderLodRoot->getNumChildren() > 0)
	renderScene = this->d->renderLodRoot->getChild(0);
    else
	renderScene = this->d->viewport->getSceneGraph();
    if (renderScene)
	renderScene->ref();

    attachment->ref();
    this->d->viewAttachment->unref();
    this->d->viewAttachment = attachment;

    if (root && !this->d->viewAttachment->hasSceneRoot())
	this->d->viewAttachment->setSceneRoot(root);
    if (root)
	root->unref();

    this->setViewportSceneGraphWithLod(renderScene);
    if (renderScene)
	renderScene->unref();
    this->clearRtPickCaches();
    this->syncRenderManager();
    this->requestRender("view-attachment");
}

BObolViewAttachment *
BObolViewController::getViewAttachment(void) const
{
    return this->d->viewAttachment;
}

BObolViewLodState *
BObolViewController::getViewLodState(void) const
{
    return this->d->viewAttachment->getViewLodState();
}

void
BObolViewController::clearViewLodState(void)
{
    this->cancelActiveLodGeneration();
    this->d->viewAttachment->clearViewLodState();
}

void
BObolViewController::renderBackground(void) const
{
    SoDB::ContextManager *contextManager = this->getRenderContextManager();
    if (!contextManager)
	return;
    /* Wrapper managers such as the Qt and Tk providers can dispatch to live
     * system GL or an offscreen fallback on the same thread.  Resolve the
     * function table for the currently active context instead of caching it
     * solely by wrapper address. */
    ControllerGLFunctions functions;
    functions.load(contextManager);
    if (!functions.complete())
	return;
    ControllerGLFunctions *gl = &functions;
    const SbColor &bottom = this->d->backgroundBottom;
    const SbColor &top = this->d->backgroundTop;
    if (bottom == top) {
	gl->clearColor(bottom[0], bottom[1], bottom[2], 1.0f);
	gl->clear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	return;
    }

    GLint matrixMode = GL_MODELVIEW;
    gl->getIntegerv(GL_MATRIX_MODE, &matrixMode);
    gl->pushAttrib(GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
	GL_CURRENT_BIT);
    /* A gradient quad is not a substitute for clearing the render target.
     * Rasterization coverage at the viewport boundary is implementation
     * dependent enough that the top scanline can retain geometry from the
     * previous QOpenGLWidget frame.  Clear the whole draw buffer first, and
     * do not inherit a scissor left by either Coin traversal or the Qt
     * compositor.  The pushed enable state restores the caller's scissor
     * policy after the background pass. */
    gl->disable(GL_SCISSOR_TEST);
    gl->clearColor(bottom[0], bottom[1], bottom[2], 1.0f);
    gl->clear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    gl->disable(GL_LIGHTING);
    gl->disable(GL_DEPTH_TEST);
    gl->depthMask(GL_FALSE);
    gl->matrixMode(GL_PROJECTION);
    gl->pushMatrix();
    gl->loadIdentity();
    gl->matrixMode(GL_MODELVIEW);
    gl->pushMatrix();
    gl->loadIdentity();
    gl->begin(GL_QUADS);
    gl->color3f(bottom[0], bottom[1], bottom[2]);
    gl->vertex2f(-1.0f, -1.0f);
    gl->vertex2f(1.0f, -1.0f);
    gl->color3f(top[0], top[1], top[2]);
    gl->vertex2f(1.0f, 1.0f);
    gl->vertex2f(-1.0f, 1.0f);
    gl->end();
    gl->matrixMode(GL_MODELVIEW);
    gl->popMatrix();
    gl->matrixMode(GL_PROJECTION);
    gl->popMatrix();
    gl->matrixMode(matrixMode);
    gl->popAttrib();
}

SbBool
BObolViewController::syncCameraFromViewContext(const void *viewCtx,
	SbBool createCamera, SbBool *changedOut)
{
    if (changedOut)
	*changedOut = FALSE;
    const struct bv *view =
	bv_context_view_const((const struct bv_context *)viewCtx);
    if (!view)
	return FALSE;

    const double perspectiveDegrees =
	bv_perspective_get(view);
    const SbBool wantPerspective = perspectiveDegrees > SMALL_FASTF ?
				   TRUE : FALSE;

    SbBool cameraReplaced = FALSE;
    SoCamera *camera = this->d->activeCamera;
    if (!camera ||
	(wantPerspective &&
	 !camera->isOfType(SoPerspectiveCamera::getClassTypeId())) ||
	(!wantPerspective &&
	 !camera->isOfType(SoOrthographicCamera::getClassTypeId()))) {
	if (!createCamera)
	    return FALSE;

	camera = wantPerspective ?
		 static_cast<SoCamera *>(new SoPerspectiveCamera) :
		 static_cast<SoCamera *>(new SoOrthographicCamera);
	this->setCamera(camera);
	cameraReplaced = TRUE;
    }

    int viewWidth = bv_width_get(view);
    int viewHeight = bv_height_get(view);
    SbVec2s window = this->d->viewportRegion.getWindowSize();
    if (window[0] <= 1 && window[1] <= 1 &&
	viewWidth > 0 && viewHeight > 0) {
	this->setViewportSize(static_cast<unsigned int>(viewWidth),
			      static_cast<unsigned int>(viewHeight));
	window = this->d->viewportRegion.getWindowSize();
    }

    double aspect = controller_aspect_from_region(this->d->viewportRegion);
    if (aspect <= SMALL_FASTF && viewWidth > 0 && viewHeight > 0)
	aspect = static_cast<double>(viewWidth) /
		 static_cast<double>(viewHeight);
    if (aspect <= SMALL_FASTF)
	aspect = 1.0;

    mat_t viewRotation;
    MAT_IDN(viewRotation);
    (void)bv_rotation_get(viewRotation, view);

    vect_t center;
    VSETALL(center, 0.0);
    (void)bv_center_get(center, view);

    const double horizontalSizeRaw = bv_size_get(view);
    double horizontalSize = horizontalSizeRaw;
    if (horizontalSize <= SMALL_FASTF) {
	const double scale = bv_scale_get(view);
	horizontalSize = scale > SMALL_FASTF ? scale * 2.0 : 2.0;
    }
    const double verticalSize = horizontalSize / aspect;

    double heightAngle = perspectiveDegrees * DEG2RAD;
    if (heightAngle <= SMALL_FASTF)
	heightAngle = 2.0 * std::atan(0.1);
    if (heightAngle < 0.001)
	heightAngle = 0.001;
    if (heightAngle > 3.0)
	heightAngle = 3.0;

    const SbBool orthographic =
	camera->isOfType(SoOrthographicCamera::getClassTypeId());
    double distance = orthographic ?
		      horizontalSize * 5.0 :
		      (verticalSize * 0.5) / std::tan(heightAngle * 0.5);
    if (distance <= horizontalSize * 0.5)
	distance = horizontalSize * 0.5;
    if (distance <= SMALL_FASTF)
	distance = 1.0;

    const SbRotation orientation =
	bobol_camera_orientation_from_brl_rotation(viewRotation);
    const double viewZ[3] = {
	viewRotation[8], viewRotation[9], viewRotation[10]
    };

    /* Track the camera-driven headlight: rewrite its world-space direction from
     * the current camera orientation so it follows the viewer.  If the camera
     * orientation changed, the camera->orientation assignment below sets
     * `changed` and drives the render; when idle the direction is unchanged so
     * this does not churn. */
    this->d->lastCameraOrientation = orientation;
    this->applyTrackedHeadlight();

    SoClipPlane *clipMinimumNode =
	controller_clip_plane(this->d->viewport, TRUE);
    SoClipPlane *clipMaximumNode =
	controller_clip_plane(this->d->viewport, FALSE);
    const SbBool clipping = bv_zclip_get(view) ? TRUE : FALSE;
    if (clipMinimumNode && clipMaximumNode) {
	const double viewScale = horizontalSize * 0.5;
	const SbVec3f minimumNormal(
	    static_cast<float>(viewZ[X]),
	    static_cast<float>(viewZ[Y]),
	    static_cast<float>(viewZ[Z]));
	const SbVec3f maximumNormal = -minimumNormal;
	const SbVec3f minimumPoint(
	    static_cast<float>(center[X] + viewZ[X] *
		this->d->clipMinimum * viewScale),
	    static_cast<float>(center[Y] + viewZ[Y] *
		this->d->clipMinimum * viewScale),
	    static_cast<float>(center[Z] + viewZ[Z] *
		this->d->clipMinimum * viewScale));
	const SbVec3f maximumPoint(
	    static_cast<float>(center[X] + viewZ[X] *
		this->d->clipMaximum * viewScale),
	    static_cast<float>(center[Y] + viewZ[Y] *
		this->d->clipMaximum * viewScale),
	    static_cast<float>(center[Z] + viewZ[Z] *
		this->d->clipMaximum * viewScale));
	clipMinimumNode->plane = SbPlane(minimumNormal, minimumPoint);
	clipMaximumNode->plane = SbPlane(maximumNormal, maximumPoint);
	clipMinimumNode->on = clipping;
	clipMaximumNode->on = clipping;
    }

    const float desiredAspect = static_cast<float>(aspect);
    const SbVec3f desiredPosition(
	static_cast<float>(center[X] + viewZ[X] * distance),
	static_cast<float>(center[Y] + viewZ[Y] * distance),
	static_cast<float>(center[Z] + viewZ[Z] * distance));
    const float desiredFocal = static_cast<float>(distance);
    const float desiredNear = static_cast<float>(
	std::max(distance * 0.001, 1.0e-6));
    const float desiredFar = static_cast<float>(
	std::max(distance + horizontalSize * 100.0, distance + 1.0));
    const auto float_changed = [](float current, float desired) {
	const float scale = std::max(1.0f, std::fabs(desired));
	return std::fabs(current - desired) > 1.0e-6f * scale;
    };

    SbBool changed = cameraReplaced;
    if (camera->viewportMapping.getValue() != SoCamera::LEAVE_ALONE) {
	camera->viewportMapping = SoCamera::LEAVE_ALONE;
	changed = TRUE;
    }
    if (float_changed(camera->aspectRatio.getValue(), desiredAspect)) {
	camera->aspectRatio = desiredAspect;
	changed = TRUE;
    }
    if (camera->position.getValue() != desiredPosition) {
	camera->position = desiredPosition;
	changed = TRUE;
    }
    if (camera->orientation.getValue() != orientation) {
	camera->orientation = orientation;
	changed = TRUE;
    }
    if (float_changed(camera->focalDistance.getValue(), desiredFocal)) {
	camera->focalDistance = desiredFocal;
	changed = TRUE;
    }
    if (float_changed(camera->nearDistance.getValue(), desiredNear)) {
	camera->nearDistance = desiredNear;
	changed = TRUE;
    }
    if (float_changed(camera->farDistance.getValue(), desiredFar)) {
	camera->farDistance = desiredFar;
	changed = TRUE;
    }

    if (camera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
	SoPerspectiveCamera *perspectiveCamera =
	    static_cast<SoPerspectiveCamera *>(camera);
	const float desiredAngle = static_cast<float>(heightAngle);
	if (float_changed(perspectiveCamera->heightAngle.getValue(),
		desiredAngle)) {
	    perspectiveCamera->heightAngle = desiredAngle;
	    changed = TRUE;
	}
    } else if (orthographic) {
	SoOrthographicCamera *orthographicCamera =
	    static_cast<SoOrthographicCamera *>(camera);
	const float desiredHeight = static_cast<float>(verticalSize);
	if (float_changed(orthographicCamera->height.getValue(),
		desiredHeight)) {
	    orthographicCamera->height = desiredHeight;
	    changed = TRUE;
	}
    }

    if (changed) {
	this->syncLodViewSignature(TRUE);
	this->requestRender("rt-view-camera");
    }
    if (changedOut)
	*changedOut = changed;
    return TRUE;
}

SbBool
BObolViewController::getViewInfo(struct bv_view_info *info) const
{
    if (!info)
	return FALSE;

    bv_view_info_init(info);
    if (this->d->viewAttachment) {
	struct bv_lod_policy policy;
	bv_lod_policy_init(&policy);
	this->d->viewAttachment->getLodPolicy(&policy);
	info->lod.scale = policy.scale;
	info->lod.curve_scale = policy.curve_scale;
	info->lod.point_scale = policy.point_scale;
	info->lod.bot_threshold = policy.bot_threshold;
    }

    SbVec2s window = this->d->viewportRegion.getWindowSize();
    info->width = window[0] > 0 ? window[0] : 1;
    info->height = window[1] > 0 ? window[1] : 1;

    if (!this->d->activeCamera) {
	bv_view_info_sanitize(info);
	return FALSE;
    }

    if (this->d->activeCamera->isOfType(SoOrthographicCamera::getClassTypeId())) {
	SoOrthographicCamera *camera =
	    static_cast<SoOrthographicCamera *>(this->d->activeCamera);
	double aspect = controller_aspect_from_region(this->d->viewportRegion);
	if (aspect <= SMALL_FASTF)
	    aspect = camera->aspectRatio.getValue();
	if (aspect <= SMALL_FASTF)
	    aspect = 1.0;
	info->size = camera->height.getValue() * aspect;
    } else if (this->d->activeCamera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
	SoPerspectiveCamera *camera =
	    static_cast<SoPerspectiveCamera *>(this->d->activeCamera);
	double focal = this->d->activeCamera->focalDistance.getValue();
	double angle = camera->heightAngle.getValue();
	double aspect = controller_aspect_from_region(this->d->viewportRegion);
	if (aspect <= SMALL_FASTF)
	    aspect = this->d->activeCamera->aspectRatio.getValue();
	if (aspect <= SMALL_FASTF)
	    aspect = 1.0;
	if (focal <= 0.0)
	    focal = 1.0;
	if (angle <= 0.0)
	    angle = 2.0 * std::atan(0.1);
	info->size = 2.0 * focal * std::tan(angle * 0.5) * aspect;
    } else {
	info->size = this->d->activeCamera->focalDistance.getValue();
    }

    bv_view_info_sanitize(info);
    return TRUE;
}

SbBool
BObolViewController::realizePending(void)
{
    const SbBool ret = this->d->sceneController.realizePending();
    this->requestRender(ret ? "realize" : "realize-failed");
    return ret;
}

SbBool
BObolViewController::isForceRealizeDisplay(void) const
{
    /* No attachment / policy yet -> keep the classic force-realize behavior so
     * nothing renders half-formed before a policy is established. */
    if (!this->d->viewAttachment)
	return TRUE;
    struct bv_lod_policy policy;
    bv_lod_policy_init(&policy);
    this->d->viewAttachment->getLodPolicy(&policy);
    return (policy.policy == BV_LOD_OFF ||
	    (!policy.mesh_enabled && !policy.csg_enabled)) ? TRUE : FALSE;
}

unsigned int
BObolViewController::getLastVisitedSourceCount(void) const
{
    return this->d->sceneController.getLastVisitedSourceCount();
}

unsigned int
BObolViewController::getLastRealizedSourceCount(void) const
{
    return this->d->sceneController.getLastRealizedSourceCount();
}

unsigned int
BObolViewController::getLastFailedSourceCount(void) const
{
    return this->d->sceneController.getLastFailedSourceCount();
}

const SbString &
BObolViewController::getLastDiagnostics(void) const
{
    return this->d->sceneController.getLastDiagnostics();
}

void
BObolViewController::setEndpointGraphicalRenderingEnabled(SbBool enabled)
{
    const int requested = enabled ? 1 : 0;
    if (this->d->endpointGraphicalRenderingEnabled.exchange(requested,
	    std::memory_order_acq_rel) == requested)
	return;
    this->syncRenderManager();
    if (requested)
	this->requestRender("render-engine");
    else
	this->clearRenderRequest();
}

void
BObolViewController::invalidateRendererPerformanceHistory(void)
{
    /* Renderer capacity is part of an exact-view quality proof even though it
	* is deliberately not part of camera identity.  Retain immutable scene
	* residency and the coherent currently presented cut, but invalidate every
	* timing-derived proof and arrange one measured-frame rescan. */
    this->d->resetRendererPerformanceEvidence();
    if (this->d->lodAutoSubmit)
	this->markProgressiveWorkPending();
    this->requestRender("renderer-performance");
}

void
BObolViewController::requestRender(const char *reason)
{
    SbBool wakeEndpoint = FALSE;
    uint64_t requestSerial = 0;
    {
	std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
	wakeEndpoint = this->d->renderRequested ? FALSE : TRUE;
	this->d->renderRequested = TRUE;
	this->d->renderLodCapacityRelevant = TRUE;
	this->d->renderReason = reason ? reason : "";
	if (++this->d->hostWorkRevision == 0)
	    ++this->d->hostWorkRevision;
	requestSerial = ++this->d->renderRequestSerial;
	if (requestSerial == 0)
	    requestSerial = ++this->d->renderRequestSerial;
    }
    /*
     * A render request is itself the complete host contract.  Requiring every
     * caller to remember a separate notifyFrameRequest() made camera sync and
     * final calibration replays easy to strand after the progressive timer
     * became idle.  Wake only on the empty-to-pending edge; repeated reasons
     * still update diagnostics without flooding the Qt event queue.
     */
    if (wakeEndpoint && controller_lod_trace_enabled(
	    "BOBOL_LOD_TRACE_RENDER_REQUEST",
	    this->d->lodViewRevision.value()))
	bu_log("BObol LoD render request serial=%llu view=%llu "
	       "policy=%llu reason=%s progressive=%d\n",
	       static_cast<unsigned long long>(requestSerial),
	       static_cast<unsigned long long>(
		   this->d->lodViewRevision.value()),
	       static_cast<unsigned long long>(
		   this->d->lodPolicyRevision.value()),
	       reason ? reason : "",
	       this->hasProgressiveWorkPending() ? 1 : 0);
    if (wakeEndpoint)
	this->notifyFrameRequest(reason);
}

void
BObolViewController::requestPresentationRender(const char *reason)
{
    SbBool wakeEndpoint = FALSE;
    {
	std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
	wakeEndpoint = this->d->renderRequested ? FALSE : TRUE;
	if (!this->d->renderRequested)
	    this->d->renderLodCapacityRelevant = FALSE;
	this->d->renderRequested = TRUE;
	this->d->renderReason = reason ? reason : "";
	if (++this->d->hostWorkRevision == 0)
	    ++this->d->hostWorkRevision;
	if (++this->d->renderRequestSerial == 0)
	    ++this->d->renderRequestSerial;
    }
    if (wakeEndpoint)
	this->notifyFrameRequest(reason);
}

void
BObolViewController::setFrameRequestCallback(
    BObolFrameRequestCallback callback, void *userData)
{
    ControllerFrameRequestState *replacement = callback ?
	new (std::nothrow) ControllerFrameRequestState(callback, userData) : NULL;
    if (callback && !replacement)
	return;

    ControllerFrameRequestState *previous = NULL;
    {
	std::lock_guard<std::mutex> lock(this->d->frameRequestMutex);
	previous = static_cast<ControllerFrameRequestState *>(
	    this->d->frameRequestUserData);
	this->d->frameRequestCallback = callback;
	this->d->frameRequestUserData = replacement;
    }
    if (previous) {
	previous->close();
	delete previous;
    }

    /*
     * Installing a host is a level-triggered attachment, not an edge-only
     * subscription.  Providers and draw commands may legitimately request a
     * frame before Qt (or another graphical host) finishes binding its
     * callback.  Replaying the already-pending level here prevents that work
     * from remaining invisible until an unrelated expose, input event, or
     * explicit synchronous pump happens to arrive.
     *
     * Do this after the old callback has quiesced and outside the state mutex:
     * notifyFrameRequest() owns the normal dispatch/lifetime protocol and a
     * host callback may immediately re-enter the controller.
     */
    if (callback && this->getHostWorkSnapshot().flags !=
	BOBOL_HOST_WORK_NONE)
	this->notifyFrameRequest("host-attached-pending");
}

void
BObolViewController::clearFrameRequestCallback(void *userData)
{
    ControllerFrameRequestState *state = NULL;
    {
	std::lock_guard<std::mutex> lock(this->d->frameRequestMutex);
	state = static_cast<ControllerFrameRequestState *>(
	    this->d->frameRequestUserData);
	if (!state || state->userData != userData)
	    return;
	this->d->frameRequestCallback = NULL;
	this->d->frameRequestUserData = NULL;
    }
    state->close();
    delete state;
}

void
BObolViewController::setPresentationSyncCallback(
    BObolPresentationSyncCallback callback, void *userData)
{
    std::lock_guard<std::mutex> lock(this->d->presentationSyncMutex);
    this->d->presentationSyncCallback = callback;
    this->d->presentationSyncUserData = callback ? userData : NULL;
}

void
BObolViewController::clearPresentationSyncCallback(void *userData)
{
    std::lock_guard<std::mutex> lock(this->d->presentationSyncMutex);
    if (this->d->presentationSyncUserData != userData)
	return;
    this->d->presentationSyncCallback = NULL;
    this->d->presentationSyncUserData = NULL;
}

void
BObolViewController::synchronizePresentation(void)
{
    const uint64_t started = this->beginRenderTiming();
    BObolPresentationSyncCallback callback = NULL;
    void *userData = NULL;
    {
	std::lock_guard<std::mutex> lock(this->d->presentationSyncMutex);
	callback = this->d->presentationSyncCallback;
	userData = this->d->presentationSyncUserData;
    }
    if (callback)
	(*callback)(userData);
    controller_synchronize_compact_cad_presentations(this);
    const uint64_t completed = this->beginRenderTiming();
    this->d->lastPresentationSyncTimeNanoseconds =
	(completed > started) ? completed - started : 0;
}

void
BObolViewController::notifyFrameRequest(const char *reason)
{
    if (!this->d->endpointGraphicalRenderingEnabled.load(
	    std::memory_order_acquire))
	return;
    BObolFrameRequestCallback callback = NULL;
    void *userData = NULL;
    ControllerFrameRequestState *state = NULL;
    {
	std::lock_guard<std::mutex> lock(this->d->frameRequestMutex);
	state = static_cast<ControllerFrameRequestState *>(
	    this->d->frameRequestUserData);
	if (this->d->frameRequestCallback && state && state->beginDispatch()) {
	    callback = state->callback;
	    userData = state->userData;
	}
    }
    if (!callback)
	return;

    (*callback)(userData, reason ? reason : "");
    state->finishDispatch();
}

BObolHostWorkSnapshot
BObolViewController::getHostWorkSnapshot(void) const
{
    BObolHostWorkSnapshot snapshot;
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    snapshot.revision = this->d->hostWorkRevision;
    snapshot.renderRevision = this->d->renderRequestSerial;
    if (this->d->progressiveWorkPending)
	snapshot.flags |= BOBOL_HOST_WORK_PUMP;
    if (this->d->renderRequested) {
	snapshot.flags |= BOBOL_HOST_WORK_RENDER;
	if (this->d->renderLodCapacityRelevant)
	    snapshot.flags |= BOBOL_HOST_WORK_CAPACITY_SAMPLE;
    }
    return snapshot;
}

void
BObolViewController::clearRenderRequest(void)
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    const SbBool changed = this->d->renderRequested ||
	this->d->renderLodCapacityRelevant ||
	this->d->renderReason.getLength() > 0 ? TRUE : FALSE;
    this->d->renderRequested = FALSE;
    this->d->renderLodCapacityRelevant = FALSE;
    this->d->renderReason = "";
    if (changed && ++this->d->hostWorkRevision == 0)
	++this->d->hostWorkRevision;
    if (++this->d->renderRequestSerial == 0)
	++this->d->renderRequestSerial;
}

SbBool
BObolViewController::consumeRenderRequest(SbString *reason,
	SbBool *lodCapacityRelevant)
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    const SbBool ret = this->d->renderRequested;
    if (reason)
	*reason = this->d->renderReason;
    if (lodCapacityRelevant)
	*lodCapacityRelevant = this->d->renderLodCapacityRelevant;
    this->d->renderRequested = FALSE;
    this->d->renderLodCapacityRelevant = FALSE;
    this->d->renderReason = "";
    if (ret && ++this->d->hostWorkRevision == 0)
	++this->d->hostWorkRevision;
    if (++this->d->renderRequestSerial == 0)
	++this->d->renderRequestSerial;
    return ret;
}

void
BObolViewController::clearRenderRequestIfUnchanged(uint64_t serial)
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    if (this->d->renderRequestSerial != serial)
	return;
    this->d->renderRequested = FALSE;
    this->d->renderLodCapacityRelevant = FALSE;
    this->d->renderReason = "";
    if (++this->d->hostWorkRevision == 0)
	++this->d->hostWorkRevision;
    if (++this->d->renderRequestSerial == 0)
	++this->d->renderRequestSerial;
}

uint64_t
BObolViewController::renderRequestSerialGet(void) const
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    return this->d->renderRequestSerial;
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

static const char *
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
    const uint64_t consumerId = static_cast<uint64_t>(
	reinterpret_cast<uintptr_t>(controller));
    const bool generationReady = generation != 0 &&
	service->queuedResultCountForGeneration(generation) > 0;
    const bool compactionReady =
	service->queuedResidentMeshCompactionResultCountForDiagnostics(
	    consumerId) > 0;
    if (!generationReady && !compactionReady)
	return;

    if (generationReady) {
	controller->d->lodResultsPending.store(1);
	int64_t expected = 0;
	(void)controller->d->lodResultsFirstReadyMicroseconds.
	    compare_exchange_strong(expected, bu_gettime());
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
	    static_cast<uint64_t>(reinterpret_cast<uintptr_t>(this)));
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
    this->d->lodResultsPending.store(0);
    this->d->lodResultsFirstReadyMicroseconds.store(0);
    this->d->lodActiveGeneration = 0;
    this->d->lodSubmissionSourceIndex = 0;
    this->d->lodSubmissionEntryOffset = 0;
    this->d->clearLodSubmissionPlan();
    this->d->lodSubmissionPending = FALSE;
    this->d->lodSubmissionRescanPending = FALSE;
    this->d->lodViewDemandPolicy.reset();
    this->d->lodInteractionStartScaleSignatureValid = FALSE;
    this->d->lodInteractionStartedFromReadyView = FALSE;
    this->d->lodStableBudgetBeforeInteraction = 0;
    this->d->lodStableBudgetBeforeInteractionValid = FALSE;
    this->d->lodRetainPoseOccurrenceCuts = FALSE;
    this->d->lodRetainedImportanceCensusPending = FALSE;
    this->d->lodResidentGrowthPolicy.reset();
    this->d->lodResidentGrowthResidencyDrainActive = FALSE;
    this->d->lodResidentAdmissionRetryActive = FALSE;
    this->d->lodDiscretePopulationTrialAvailable = FALSE;
    this->d->lodStaticOverscanActive = FALSE;
    this->d->lodStaticOverscanRejected = FALSE;
    this->d->lodCoveragePolicy.reset();
    if (service)
	this->d->lodCoveragePolicy.activate(true);
    this->d->lodRetainedRefinementPending = FALSE;
    this->d->lodRetainedResidencyPending = FALSE;
    this->d->lodStructuralPresentationRepairPending = FALSE;
    this->d->lodStructuralRepairFrontierCount = 0;
    this->d->lodStructuralCoverageCostReservation = 0;
    this->d->lodRetainedRefinementCutAdvanced = FALSE;
    this->d->lodRetainedRefinementBudgetBlocked = FALSE;
    this->d->lodBudgetPolicy.resetCalibration();
    this->d->lodBudgetPolicy.clearRetainedQualityFloorBudget();
    this->d->lodPassAdmittedWork = FALSE;
    this->d->lodBudgetPolicy.resetPass();
    this->d->lodBudgetPolicy.resetOverloadRecovery();
    this->d->lodSettleAfterRenderSerial = 0;
    this->d->lodLastSubmittedViewRevision.reset();
    this->d->lodLastSubmittedPolicyRevision.reset();
    this->d->lodLastSubmittedSources.clear();
    this->d->lodSubmissionDeltaActive = FALSE;
    this->d->lodSubmissionDeltaSources.clear();
    this->d->lodSubmissionDeltaPlans.clear();
    this->d->lodResourcePolicy.resetForServiceChange();
    this->d->lodHeadroomPolicy.cancelRetry();
    this->d->lodResidentAdmissionRevision =
	service ? service->residentMeshAdmissionRevision() : 0;
    this->d->lodCompactionPolicy.resetForServiceChange(
	service != NULL, bu_gettime(), 750000);
    this->d->lodFrameObligation.reset();
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodPublicationPolicy.reset();
    this->d->lodInterruptedPresentationReplayPending = FALSE;

    if (this->d->lodService)
	this->d->lodResultSubscriberId =
	    this->d->lodService->subscribeResultReady(
		BObolViewController::lodResultReadyCB, this);
    this->d->reconcilePhase(BObolLodStateMachine::Event::SERVICE_CHANGED);
}

uint64_t
BObolViewController::beginLodGeneration(void)
{
    if (!this->d->lodService || !this->d->lodService->isRunning())
	return 0;
    this->cancelActiveLodGeneration();
    this->d->lodActiveGeneration =
	this->d->lodService->beginGeneration();
    this->d->lodResultsPending.store(0);
    this->d->lodResultsFirstReadyMicroseconds.store(0);
    this->d->reconcilePhase(
	BObolLodStateMachine::Event::WORK_SCHEDULED);
    return this->d->lodActiveGeneration;
}

void
BObolViewController::resetDiscoveryPointProxyFloor(SbBool requestFrame)
{
    const float oldEffective = std::max(
	this->d->lodPresentationPointProxyPixelThreshold,
	this->d->lodDiscoveryPointProxyPixelThreshold);
    this->d->lodDiscoveryPointProxyPixelThreshold = 1.0f;
    this->d->lodDiscoveryPointProxyPolicy.reset();
    this->d->lodAdmissionPointProxyFramePending = FALSE;
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
    if (this->d->lodService && this->d->lodActiveGeneration != 0)
	this->d->lodService->cancelGeneration(this->d->lodActiveGeneration);
    this->d->lodActiveGeneration = 0;
    this->d->lodSubmissionSourceIndex = 0;
    this->d->lodSubmissionEntryOffset = 0;
    this->d->clearLodSubmissionPlan();
    this->d->lodSubmissionPending = FALSE;
    this->d->lodSubmissionRescanPending = FALSE;
    this->d->lodViewDemandPolicy.reset();
    this->d->lodDiscretePopulationTrialAvailable = FALSE;
    this->d->lodResidentGrowthPolicy.reset();
    this->d->lodResidentGrowthResidencyDrainActive = FALSE;
    this->d->lodResidentAdmissionRetryActive = FALSE;
    this->d->lodRetainedRefinementPending = FALSE;
    this->d->lodRetainedResidencyPending = FALSE;
    this->d->lodStructuralPresentationRepairPending = FALSE;
    this->d->lodStructuralRepairFrontierCount = 0;
    this->d->lodStructuralCoverageCostReservation = 0;
    this->d->lodRetainedRefinementCutAdvanced = FALSE;
    this->d->lodRetainedRefinementBudgetBlocked = FALSE;
    this->d->lodBudgetPolicy.resetCalibration();
    this->d->lodBudgetPolicy.clearRetainedQualityFloorBudget();
    this->d->lodPassAdmittedWork = FALSE;
    this->d->lodBudgetPolicy.resetPass();
    this->d->lodSettleAfterRenderSerial = 0;
    this->d->lodResultsPending.store(0);
    this->d->lodResultsFirstReadyMicroseconds.store(0);
    this->d->lodLastSubmittedViewRevision.reset();
    this->d->lodLastSubmittedPolicyRevision.reset();
    this->d->lodLastSubmittedSources.clear();
    this->d->lodSubmissionDeltaActive = FALSE;
    this->d->lodSubmissionDeltaSources.clear();
    this->d->lodSubmissionDeltaPlans.clear();
    this->d->lodFrameObligation.reset();
    this->d->lodPublicationPolicy.reset();
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodInterruptedPresentationReplayPending = FALSE;
    this->d->lodInteractiveProgressiveCeiling = -1;
    this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
    this->d->lodPointProxyCalibrationPolicy.reset();
    this->d->lodStructuralAdmissionPolicy.reset();
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodPointProxyTriangleRecoveryPending = FALSE;
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
    this->d->reconcilePhase(
	BObolLodStateMachine::Event::GENERATION_CANCELLED);
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
    this->d->lodStaticOverscanActive = FALSE;
    this->d->lodStaticOverscanRejected = FALSE;
    struct bv_lod_policy lodPolicy;
    bv_lod_policy_init(&lodPolicy);
    this->d->viewAttachment->getLodPolicy(&lodPolicy);
    this->d->lodCoveragePolicy.setRequired(
	this->d->lodAutoSubmit && this->d->lodService &&
	lodPolicy.policy != BV_LOD_OFF && lodPolicy.mesh_enabled);
    if (this->d->lodAutoSubmit) {
	this->requestRender("lod-auto-submit");
    } else {
	this->d->lodHeadroomPolicy.cancelRetry();
	this->d->lodResidentGrowthPolicy.reset();
	this->d->lodResidentGrowthResidencyDrainActive = FALSE;
	this->d->lodResidentAdmissionRetryActive = FALSE;
	this->d->lodInteractiveProgressiveCeiling = -1;
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->lodPointProxyCalibrationPolicy.reset();
	this->d->lodStructuralAdmissionPolicy.reset();
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	this->d->lodPointProxyTriangleRecoveryPending = FALSE;
	this->resetDiscoveryPointProxyFloor(FALSE);
	this->d->lodPresentationPolicy.reset();
	this->d->lodViewDemandPolicy.reset();
	this->d->lodDiscretePopulationTrialAvailable = FALSE;
	if (this->d->viewAttachment->getViewLodState()) {
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    viewState->setCadPresentationProgressiveCutCeiling(-1);
	    viewState->setCadPresentationPointProxyPixelThreshold(1.0f);
	    viewState->setCadPresentationCameraMotionFrameReuse(FALSE);
	}
    }
    this->d->reconcilePhase(
	BObolLodStateMachine::Event::AUTO_SUBMIT_CHANGED);
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

    if (this->d->lodUseForcedCut && this->d->lodForcedCut == cut)
	return;

    this->d->resetLodViewQualityHistory();
    this->d->lodUseForcedCut = TRUE;
    this->d->lodForcedCut = cut;
    this->advanceLodPolicyRevision();
    this->requestRender("lod-policy");
}

void
BObolViewController::clearLodForcedCut(void)
{
    if (!this->d->lodUseForcedCut)
	return;

    this->d->resetLodViewQualityHistory();
    this->d->lodUseForcedCut = FALSE;
    this->d->lodForcedCut = 0;
    this->advanceLodPolicyRevision();
    this->requestRender("lod-policy");
}

SbBool
BObolViewController::hasLodForcedCut(void) const
{
    return this->d->lodUseForcedCut;
}

int
BObolViewController::getLodForcedCut(void) const
{
    return this->d->lodForcedCut;
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

int
BObolViewController::consumeExportSourceFullDetail(
    SoBRLExportAction &exportAction,
    uint64_t generation,
    int *submittedRequestCount)
{
    return controller_consume_source_full_detail(this, exportAction,
	    generation, submittedRequestCount);
}

int
BObolViewController::consumeMeasureSourceFullDetail(
    SoBRLMeasureAction &measureAction,
    uint64_t generation,
    int *submittedRequestCount)
{
    return controller_consume_source_full_detail(this, measureAction,
	    generation, submittedRequestCount);
}

int
BObolViewController::consumeSnapSourceFullDetail(
    SoBRLSnapAction &snapAction,
    uint64_t generation,
    int *submittedRequestCount)
{
    return controller_consume_source_full_detail(this, snapAction,
	    generation, submittedRequestCount);
}

void
BObolViewController::clearRtPickCaches(void)
{
    for (size_t i = 0; i < this->d->rtPickCaches.size(); i++)
	delete this->d->rtPickCaches[i];
    this->d->rtPickCaches.clear();
    this->d->rtPickCachePaths.clear();
    this->d->rtPickCacheDatabases.clear();
    this->d->rtPickCacheSourceRevisions.clear();
}

int
BObolViewController::prepareRtPickCaches(void)
{
    std::vector<SbString> sourcePaths;
    std::vector<struct db_i *> sourceDatabases;
    std::vector<uint32_t> sourceRevisions;

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_sources(this);
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (!source || !source->getDatabase() ||
	    source->path.getValue().getLength() == 0)
	    continue;
	sourcePaths.push_back(source->path.getValue());
	sourceDatabases.push_back(source->getDatabase());
	sourceRevisions.push_back(
	    static_cast<uint32_t>(source->sourceRevision.getValue()));
    }

    SbBool sameSignature =
	sourcePaths.size() == this->d->rtPickCachePaths.size() &&
	sourceDatabases.size() == this->d->rtPickCacheDatabases.size() &&
	sourceRevisions.size() == this->d->rtPickCacheSourceRevisions.size() &&
	this->d->rtPickCaches.size() == this->d->rtPickCachePaths.size();
    if (sameSignature) {
	for (size_t i = 0; i < sourcePaths.size(); i++) {
	    if (sourceDatabases[i] != this->d->rtPickCacheDatabases[i] ||
		sourceRevisions[i] != this->d->rtPickCacheSourceRevisions[i] ||
		bu_strcmp(sourcePaths[i].getString(),
		       this->d->rtPickCachePaths[i].getString()) != 0 ||
		!this->d->rtPickCaches[i] ||
		!this->d->rtPickCaches[i]->isReady()) {
		sameSignature = FALSE;
		break;
	    }
	}
    }

    if (sameSignature)
	return static_cast<int>(this->d->rtPickCaches.size());

    this->clearRtPickCaches();
    for (size_t i = 0; i < sourcePaths.size(); i++) {
	std::vector<SbString> objectPaths;
	objectPaths.push_back(sourcePaths[i]);

	BObolRtPickCache *cache = new BObolRtPickCache;
	if (!cache->prepare(sourceDatabases[i], objectPaths)) {
	    delete cache;
	    continue;
	}

	this->d->rtPickCaches.push_back(cache);
	this->d->rtPickCachePaths.push_back(sourcePaths[i]);
	this->d->rtPickCacheDatabases.push_back(sourceDatabases[i]);
	this->d->rtPickCacheSourceRevisions.push_back(sourceRevisions[i]);
    }

    return static_cast<int>(this->d->rtPickCaches.size());
}

int
BObolViewController::getRtPickCacheCount(void) const
{
    return static_cast<int>(this->d->rtPickCaches.size());
}

BObolRtPickCache *
BObolViewController::getRtPickCache(int index) const
{
    if (index < 0 || static_cast<size_t>(index) >= this->d->rtPickCaches.size())
	return NULL;
    return this->d->rtPickCaches[static_cast<size_t>(index)];
}

uint32_t
BObolViewController::getRtPickCacheSourceRevision(int index) const
{
    if (index < 0 ||
	static_cast<size_t>(index) >=
	this->d->rtPickCacheSourceRevisions.size())
	return 0;
    return this->d->rtPickCacheSourceRevisions[static_cast<size_t>(index)];
}

int
BObolViewController::pickSourceMeshExactRay(
    BObolSourceMeshPickResult &pick,
    const SbVec3f &rayOrigin,
    const SbVec3f &rayDirection,
    uint64_t generation,
    int *submittedRequestCount)
{
    pick.clear();
    if (submittedRequestCount)
	*submittedRequestCount = 0;

    if (!this->d->viewport || !this->d->viewport->getRoot())
	return 0;

    BObolLodService *service = this->getLodService();
    if (!service || !service->isRunning())
	return 0;

    SoBRLSourceMeshPickAction sourcePickAction;
    sourcePickAction.setRay(rayOrigin, rayDirection);
    sourcePickAction.apply(this->d->viewport->getRoot());
    SbPlane clipPlanes[2];
    const size_t clipPlaneCount = this->getActiveClipPlanes(clipPlanes);

    const int requestCount =
	sourcePickAction.getSourceBackedFullDetailRequestCount();
    if (requestCount <= 0)
	return 0;

    std::vector<BObolLodRequest> expectedRequests;
    std::vector<BObolSourceMeshRequest> expectedSourceRequests;
    std::vector<int> expectedRequestIndices;
    std::vector<BObolSourceMeshRequest> submitSourceRequests;
    std::vector<SoBRLDatabaseSource *> requestSources;
    std::vector<int> submitRequestIndices;
    const int databaseSourceCount = controller_database_source_count(this);
    for (int i = 0; i < requestCount; i++) {
	const BObolSourceMeshRequest &sourceRequest =
	    sourcePickAction.getSourceBackedFullDetailRequest(i);
	SoBRLDatabaseSource *source =
	    controller_database_source_for_request(this, sourceRequest);
	if (source) {
	    BObolLodRequest requestTemplate;
	    controller_source_request_template(requestTemplate, source);

	    BObolLodRequest request;
	    if (sourcePickAction.makeSourceBackedFullDetailLodRequest(i,
		    request, &requestTemplate)) {
		expectedRequests.push_back(request);
		expectedSourceRequests.push_back(sourceRequest);
		expectedRequestIndices.push_back(i);
		submitSourceRequests.push_back(sourceRequest);
		requestSources.push_back(source);
		submitRequestIndices.push_back(i);
	    }
	}

	if (databaseSourceCount <= 1) {
	    BObolLodRequest request;
	    if (sourcePickAction.makeSourceBackedFullDetailLodRequest(i,
		    request)) {
		expectedRequests.push_back(request);
		expectedSourceRequests.push_back(sourceRequest);
		expectedRequestIndices.push_back(i);
	    }
	}
    }
    if (expectedRequests.empty() && submitSourceRequests.empty())
	return 0;

    std::vector<SbBool> requestMatched(
	static_cast<size_t>(requestCount), FALSE);
    if (!expectedRequests.empty()) {
	std::vector<BObolLodResult> sourceResults;
	service->drainMatchingResults(sourceResults, expectedRequests);
	if (!sourceResults.empty()) {
	    std::vector<SbBool> used(sourceResults.size(), FALSE);
	    std::vector<SbBool> requestConsumed(
		static_cast<size_t>(requestCount), FALSE);
	    for (size_t i = 0; i < expectedRequests.size(); i++) {
		const int requestIndex = expectedRequestIndices[i];
		if (requestIndex < 0 ||
		    static_cast<size_t>(requestIndex) >=
		    requestConsumed.size() ||
		    requestConsumed[static_cast<size_t>(requestIndex)])
		    continue;

		for (size_t j = 0; j < sourceResults.size(); j++) {
		    if (used[j] ||
			!bobol_lod_result_matches_request(
			    sourceResults[j], expectedRequests[i]))
			continue;

		    used[j] = TRUE;
		    requestMatched[static_cast<size_t>(requestIndex)] =
			TRUE;
		    BObolSourceMeshPickResult candidate;
		    if (bobol_pick_source_full_detail_result(candidate,
			    expectedSourceRequests[i], sourceResults[j],
			    sourcePickAction.getRayOrigin(),
			    sourcePickAction.getRayDirection(), clipPlanes,
			    clipPlaneCount)) {
			requestConsumed[static_cast<size_t>(requestIndex)] =
			    TRUE;
			if (!pick.hit || candidate.distance < pick.distance)
			    pick = candidate;
		    }
		    break;
		}
	    }
	}
    }

    int submitted = 0;
    for (size_t i = 0; i < submitSourceRequests.size(); i++) {
	const int requestIndex = submitRequestIndices[i];
	if (requestIndex >= 0 &&
	    static_cast<size_t>(requestIndex) < requestMatched.size() &&
	    requestMatched[static_cast<size_t>(requestIndex)])
	    continue;

	BObolLodRequest requestTemplate;
	controller_source_request_template(requestTemplate, requestSources[i]);
	if (bobol_lod_submit_rt_source_full_detail_request(service,
		generation, submitSourceRequests[i],
		requestSources[i]->getDatabase(), &requestTemplate,
		this->getMaxExactFullDetailFaceCount(),
		this->getMaxExactFullDetailPointCount()) != 0)
	    submitted++;
    }
    if (submittedRequestCount)
	*submittedRequestCount = submitted;

    if (pick.hit)
	return 1;
    pick.clear();
    return 0;
}

int
BObolViewController::pickRtExactRay(
    std::vector<BObolRtPickResult> &results,
    const SbVec3f &rayOrigin,
    const SbVec3f &rayDirection,
    SbBool pickAll)
{
    results.clear();

    SbVec3f direction = rayDirection;
    if (direction.length() <= 0.0f)
	return 0;
    direction.normalize();

    const int cacheCount = this->prepareRtPickCaches();
    SbPlane clipPlanes[2];
    const size_t clipPlaneCount = this->getActiveClipPlanes(clipPlanes);
    for (int i = 0; i < cacheCount; i++) {
	BObolRtPickCache *cache = this->getRtPickCache(i);
	if (!cache || !cache->isReady())
	    continue;

	BObolRtPickResult rtPick;
	if (!cache->pickRay(rtPick, rayOrigin, direction, clipPlanes,
		clipPlaneCount) || !rtPick.hit)
	    continue;
	if (rt_pick_result_path_recorded(results, rtPick))
	    continue;

	insert_rt_pick_result(results, rtPick, pickAll ? TRUE : FALSE);
    }

    return static_cast<int>(results.size());
}

void
BObolViewController::setMeshResidencyBudget(
    size_t maxBytes,
    SbBool evictDisplayPayloads)
{
    this->d->meshResidencyBudgetEnabled = TRUE;
    this->d->maxResidentMeshBytes = maxBytes;
    this->d->meshResidencyEvictDisplayPayloads =
	evictDisplayPayloads ? TRUE : FALSE;
}

void
BObolViewController::clearMeshResidencyBudget(void)
{
    this->d->meshResidencyBudgetEnabled = FALSE;
    this->d->maxResidentMeshBytes = 0;
    this->d->meshResidencyEvictDisplayPayloads = TRUE;
}

SbBool
BObolViewController::hasMeshResidencyBudget(void) const
{
    return this->d->meshResidencyBudgetEnabled;
}

size_t
BObolViewController::getMaxResidentMeshBytes(void) const
{
    return this->d->maxResidentMeshBytes;
}

SbBool
BObolViewController::isMeshResidencyDisplayEvictionEnabled(void) const
{
    return this->d->meshResidencyEvictDisplayPayloads;
}

size_t
BObolViewController::evictMeshPayloadsToBudget(
    size_t maxBytes,
    SbBool evictDisplayPayloads)
{
    this->d->lastMeshBudgetInitialResidentBytes = 0;
    this->d->lastMeshBudgetFinalResidentBytes = 0;
    this->d->lastMeshBudgetFreedFullDetailBytes = 0;
    this->d->lastMeshBudgetFreedDisplayBytes = 0;
    this->d->lastMeshBudgetVisitedMeshCount = 0;
    this->d->lastMeshBudgetEvictedFullDetailMeshCount = 0;
    this->d->lastMeshBudgetEvictedDisplayMeshCount = 0;

    SoNode *root = this->getRenderSceneRoot();
    if (!root)
	root = this->getSceneRoot();
    if (!root)
	return 0;

    SoBRLMeshResidencyAction action;
    action.setMaxResidentMeshBytes(maxBytes);
    action.setEvictDisplayPayloads(evictDisplayPayloads);
    action.apply(root);

    const size_t viewLodBytes = this->d->viewAttachment->getViewLodState() ?
				this->d->viewAttachment->getViewLodState()->estimateDisplayMeshBytes() : 0;
    this->d->lastMeshBudgetInitialResidentBytes =
	action.getInitialResidentMeshBytes() + viewLodBytes;
    this->d->lastMeshBudgetFinalResidentBytes =
	action.getFinalResidentMeshBytes() + viewLodBytes;
    this->d->lastMeshBudgetFreedFullDetailBytes =
	action.getFreedFullDetailBytes();
    this->d->lastMeshBudgetFreedDisplayBytes =
	action.getFreedDisplayBytes();
    this->d->lastMeshBudgetVisitedMeshCount = action.getVisitedMeshCount();
    this->d->lastMeshBudgetEvictedFullDetailMeshCount =
	action.getEvictedFullDetailMeshCount();
    this->d->lastMeshBudgetEvictedDisplayMeshCount =
	action.getEvictedDisplayMeshCount();

    if (this->d->lastMeshBudgetFinalResidentBytes > maxBytes &&
	this->d->viewAttachment->getViewLodState()) {
	std::vector<SoBRLMeshShape *> shapes =
	    controller_render_mesh_shapes(this);
	for (size_t i = 0; i < shapes.size(); i++) {
	    if (this->d->lastMeshBudgetFinalResidentBytes <= maxBytes)
		break;
	    if (!this->d->viewAttachment->getViewLodState()->findMesh(shapes[i]))
		continue;
	    size_t freed = shapes[i]->evictDisplayMeshPreservingSourceMetrics();
	    if (freed == 0)
		continue;
	    this->d->lastMeshBudgetFreedFullDetailBytes += freed;
	    this->d->lastMeshBudgetEvictedFullDetailMeshCount++;
	    this->d->lastMeshBudgetFinalResidentBytes =
		this->d->lastMeshBudgetFinalResidentBytes > freed ?
		this->d->lastMeshBudgetFinalResidentBytes - freed : 0;
	}
    }

    if (evictDisplayPayloads &&
	this->d->lastMeshBudgetFinalResidentBytes > maxBytes &&
	this->d->viewAttachment->getViewLodState()) {
	unsigned int evicted = 0;
	size_t freed = this->d->viewAttachment->getViewLodState()->
	    evictDisplayMeshPayloads(&evicted);
	if (freed > 0) {
	    this->d->lastMeshBudgetFreedDisplayBytes += freed;
	    this->d->lastMeshBudgetEvictedDisplayMeshCount += evicted;
	    this->d->lastMeshBudgetFinalResidentBytes =
		this->d->lastMeshBudgetFinalResidentBytes > freed ?
		this->d->lastMeshBudgetFinalResidentBytes - freed : 0;
	}
    }

    size_t freedBytes = this->getLastMeshBudgetFreedResidentBytes();
    if (freedBytes > 0)
	this->requestRender("lod-memory-budget");
    return freedBytes;
}

size_t
BObolViewController::enforceMeshResidencyBudget(void)
{
    if (!this->d->meshResidencyBudgetEnabled)
	return 0;

    return this->evictMeshPayloadsToBudget(
	       this->d->maxResidentMeshBytes,
	       this->d->meshResidencyEvictDisplayPayloads);
}

size_t
BObolViewController::getLastMeshBudgetInitialResidentBytes(void) const
{
    return this->d->lastMeshBudgetInitialResidentBytes;
}

size_t
BObolViewController::getLastMeshBudgetFinalResidentBytes(void) const
{
    return this->d->lastMeshBudgetFinalResidentBytes;
}

size_t
BObolViewController::getLastMeshBudgetFreedResidentBytes(void) const
{
    return this->d->lastMeshBudgetInitialResidentBytes >
	   this->d->lastMeshBudgetFinalResidentBytes ?
	   this->d->lastMeshBudgetInitialResidentBytes -
	   this->d->lastMeshBudgetFinalResidentBytes : 0;
}

size_t
BObolViewController::getLastMeshBudgetFreedFullDetailBytes(void) const
{
    return this->d->lastMeshBudgetFreedFullDetailBytes;
}

size_t
BObolViewController::getLastMeshBudgetFreedDisplayBytes(void) const
{
    return this->d->lastMeshBudgetFreedDisplayBytes;
}

unsigned int
BObolViewController::getLastMeshBudgetVisitedMeshCount(void) const
{
    return this->d->lastMeshBudgetVisitedMeshCount;
}

unsigned int
BObolViewController::getLastMeshBudgetEvictedFullDetailMeshCount(void) const
{
    return this->d->lastMeshBudgetEvictedFullDetailMeshCount;
}

unsigned int
BObolViewController::getLastMeshBudgetEvictedDisplayMeshCount(void) const
{
    return this->d->lastMeshBudgetEvictedDisplayMeshCount;
}

SbBool
BObolViewController::hasPendingLodResults(void) const
{
    return this->d->lodResultsPending.load() != 0 ? TRUE : FALSE;
}

SbBool
BObolViewController::hasPendingLodSubmissions(void) const
{
    return this->d->lodSubmissionPending;
}

SbBool
BObolViewController::hasPendingLodRefinementFrame(void) const
{
    /* A budget calibration probe is just as much a refinement barrier as a
     * newly selected PoP cut.  Reporting idle in the interval between
     * requesting that frame and completeRenderTiming() made warm scenes stop
     * at cache-temperature-dependent coarse cuts. */
    return this->d->lodFrameObligation.pending() ||
	this->d->lodBudgetPolicy.rescanAfterFrame() ||
	(this->d->lodPresentationPolicy.handoffPending() &&
	 !this->d->lodSubmissionPending) ||
	this->d->lodAdmissionPointProxyFramePending ||
	this->d->lodStablePointProxyCalibrationPending ||
	this->d->lodHeadroomPolicy.retryPending();
}

size_t
BObolViewController::processPendingLodResults(size_t maxResults,
	uint64_t maxMicroseconds)
{
    if (!this->d->lodService)
	return 0;

    const auto clear_lod_wakeup_if_idle = [this]() {
	if (this->d->lodSubmissionPending ||
	    this->d->lodResultsPending.load() != 0 ||
	    this->d->lodResidentGrowthPolicy.pending() ||
	    (this->d->lodService &&
	     this->d->lodService->queuedResultCountForGeneration(
		 this->d->lodActiveGeneration) > 0))
	    return;

	/* Merely having a registered provider does not mean that provider has
	 * work.  Keeping the latch set for the lifetime of every GED provider
	 * made a late result-ready callback leave an otherwise idle warm scene
	 * polling forever.  Provider callbacks report their own hasMore state
	 * from advanceProgressiveWork().
	 *
	 * A result-ready callback may race this drain.  Clear first, then recheck
	 * so a concurrent callback cannot lose its frame wakeup. */
	this->clearProgressiveWorkPending();
	if (this->d->lodSubmissionPending ||
	    this->d->lodResultsPending.load() != 0 ||
	    this->d->lodResidentGrowthPolicy.pending() ||
	    (this->d->lodService &&
	     this->d->lodService->queuedResultCountForGeneration(
		 this->d->lodActiveGeneration) > 0))
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
     * A first useful mesh and an existing-mesh refinement barrier still yield
     * after one quantum so the user sees that state before later generations.
     * Explicit callers which pass a zero time budget retain the unbounded
     * drain path above.
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
	if (firstUsefulCoverage ||
	    this->d->lodFrameObligation.pending())
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
	int reset)
{
    if (!this->d->lodService || !this->d->lodService->isRunning())
	return 0;
    if (!this->d->activeCamera ||
	(!this->getSceneRoot() && !this->getRenderSceneRoot()))
	return 0;

    this->syncLodViewSignature(TRUE);

    /* The discovery threshold is derived from an exact structural frame and
	 * therefore gates first-wave provider admission.  Stable point calibration
	 * is different: it classifies already discovered occurrences and must never
	 * block the source pass which creates or repopulates its CAD assembly.  It
	 * continues to gate convergence and is consumed only by a completed frame. */
    if (BObolLodPointProxyCalibrationPolicy::blocksSourceAdmission(
	    this->d->lodAdmissionPointProxyFramePending != FALSE,
	    this->d->lodStablePointProxyCalibrationPending != FALSE))
	return 0;

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
	this->d->lodSubmissionPending = FALSE;
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionDeltaActive = FALSE;
	this->d->lodSubmissionDeltaSources.clear();
	this->d->lodSubmissionDeltaPlans.clear();
	/* The empty-contract edge retires the source identity domain, not merely
	 * the current submission cursor.  A Z/redraw may pass through this state
	 * before publishing a replacement source; preserving the old dense census
	 * then adds the replacement population to an unreachable source entry. */
	this->d->lodCoveragePolicy.reset();
	this->d->clearLodConvergenceCandidates();
	this->d->lodProjectedDemandCache.clear();
	this->d->lodResidentGrowthPolicy.reset();
	this->d->lodResidentGrowthResidencyDrainActive = FALSE;
	this->d->lodResidentAdmissionRetryActive = FALSE;
	this->d->lodViewDemandPolicy.clearDemandRefresh();
	this->d->lodLastSubmittedSources.clear();
	return 0;
    }

    if (this->d->lodLastSubmittedViewRevision == this->d->lodViewRevision &&
	this->d->lodLastSubmittedPolicyRevision == this->d->lodPolicyRevision &&
	signatures.sameInventories(this->d->lodLastSubmittedSources)) {
	if (this->d->lodSubmissionPending) {
	    /* A pending cursor is not meaningful without an owned generation.
	     * Scene replacement normally clears both together, but a late
	     * resident-growth/coverage edge can legitimately re-arm submission
	     * after cancellation.  The unchanged-signature fast path formerly
	     * returned forever in that state: no task, result, frame, or camera
	     * event remained to create another generation. */
	    if (this->d->lodActiveGeneration == 0) {
		this->d->lodActiveGeneration =
		    this->d->lodService->beginGeneration();
		this->d->lodResultsPending.store(0);
		this->d->lodResultsFirstReadyMicroseconds.store(0);
		if (this->d->lodActiveGeneration == 0)
		    return 0;
		this->d->reconcilePhase(
		    BObolLodStateMachine::Event::WORK_SCHEDULED);
	    }
	    return this->submitLodRequests(this->d->lodService,
		this->d->lodActiveGeneration,
		this->d->lodSubmissionRefreshMissing,
		this->d->lodSubmissionReset);
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
	this->d->lodResidentGrowthResidencyDrainActive = FALSE;
    /* The first source contract is submitted immediately.  Once that
     * time-to-first-mesh edge exists, allow an append-only producer to build
     * a bounded inventory wave before starting another owner-thread LoD
     * transaction.  All semantic invalidations and the producer's final edge
     * bypass this gate.  Because lodLastSubmittedSources remains unchanged
     * while deferred, the source journal supplies the complete accumulated
     * delta when the deadline expires. */
    const bool deferInventoryDelta =
	!sourceSetChanged && !viewOrPolicyChanged && inventoryChanged &&
	this->d->lodInventoryDeltaPolicy.defer(
	    true, this->d->progressiveProviderPendingCount > 0,
	    this->d->lodSubmissionPending != FALSE,
	    !this->d->lodLastSubmittedSources.empty(),
	    this->d->lodInteractive != FALSE, bu_gettime());
    if (deferInventoryDelta) {
	this->markProgressiveWorkPending();
	return 0;
    }
    this->d->lodInventoryDeltaPolicy.committed();
    if (sourceSetChanged || inventoryChanged || viewOrPolicyChanged)
	this->d->lodResidentAdmissionRetryActive = FALSE;
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
	this->d->lodSubmissionPending &&
	this->d->lodSubmissionDeltaActive &&
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
	    const auto existingPlan = std::find_if(
		this->d->lodSubmissionDeltaPlans.begin(),
		this->d->lodSubmissionDeltaPlans.end(),
		[changedSource](const auto &entry) {
		    return entry.first == changedSource;
		});
	    if (existingPlan ==
		    this->d->lodSubmissionDeltaPlans.end()) {
		/* No selective plan means this source is already being scanned
		 * in full.  Its identity plan is extended from compactCount in
		 * submitLodRequests below. */
		const bool alreadyTargeted =
		    std::find(this->d->lodSubmissionDeltaSources.begin(),
			this->d->lodSubmissionDeltaSources.end(),
			changedSource) !=
		    this->d->lodSubmissionDeltaSources.end();
		if (!alreadyTargeted)
		    pendingDeltaNeedsFullRescan = true;
		else
		    extendedPendingDelta = true;
		continue;
	    }
	    existingPlan->second.insert(existingPlan->second.end(),
		changedEntries.begin(), changedEntries.end());
	    if (this->d->lodSubmissionPlanValid &&
		this->d->lodSubmissionPlanSource == changedSource)
		this->d->lodSubmissionPlanEntries.insert(
		    this->d->lodSubmissionPlanEntries.end(),
		    changedEntries.begin(), changedEntries.end());
	    extendedPendingDelta = true;
	}
    }
    if ((sourceSetChanged || inventoryChanged) &&
	!viewOrPolicyChanged &&
	!this->d->lodSubmissionPending &&
	!this->d->lodLastSubmittedSources.empty()) {
	if (!this->d->lodSubmissionDeltaActive) {
	    this->d->lodSubmissionDeltaSources.clear();
	    this->d->lodSubmissionDeltaPlans.clear();
	}
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
		std::find(this->d->lodSubmissionDeltaSources.begin(),
		    this->d->lodSubmissionDeltaSources.end(),
		    changedSource) !=
		this->d->lodSubmissionDeltaSources.end();
	    if (!alreadyTargeted)
		this->d->lodSubmissionDeltaSources.push_back(changedSource);
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
		    auto existingPlan = std::find_if(
			this->d->lodSubmissionDeltaPlans.begin(),
			this->d->lodSubmissionDeltaPlans.end(),
			[changedSource](const auto &entry) {
			    return entry.first == changedSource;
			});
		    /* An existing target without a delta plan already requires
		     * a full source scan; a later selective update must not
		     * accidentally narrow it. */
		    if (existingPlan !=
			    this->d->lodSubmissionDeltaPlans.end()) {
			existingPlan->second.insert(
			    existingPlan->second.end(),
			    changedEntries.begin(), changedEntries.end());
			std::sort(existingPlan->second.begin(),
			    existingPlan->second.end());
			existingPlan->second.erase(std::unique(
			    existingPlan->second.begin(),
			    existingPlan->second.end()),
			    existingPlan->second.end());
		    } else if (!alreadyTargeted) {
			this->d->lodSubmissionDeltaPlans.push_back(
			    std::make_pair(changedSource,
				std::move(changedEntries)));
		    }
		} else {
		    sourceDeltaInvalidatesCoverage = true;
		    this->d->lodSubmissionDeltaPlans.erase(
			std::remove_if(
			    this->d->lodSubmissionDeltaPlans.begin(),
			    this->d->lodSubmissionDeltaPlans.end(),
			    [changedSource](const auto &entry) {
				return entry.first == changedSource;
			    }),
			this->d->lodSubmissionDeltaPlans.end());
		}
	    } else {
		this->d->lodSubmissionDeltaPlans.erase(
		    std::remove_if(
			this->d->lodSubmissionDeltaPlans.begin(),
			this->d->lodSubmissionDeltaPlans.end(),
			[changedSource](const auto &entry) {
			    return entry.first == changedSource;
			}),
		    this->d->lodSubmissionDeltaPlans.end());
	    }
	}
	useSourceDelta = !this->d->lodSubmissionDeltaSources.empty();
	this->d->lodSubmissionDeltaActive =
	    useSourceDelta ? TRUE : FALSE;
    } else if ((sourceSetChanged || inventoryChanged) &&
	!this->d->lodSubmissionPending) {
	this->d->lodSubmissionDeltaActive = FALSE;
	this->d->lodSubmissionDeltaSources.clear();
	this->d->lodSubmissionDeltaPlans.clear();
    }
    const bool hasExactSourceDelta =
	useSourceDelta || this->d->lodSubmissionDeltaActive;
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
	this->d->lodRetainPoseOccurrenceCuts = FALSE;
	this->d->lodRetainedImportanceCensusPending = FALSE;
	this->d->lodBudgetPolicy.clearRetainedQualityFloorBudget();
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
	this->d->lodPointProxyCalibrationPolicy.reset();
	const SbBool pointRelaxationRequired =
	    this->d->lodPresentationPointProxyPixelThreshold > 1.01f ?
		TRUE : FALSE;
	if (this->d->lodPointProxyTriangleRecoveryPending) {
	    /* Structural streaming extended the population during the bounded
	     * recovery.  Keep its terminal obligation and its evidence provenance
	     * when restarting the complete-source pass.  Converting a one-pass
	     * handoff normalization into a persistent measured ceiling freezes a
	     * warm cache at whichever tiny prefix happened to arrive first. */
	    if (this->d->lodBudgetPolicy.retainedRecoveryCeilingActive())
		this->d->lodBudgetPolicy.requestRetainedRecovery(
		    this->d->lodBudgetPolicy.currentBudget());
	    else
		this->d->lodBudgetPolicy.requestRetainedNormalization(
		    this->d->lodBudgetPolicy.currentBudget());
	    this->d->lodStablePointProxyCalibrationPending = FALSE;
	} else {
	    this->d->lodStablePointProxyCalibrationPending =
		pointRelaxationRequired;
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
	!this->d->lodSubmissionPending) {
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
	this->d->lodSubmissionRescanPending =
	    useSourceDelta && sourceCoverageInvalidated ? TRUE : FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedResidencyPending = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodBudgetPolicy.resetCalibration();
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	if (useSourceDelta)
	    this->d->lodCoveragePolicy.clearPassCounters();
	if (sourceSetChanged || inventoryChanged)
	    this->d->lodBudgetPolicy.resetOverloadRecovery();
    } else if (!extendedPendingDelta) {
	this->d->lodSubmissionRescanPending = TRUE;
    }
    if (pendingDeltaNeedsFullRescan)
	this->d->lodSubmissionRescanPending = TRUE;
    this->d->lodSubmissionPending = TRUE;
    this->d->lodBudgetPolicy.clearBudgetLimit();
    this->d->lodSubmissionRefreshMissing = refreshMissing;
    this->d->lodSubmissionReset = reset;
    int submitted = this->submitLodRequests(this->d->lodService, generation,
					    refreshMissing, reset);
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
	int reset)
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

    if (!this->d->lodSubmissionPending) {
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedResidencyPending = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	this->d->lodCoveragePolicy.clearPassCounters();
	this->d->lodSubmissionPending = TRUE;
	this->d->lodSubmissionRefreshMissing = refreshMissing;
	this->d->lodSubmissionReset = reset;
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
     * beginPass block: one retained transaction spans many bounded source
     * windows, and its completion/handoff uses the same producer state. */
    const bool retainedProviderInventorySettled =
	this->d->progressiveProviderPendingCount == 0 &&
	!controller_lod_compact_inventory_incomplete(sources);
    const bool retainedServiceStreamIdle = !service ||
	(service->activeTaskCountForGeneration(generation) == 0 &&
	 service->queuedResultCountForGeneration(generation) == 0);
    const bool retainedResultDeliveryIdle =
	this->d->lodResultsPending.load() == 0 &&
	!this->d->lodPublicationPolicy.pending();
    const bool retainedPopulationSettled =
	BObolLodResidentGrowthPolicy::allocationPopulationSettled(
	    retainedProviderInventorySettled, retainedServiceStreamIdle,
	    retainedResultDeliveryIdle,
	    this->d->lodResidentGrowthPolicy.pending(),
	    this->d->lodResidentGrowthResidencyDrainActive != FALSE);
    if (!this->d->lodBudgetPolicy.passInitialized()) {
	const size_t presentationReconciliationBudget =
	    this->d->lodPresentationPolicy.handoffReconciliationBudget();
	if (presentationReconciliationBudget)
	    this->d->lodBudgetPolicy.requestPresentationReconciliation(
		presentationReconciliationBudget);
	this->d->lodMissingMeshBudgetBlockedCount = 0;
	/* A retained minimax allocation is one complete occurrence-plan
	 * transaction, even when result publication or an endpoint deadline
	 * inserts a presentation between bounded planning slices.  Those frames
	 * reset the per-frame allowance, but they must not silently turn the
	 * unconsumed tail of the plan into an ordinary first-come refinement
	 * pass.  Re-arm the same measured-budget allocation until completedPass
	 * retires the mode below. */
	if (this->d->lodSubmissionRetainedAdmissionMode &&
	    this->d->lodSubmissionPending)
	    this->d->lodBudgetPolicy.requestRetainedReallocation();
	const BObolViewLodState *lodState =
	    this->d->viewAttachment->getViewLodState();
	const size_t activeFaces = lodState ? lodState->activeFaceCount() : 0;
	size_t activeCost = 0;
	size_t minimumActiveCost = 0;
	controller_lod_effective_population_cost(lodState, activeCost,
	    minimumActiveCost);
	const SbBool scaleQualityProbe =
	    this->d->lodInteractive &&
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
	    !this->d->lodInteractive &&
	    (this->d->lodStructuralPresentationRepairPending ||
	     this->d->lodStaticOverscanActive ||
	     this->d->lodStaticOverscanRejected) ? TRUE : FALSE;
	const float targetFps =
	    scaleQualityProbe ?
		BObolLodViewDemandPolicy::
		    qualityTargetFramesPerSecond() :
		(this->d->lodInteractive ?
		    this->d->lodInteractiveTargetFps :
		    (hardDeadlinePresentation ?
			this->d->staticQualityTargetFps() :
			this->d->quietAllocationTargetFps()));
	const long double calibratedCostPerSecond =
	    scaleQualityProbe ?
		this->d->lodStableCalibratedRenderCostPerSecond :
		(this->d->lodInteractive ?
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
	BObolLodBudgetPolicy::Inputs budgetInputs;
	budgetInputs.activeCost = activeCost;
	budgetInputs.minimumActiveCost = minimumActiveCost;
	budgetInputs.targetFps = targetFps;
	budgetInputs.calibratedCostPerSecond = calibratedCostPerSecond;
	budgetInputs.observedStableNanoseconds = observedStableNanoseconds;
	budgetInputs.lastRenderNanoseconds =
	    this->d->lastRenderTimeNanoseconds;
	budgetInputs.smoothedRenderNanoseconds =
	    this->d->smoothedRenderTimeNanoseconds;
	budgetInputs.interactive = this->d->lodInteractive != FALSE;
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
	budgetInputs.structuralCoverageRepair =
	    this->d->lodStructuralPresentationRepairPending != FALSE;
	budgetInputs.forceTerminal =
	    this->d->forceTerminalLodRefinement != FALSE;
	budgetInputs.releaseCutFloor =
	    this->d->lodReleaseCutFloorActive != FALSE;
	budgetInputs.stablePresentationHandoff =
		    this->d->lodPresentationPolicy.handoffPending();
	budgetInputs.stablePresentationCostFloor =
		    this->d->lodPresentationPolicy.handoffCostFloor();
    const BObolLodBudgetPolicy::Decision budget =
	    this->d->lodBudgetPolicy.beginPass(budgetInputs);
    if (this->d->lodStructuralPresentationRepairPending &&
	!this->d->lodStructuralCoverageCostReservation)
	this->d->lodStructuralCoverageCostReservation =
	    BObolLodStructuralAdmissionPolicy::perOccurrenceReservation(
		budget.totalBudget, this->d->lodBudgetPolicy.passActiveCost(),
		this->d->lodStructuralRepairFrontierCount);
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
	    !this->d->lodSubmissionRetainedAdmissionMode) {
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRetainedAdmissionMode = TRUE;
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	} else if (!beginRetainedAdmission &&
	    this->d->lodSubmissionRetainedAdmissionMode) {
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRetainedAdmissionMode = FALSE;
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

    const bool effectiveRetainedAdmission =
	this->d->lodBudgetPolicy.retainedAdmission() &&
	retainedPopulationSettled;

    /* The ordinary quiet path reaches the raster-stable tier through measured
     * 1.0 -> 0.75/0.5 -> 0.25 pixel admissions.  A deterministic terminal
     * caller has deliberately removed those capacity ceilings, so select the
     * same final physical target directly.  This remains view-dependent and
     * does not imply loading the hierarchy's full topology. */
    float scenePixelError = this->d->forceTerminalLodRefinement ?
	std::min(this->d->lodTargetPixelError, 0.25f) :
	this->d->lodTargetPixelError;
    const BObolViewLodState *sceneLodState =
	this->d->viewAttachment->getViewLodState();
    const size_t sceneActiveFaces = sceneLodState ?
	sceneLodState->activeFaceCount() : 0;
    /* beginPass() owns one complete population census.  Reuse its frozen
     * currencies across every bounded compact window; rescanning all retained
     * occurrences here made a 150k scene spend roughly three quarters of its
     * CPU in point-proxy cost discovery instead of publishing geometry. */
    const size_t sceneActiveCost =
	this->d->lodBudgetPolicy.passActiveCost();
    const size_t sceneMinimumActiveCost =
	this->d->lodBudgetPolicy.passMinimumActiveCost();
    if (!this->d->lodUseForcedCut &&
	this->d->lodBudgetPolicy.currentBudget() != SIZE_MAX &&
	this->d->lodBudgetPolicy.currentBudget() > 0 &&
	sceneActiveCost > this->d->lodBudgetPolicy.currentBudget()) {
	const long double over =
	    static_cast<long double>(sceneActiveCost) /
	    static_cast<long double>(this->d->lodBudgetPolicy.currentBudget());
	scenePixelError *= static_cast<float>(std::sqrt(over));
    }

    /* Active motion uses the renderer's O(1) aggregate ceiling and leaves
     * occurrence cuts intact while zoom prefetches missing resident suffixes.
     * A complete minimax reallocation is stable-view policy: running it in
     * the first wheel callback blocked input for hundreds of milliseconds at
     * 50k leaves and duplicated work the quiet demand census must perform
     * against the final camera anyway. */
    const SbBool retainedAllocationPass =
	effectiveRetainedAdmission &&
	!this->d->lodInteractive ? TRUE : FALSE;
    if (retainedAllocationPass &&
	this->d->lodSubmissionSourceIndex == 0 &&
	this->d->lodSubmissionEntryOffset == 0) {
	size_t maximumProtectedBudget =
	    this->d->lodBudgetPolicy.currentBudget();
	size_t maximumMarginalBudget = maximumProtectedBudget;
	if (!this->d->lodInteractive &&
	    !this->d->lodStaticOverscanRejected &&
	    !this->d->lodBudgetPolicy.retainedQualityFloorRejected() &&
	    this->d->lodStableCalibratedRenderCostPerSecond > 0.0L) {
	    const long double hardQualityFps =
		static_cast<long double>(BObolLodViewDemandPolicy::
		    qualityTargetFramesPerSecond());
	    const long double affordable = hardQualityFps > 0.0L ?
		this->d->lodStableCalibratedRenderCostPerSecond * 0.80L /
		    hardQualityFps : 0.0L;
	    const size_t hardQualityBudget =
		affordable >= static_cast<long double>(SIZE_MAX) ? SIZE_MAX :
		(affordable > 0.0L ? static_cast<size_t>(affordable) : 0);
	    maximumProtectedBudget = std::max(
		maximumProtectedBudget, hardQualityBudget);
	}
	/* A failed protected-floor frame is evidence against that exact complete
	 * population, not against every smaller static-quality improvement.  Keep
	 * the floor disabled, but let the occurrence-local allocator spend a
	 * safety-margin-capped fraction of the measured static capacity.  The
	 * deadline ceiling prevents it from retrying the rejected population. */
	if (!this->d->lodInteractive &&
	    (this->d->lodStaticOverscanRejected ||
	     this->d->lodBudgetPolicy.retainedQualityFloorRejected()) &&
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
		this->d->lodBudgetPolicy.deadlineCapacityCeiling();
	    safeBudget = BObolLodStaticQualityPolicy::capMarginalBudget(
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
	inputs.sceneBudget = this->d->lodBudgetPolicy.currentBudget();
	inputs.maximumMarginalBudget = maximumMarginalBudget;
	inputs.allowProtectedFloor = !this->d->lodInteractive &&
	    !this->d->lodStaticOverscanRejected &&
	    !this->d->lodBudgetPolicy.retainedQualityFloorRejected() &&
	    this->d->lodPresentationPolicy.
		handoffReconciliationBudget() == 0;
	inputs.maximumProtectedBudget = maximumProtectedBudget;
	inputs.viewRevision = this->d->lodViewRevision.value();
	inputs.policyRevision = this->d->lodPolicyRevision.value();
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
	    priorAllocation.viewRevision == inputs.viewRevision &&
	    priorAllocation.policyRevision == inputs.policyRevision &&
	    std::isfinite(priorAllocation.pointProxyPixelThreshold) &&
	    std::fabs(priorAllocation.pointProxyPixelThreshold -
		inputs.pointProxyPixelThreshold) <= 1.0e-6f &&
	    priorAllocation.requestedSceneBudget == inputs.sceneBudget &&
	    priorAllocation.externalPresentationCost ==
		inputs.externalPresentationCost &&
	    priorAllocation.maximumMarginalBudget ==
		inputs.maximumMarginalBudget &&
	    priorAllocation.maximumProtectedBudget ==
		inputs.maximumProtectedBudget &&
	    priorAllocation.allowProtectedFloor == inputs.allowProtectedFloor;
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
	if (allocationStatus == BOBOL_RETAINED_ALLOCATION_PENDING ||
	    allocationStatus == BOBOL_RETAINED_ALLOCATION_STALE) {
	    /* Preserve the prior coherent presentation while this unpublished
	     * plan advances.  The ordinary progressive pump owns the timer edge;
	     * no render is needed merely to calculate another slice. */
	    this->d->lodSubmissionPending = TRUE;
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
	this->d->lodBudgetPolicy.setRetainedQualityFloorBudget(
	    protectedFloorBudget, protectedFloorSignature, sceneActiveCost,
	    sceneMinimumActiveCost);
	if (getenv("BOBOL_LOD_TRACE_BUDGET") ||
	    getenv("BOBOL_LOD_TRACE_ALLOCATOR"))
	    bu_log("BObol LoD retained importance ceiling=%.6f "
		   "budget=%zu selected=%zu certified=%zu "
		   "marginal_limit=%zu protected_limit=%zu reused=%d "
		   "elapsed_us=%lld\n",
		   this->d->lodRetainedAdmissionMaximumNormalizedError,
		   this->d->lodBudgetPolicy.currentBudget(),
		   this->d->lodRetainedAllocationCertificate.
		       selectedPresentationCost,
		   this->d->lodRetainedAllocationCertificate.
		       certifiedPresentationBudget,
		   maximumMarginalBudget,
		   maximumProtectedBudget, reuseCoveredAllocation ? 1 : 0,
		   static_cast<long long>(retainedAllocationElapsed));
    }
    SbBool boundedScenePass = FALSE;
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
	if (this->d->lodSubmissionDeltaActive &&
	    std::find(this->d->lodSubmissionDeltaSources.begin(),
		this->d->lodSubmissionDeltaSources.end(), source) ==
		this->d->lodSubmissionDeltaSources.end()) {
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
		this->d->lodInteractive != FALSE) ? TRUE : FALSE;
	const SbBool poseOnlyFullyResident =
	    this->d->lodInteractive &&
	    !scaleInteraction &&
	    !this->d->lodSubmissionDeltaActive &&
	    sourceMeshRequests > 0 && viewLodState &&
	    viewLodState->cadMeshPayloadCountForSource(source) >=
		sourceMeshRequests;
	if (poseOnlyFullyResident) {
	    this->d->setLodConvergenceCandidateCount(
		source,
		viewLodState->cadMeshPayloadCountForSource(source));
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
	/* Source/inventory discovery already supplied a leaf proxy.  Complete the
	 * useful-coverage census without launching one PoP build per visible leaf;
	 * a fresh budgeted quality pass follows immediately.  A zoom census keeps
	 * existing cuts and is therefore not structural-only. */
	action.setStructuralCoverageOnly(
	    this->d->lodCoveragePolicy.active() &&
	    !this->d->lodViewDemandPolicy.scaleDemandRefreshActive());
	action.setStructuralPresentationRepair(
	    this->d->lodStructuralPresentationRepairPending);
	if (this->d->lodStructuralPresentationRepairPending)
	    action.setStructuralCoverageCostReservation(
		this->d->lodStructuralCoverageCostReservation);
	action.setSelectedOccurrenceCount(selectedOccurrenceCount);
	action.setGeneration(generation);
	action.setRevisions(this->d->lodViewRevision.value(),
	    this->d->lodPolicyRevision.value());
	action.setRefreshMissing(refreshMissing);
	action.setReset(reset);
	action.setViewLodState(viewLodState);
	action.setProjectedDemandCache(&this->d->lodProjectedDemandCache);
	action.setRefinementCutCeiling(
	    this->d->lodInteractive ?
		this->d->lodInteractiveProgressiveCeiling : -1);
	/* Coarsening a retained prefix only changes the draw count/snap level;
	 * it does not rebuild or reread the mesh.  Permit it during motion so a
	 * previously settled, expensive cut cannot pin interactive FPS. */
	const SbBool scaleDemandChanged =
	    (!this->d->lodInteractive &&
	     !this->d->lodRetainPoseOccurrenceCuts) || scaleInteraction;
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
	    !this->d->lodInteractive &&
	    !this->d->lodGestureActive &&
	    !this->d->lodResidentGrowthResidencyDrainActive &&
	    (scaleDemandChanged || applyRetainedAdmission));
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
	     !this->d->lodInteractive || scaleDemandChanged) &&
	    !this->d->lodResidentGrowthResidencyDrainActive &&
	    !this->d->lodFrameObligation.pending() ? TRUE : FALSE);
	/* A Schmitt boundary is useful while successive input events perturb the
	 * projection by fractions of a pixel.  It must not redefine the quiet
	 * view's convergence target: the first post-gesture pass recomputes the
	 * exact producer cut and either presents it or records genuine budget /
	 * residency debt. */
	action.setCutHysteresisEnabled(
	    this->d->lodInteractive || this->d->lodGestureActive);
	/* Once coverage is proven, residency follows pixel demand independently
	 * of the calibrated draw allowance.  This applies during zoom and during
	 * quiet warm-cache convergence: one worker read may append the complete
	 * demanded immutable suffix while the currently affordable cut remains on
	 * screen.  Working-set and resident-byte admission remain authoritative. */
    const SbBool quietCoveredPrefetch =
	!this->d->lodInteractive &&
	this->d->lodCoveragePolicy.coverageComplete() ? TRUE : FALSE;
    action.setAllowResidentPrefetch(
	    (scaleInteraction || quietCoveredPrefetch) &&
	    !this->d->lodFrameObligation.pending() ? TRUE : FALSE);
    action.setAllowResidentPrefetchPastAllocation(scaleInteraction);
	action.setAllowRepresentationRefinement(
	    scaleDemandChanged &&
	    !this->d->lodResidentGrowthResidencyDrainActive &&
	    !this->d->lodFrameObligation.pending() ? TRUE : FALSE);
	/* A calibrated face budget governs richness above the minimum drawable
	 * mesh floor.  Returning a visible resident PoP payload to its box is
	 * both a visual regression and a false economy: the renderer can batch
	 * tiny minimum prefixes into the aggregate point channel without
	 * discarding the useful retained asset. */
	action.setPreserveMeshCoverage(TRUE);
	action.setRefinementCostBudget(
	    this->d->lodBudgetPolicy.refinementRemaining());
	if (this->d->lodBudgetPolicy.singleOccurrenceBootstrap())
	    action.setInitialProviderCostBudget(
		this->d->lodBudgetPolicy.refinementRemaining());
	if (this->d->lodDiscretePopulationTrialAvailable) {
	    size_t discreteAllowance =
		BObolLodQualityPolicy::discreteTrialOverBudgetAllowance(
		    sceneActiveCost,
		    this->d->lodBudgetPolicy.currentBudget());
	    action.setOneOverBudgetRefinementLimit(discreteAllowance);
	}
	/* Every bounded window owns a disjoint part of the pinned occurrence
	 * plan and consumes from this carried scene-wide remainder.  The action
	 * deliberately skips its full priority recovery when a finite window is
	 * configured, avoiding an O(scene size) input stall. */
	if (applyRetainedAdmission)
	    action.setRetainedSceneUpgradeCostBudget(
		this->d->lodBudgetPolicy.retainedAdmissionRemaining());
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
	    (this->d->lodInteractive ? 3000 : 8000);
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
	const size_t compactWave = this->d->lodInteractive ?
	    interactiveCompactWave : quietCompactWave;
	const int compactCount = source->getCompactInstanceCount();
	const SbBool boundedLargeCompact =
	    compactCount > static_cast<int>(compactWave) ?
	    TRUE : FALSE;
	if (boundedLargeCompact)
	    boundedScenePass = TRUE;
	const bool usingMemoryAdmissionFrontier =
	    boundedLargeCompact &&
	    this->d->lodResidentAdmissionRetryActive &&
	    !this->d->lodSubmissionDeltaActive && viewLodState;
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
	    !this->d->lodSubmissionDeltaActive &&
	    !this->d->lodStructuralPresentationRepairPending &&
	    !this->d->lodCoveragePolicy.active() &&
	    this->d->lodCoveragePolicy.coverageComplete() &&
	    !this->d->lodBudgetPolicy.retainedAdmission() &&
	    viewLodState && sourceMeshRequests > 0 &&
	    (usingMemoryAdmissionFrontier ||
	     currentViewSourceCensusComplete);
	bool selectiveDeltaPlan = false;
	if (this->d->lodSubmissionDeltaActive &&
	    (!this->d->lodSubmissionPlanValid ||
	     this->d->lodSubmissionPlanSource != source)) {
	    for (const auto &deltaPlan :
		this->d->lodSubmissionDeltaPlans) {
		if (deltaPlan.first != source)
		    continue;
		selectiveDeltaPlan = true;
		this->d->lodSubmissionPlanEntries = deltaPlan.second;
		this->d->lodSubmissionPlanSource = source;
		this->d->lodSubmissionPlanValid = TRUE;
		break;
	    }
	} else if (this->d->lodSubmissionDeltaActive) {
	    selectiveDeltaPlan = std::any_of(
		this->d->lodSubmissionDeltaPlans.begin(),
		this->d->lodSubmissionDeltaPlans.end(),
		[source](const auto &deltaPlan) {
		    return deltaPlan.first == source;
		});
	}
	if (boundedLargeCompact && applyRetainedAdmission &&
	    (!this->d->lodSubmissionPlanValid ||
	     this->d->lodSubmissionPlanSource != source)) {
	    controller_prioritize_lod_recovery(source, viewLodState,
		this->d->lodSubmissionPlanEntries,
		selectedOccurrenceCount == 1);
	    this->d->lodSubmissionPlanSource = source;
	    this->d->lodSubmissionPlanValid = TRUE;
	    this->d->lodSubmissionPlanRetainedAdmission = TRUE;
	}
	if (usingUnsatisfiedFrontier &&
	    (!this->d->lodSubmissionPlanValid ||
	     this->d->lodSubmissionPlanSource != source)) {
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
	    controller_prioritize_lod_frontier(source, viewLodState,
		unsatisfiedKeys, this->d->lodSubmissionPlanEntries,
		selectedOccurrenceCount == 1);
	    /* Unsatisfied residency is deliberately retained while an occurrence
	     * leaves the frustum so a later camera pose can reuse its immutable
	     * mesh.  It is not actionable work in the current view, however.  The
	     * exact dense census is already the authority which enabled this sparse
	     * frontier; intersect with it here.  Otherwise an off-screen occurrence
	     * is projected, skipped, and reconstructed as the same one-entry plan
	     * on every GUI pump, preventing a quiet view from ever converging. */
	    this->d->lodSubmissionPlanEntries.erase(
		std::remove_if(this->d->lodSubmissionPlanEntries.begin(),
		    this->d->lodSubmissionPlanEntries.end(),
		    [this, source](size_t entryIndex) {
			return !this->d->isVisibleLodConvergenceCandidate(
			    source, entryIndex);
		    }),
		this->d->lodSubmissionPlanEntries.end());
	    this->d->lodSubmissionPlanSource = source;
	    this->d->lodSubmissionPlanValid = TRUE;
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
	    (!this->d->lodSubmissionPlanValid ||
	     this->d->lodSubmissionPlanSource != source)) {
	    this->d->lodSubmissionPlanEntries.resize(
		static_cast<size_t>(compactCount));
	    for (size_t planIndex = 0;
		 planIndex < this->d->lodSubmissionPlanEntries.size();
		 planIndex++)
		this->d->lodSubmissionPlanEntries[planIndex] = planIndex;
	    this->d->lodSubmissionPlanSource = source;
	    this->d->lodSubmissionPlanValid = TRUE;
	} else if (boundedLargeCompact && !usingUnsatisfiedFrontier &&
	    (!this->d->lodSubmissionDeltaActive || !selectiveDeltaPlan) &&
	    this->d->lodSubmissionPlanValid &&
	    this->d->lodSubmissionPlanSource == source &&
	    this->d->lodSubmissionPlanEntries.size() <
		static_cast<size_t>(compactCount)) {
	    /* Structural streaming appends leaves while a pinned scan is in
	     * flight.  Extend the plan tail without restarting the consumed
	     * prefix; restarting at zero on every batch starves late leaves. */
	    const size_t oldCount =
		this->d->lodSubmissionPlanEntries.size();
	    this->d->lodSubmissionPlanEntries.resize(
		static_cast<size_t>(compactCount));
	    for (size_t planIndex = oldCount;
		 planIndex < this->d->lodSubmissionPlanEntries.size();
		 planIndex++)
		this->d->lodSubmissionPlanEntries[planIndex] = planIndex;
	}
	if (this->d->lodSubmissionPlanValid &&
	    this->d->lodSubmissionPlanSource == source)
	    action.setCompactEntryPlanView(
		&this->d->lodSubmissionPlanEntries);
	action.setCompactEntryRange(this->d->lodSubmissionEntryOffset,
	    boundedLargeCompact ? compactWave : SIZE_MAX);
	if (this->d->lodUseForcedCut)
	    action.setForcedCut(this->d->lodForcedCut);
	const bool fullViewCandidatePass =
	    !usingUnsatisfiedFrontier &&
	    !this->d->lodSubmissionDeltaActive;
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
		this->d->lodSubmissionPlanEntries.size(),
		action.getVisitedMeshCount(), action.getSubmittedTaskCount(),
		action.getUpdatedCutCount(),
		action.hasDeferredCompactEntries() ? 1 : 0,
		action.getCompactEntryNext(),
		static_cast<unsigned long long>(
		    submissionTimeLimitMicroseconds),
		this->d->forceTerminalLodRefinement ? 1 : 0);
	if (action.getOneOverBudgetRefinementUsed())
	    this->d->lodDiscretePopulationTrialAvailable = FALSE;
	if (source->hasCompactInstanceIndex())
	    this->d->observeLodConvergenceCandidateVisibility(source,
		static_cast<size_t>(std::max(0, compactCount)),
		action.getCompactEntryVisibilityObservations());
	if (fullViewCandidatePass) {
	    if (this->d->lodSubmissionVisibleCountSource != source) {
		this->d->lodSubmissionVisibleCountSource = source;
		this->d->lodSubmissionVisibleCount = 0;
	    }
	    const size_t visible = action.getVisibleMeshCount();
	    this->d->lodSubmissionVisibleCount =
		visible > SIZE_MAX -
			this->d->lodSubmissionVisibleCount ?
		    SIZE_MAX :
		    this->d->lodSubmissionVisibleCount + visible;
	}
	if (this->d->lodCoveragePolicy.active()) {
	    const size_t visible = action.getVisibleMeshCount();
	    const size_t covered = action.getCoveredVisibleMeshCount();
	    this->d->lodCoveragePolicy.observe(visible, covered);
	}
	const bool sparseFrontierTrace =
	    getenv("BOBOL_LOD_TRACE_SPARSE_FRONTIER") &&
	    this->d->lodSubmissionPlanValid &&
	    this->d->lodSubmissionPlanSource == source &&
	    this->d->lodSubmissionPlanEntries.size() <= 512;
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
		       this->d->lodSubmissionPlanEntries.size(),
		       action.getDiagnostics().getString());
	}
	if (fullViewCandidatePass &&
	    !action.hasDeferredCompactEntries()) {
	    if (source->hasCompactInstanceIndex())
		this->d->completeLodConvergenceCandidateCensus(source,
		    this->d->lodSubmissionVisibleCount);
	    else
		this->d->setLodConvergenceCandidateCount(source,
		    this->d->lodSubmissionVisibleCount);
	}
	if (!this->d->lodSubmissionPlanValid &&
	    source->hasCompactInstanceIndex()) {
	    action.getCompactEntryPlan(
		this->d->lodSubmissionPlanEntries);
	    this->d->lodSubmissionPlanSource = source;
	    this->d->lodSubmissionPlanValid = TRUE;
	}

	this->d->lastLodVisitedMeshCount += action.getVisitedMeshCount();
	this->d->lastLodSubmittedTaskCount += action.getSubmittedTaskCount();
	this->d->lastLodUpdatedCutCount += action.getUpdatedCutCount();
	if (action.getSubmittedTaskCount() > 0 ||
	    action.getUpdatedCutCount() > 0)
	    this->d->lodPassAdmittedWork = TRUE;
	/* Accumulate every resident cut mutation across the complete bounded
	 * pass.  The final barrier below publishes one coherent population; no
	 * individual 3 ms planning slice is entitled to rebuild the whole CAD
	 * command record. */
	if (action.getUpdatedCutCount() > 0 &&
	    (!this->d->lodBudgetPolicy.retainedAdmission() ||
	     boundedLargeCompact))
	    this->d->lodRetainedRefinementCutAdvanced = TRUE;
	if (action.getPendingRetainedRefinementCount() > 0)
	    this->d->lodRetainedRefinementPending = TRUE;
	if (action.getPendingResidentRefinementCount() > 0)
	    this->d->lodRetainedResidencyPending = TRUE;
	if (action.getRefinementBudgetBlockedCount() > 0) {
	    /* A minimax quality ceiling is a terminal allocation for the current
	     * allowance, but it is still the witness which says that allowance
	     * prevented pixel-target convergence.  Route both an unreachable
	     * allocated cut and a deliberately coarser minimax cut through the
	     * same bounded unchanged-frame calibration.  BObolLodBudgetPolicy's
	     * three-sample series and the one-shot headroom witness make a truly
	     * saturated population terminal; suppressing the pure quality case
	     * here instead stranded conservative cold seeds forever (most visibly
	     * a one-leaf Lucy view at its first few thousand faces). */
	    this->d->lodRetainedRefinementBudgetBlocked = TRUE;
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
		   this->d->lodBudgetPolicy.refinementRemaining());
	this->d->lodBudgetPolicy.consumeRefinement(
	    action.getRefinementCostBudgetUsed());
	if (applyRetainedAdmission) {
	    this->d->lodBudgetPolicy.consumeRetainedAdmission(
		action.getRetainedSceneUpgradeCostBudgetUsed());
	}
	this->d->lastLodSkippedMeshCount += action.getSkippedMeshCount();
	const size_t missingMeshBlocked =
	    action.getMissingMeshBudgetBlockedCount();
	this->d->lodMissingMeshBudgetBlockedCount =
	    missingMeshBlocked >
		    SIZE_MAX - this->d->lodMissingMeshBudgetBlockedCount ?
		SIZE_MAX : this->d->lodMissingMeshBudgetBlockedCount +
		    missingMeshBlocked;
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
	this->d->lodResidentAdmissionRetryActive = FALSE;
    const bool completedStructuralRepair =
	completedPass && this->d->lodStructuralPresentationRepairPending;
    const size_t completedMissingMeshBudgetBlockedCount =
	completedPass ? this->d->lodMissingMeshBudgetBlockedCount : 0;
    if (completedStructuralRepair) {
	this->d->lodStructuralPresentationRepairPending = FALSE;
	this->d->lodStructuralRepairFrontierCount = 0;
	this->d->lodStructuralCoverageCostReservation = 0;
    }
    if (completedPass)
	this->d->lodMissingMeshBudgetBlockedCount = 0;
    const bool completedPassRetainedAllocation =
	completedPass &&
	this->d->lodSubmissionRetainedAdmissionMode;
    const bool completedPassChangedCut =
	completedPass &&
	this->d->lodRetainedRefinementCutAdvanced;
    const bool completedPassBudgetBlocked =
	completedPass &&
	this->d->lodRetainedRefinementBudgetBlocked;
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
		   this->d->lodSubmissionDeltaActive ? 1 : 0,
		   this->d->lodSubmissionRescanPending ? 1 : 0,
		   this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
		   this->d->lodRetainedRefinementPending ? 1 : 0,
		   this->d->lodRetainedResidencyPending ? 1 : 0,
		   completedStructuralRepair ? 1 : 0,
		   this->d->lodRetainedRefinementCutAdvanced ? 1 : 0,
		   this->d->lodBudgetPolicy.rescanAfterFrame() ? 1 : 0,
		   sceneActiveFaces, sceneActiveCost,
		   this->d->lodBudgetPolicy.currentBudget(),
		   this->d->lastLodVisitedMeshCount,
		   this->d->lastLodUpdatedCutCount);
    }
    const BObolLodCoveragePolicy::Completion coverageCompletion =
	this->d->lodCoveragePolicy.completeIfReady(
	    completedPass, this->d->lodSubmissionRescanPending);
    const bool retainedImportanceCensusCompleted =
	this->d->lodRetainedImportanceCensusPending &&
	coverageCompletion.completed && !coverageCompletion.missing;
    if (retainedImportanceCensusCompleted) {
	/* The pass which just finished recorded exact projected demand without
	 * discarding the useful retained population.  Reallocate that population
	 * once inside the already measured scene allowance.  This is neither an
	 * overload witness nor a cache/residency operation. */
	this->d->lodRetainedImportanceCensusPending = FALSE;
	this->d->lodBudgetPolicy.requestRetainedReallocation();
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD completed importance census view=%llu "
		   "policy=%llu active_cost=%zu budget=%zu\n",
		   (unsigned long long)this->d->lodViewRevision.value(),
		   (unsigned long long)this->d->lodPolicyRevision.value(),
		   sceneActiveCost,
		   this->d->lodBudgetPolicy.currentBudget());
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
    if (coverageCompletion.completed && coverageCompletion.bounded) {
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD completed coverage pass visible=%zu "
		   "covered=%zu missing=%d admitted=%d "
		   "budget_blocked=%d rescan_pending=%d\n",
		   coverageCompletion.visibleCount,
		   coverageCompletion.coveredCount,
		   coverageCompletion.missing ? 1 : 0,
		   this->d->lodPassAdmittedWork ? 1 : 0,
		   this->d->lodRetainedRefinementBudgetBlocked ? 1 : 0,
		   this->d->lodSubmissionRescanPending ? 1 : 0);
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodBudgetPolicy.clearBudgetLimit();
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	this->d->lodViewDemandPolicy.clearDemandRefresh();
	if (coverageCompletion.missing) {
	    /*
	     * Capacity for the next promotion wave must be learned from a
	     * presented frame, just like ordinary scene-budget admission.
	     * completeRenderTiming() restarts this coverage pass with the new
	     * allowance.  Until then, pending provider work may continue on the
	     * worker pool without turning the GUI loop into a busy rescan.
	     */
	    this->d->lodSubmissionPending = FALSE;
	    this->d->lodBudgetPolicy.requestRescanAfterFrame();
	    this->requestRender("lod-coverage-admission");
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
	    this->d->lodSubmissionPending = FALSE;
	    this->d->lodBudgetPolicy.requestRescanAfterFrame(true);
	    this->requestRender("lod-coverage-minimum-calibration");
	} else {
	    /* Every projected leaf has a useful structural or mesh presentation. Begin a
	     * fresh bounded pass which may now spend the remaining scene budget
	     * on screen-value-ordered PoP refinement. */
	    this->d->lodSubmissionPending = TRUE;
	    this->markProgressiveWorkPending();
	}
    } else if (BObolLodCoveragePolicy::needsDeferredQualitySuccessor(
	coverageCompletion.completed, coverageCompletion.missing,
	this->d->lodRetainedRefinementPending != FALSE,
	completedPassRetainedAllocation)) {
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
	 * A retained-allocation pass is already that authoritative successor.
	 * Routing its residual quality annotation through this branch starts an
	 * ordinary pass before the budget/handoff completion below can consume the
	 * allocation proof.  The ordinary pass then asks for another allocation,
	 * creating an allocation -> ordinary -> allocation loop over the complete
	 * scene.
	 */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = sources.empty() ? FALSE : TRUE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedResidencyPending = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	if (this->d->lodSubmissionPending)
	    this->markProgressiveWorkPending();
    } else if (retainedImportanceCensusCompleted) {
	/* Small sources do not take the bounded-coverage successor branch above,
	 * so explicitly start their one-shot importance allocation here. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = sources.empty() ? FALSE : TRUE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	if (this->d->lodSubmissionPending)
	    this->markProgressiveWorkPending();
    } else if (completedPass && this->d->lodSubmissionRescanPending) {
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
	if (this->d->lodSubmissionDeltaActive) {
	    this->d->lodSubmissionDeltaActive = FALSE;
	    this->d->lodSubmissionDeltaSources.clear();
	    this->d->lodSubmissionDeltaPlans.clear();
	}
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
	    this->d->lodSubmissionPending = FALSE;
	} else {
	    if (this->d->lodCoveragePolicy.active())
		this->d->lodCoveragePolicy.clearPassCounters();
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodSubmissionPending = sources.empty() ? FALSE : TRUE;
	}
	this->d->lodBudgetPolicy.clearBudgetLimit();
    } else if (completedPass && this->d->lodRetainedRefinementPending &&
	this->d->lodRetainedRefinementCutAdvanced &&
	!this->d->lodRetainedRefinementBudgetBlocked) {
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
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = FALSE;
	this->d->lodBudgetPolicy.clearBudgetLimit();
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodBudgetPolicy.resetProbeSeries();
    } else if (completedPass &&
	this->d->lodRetainedRefinementBudgetBlocked) {
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
	/* Stable pass decisions require a duration for this retained cut.  GPU
	 * query latency is handled by the paired throughput calibration above;
	 * reusing that duration after the occurrence population changes is not a
	 * current-cut deadline measurement. */
	const uint64_t observedStableNanoseconds =
	    this->d->lodLastRenderWasReusableCadPresentation ?
		(this->d->lastRenderTimeNanoseconds ?
		    this->d->lastRenderTimeNanoseconds :
		    this->d->smoothedRenderTimeNanoseconds) : 0;
	long double calibratedStableCostBudget = 0.0L;
	if (this->d->quietAllocationTargetFps() > 0.0f &&
	    this->d->lodStableCalibratedRenderCostPerSecond > 0.0L)
	    calibratedStableCostBudget =
		this->d->lodStableCalibratedRenderCostPerSecond /
		static_cast<long double>(
		    this->d->quietAllocationTargetFps());
	/* A probe is useful only while it can establish a larger affordable cut.
	 * Three unchanged samples reject one-time setup noise.  Capacity growth is
	 * then tested by the next distinct admitted population; repeating the same
	 * cut further cannot prove that a richer cut is affordable. */
	/* Never classify a newly reached quiet cut from a single frame.  OSMesa
	 * in particular can pay shader/cache/setup work in the first traversal,
	 * while a compositor stall can perturb System GL.  Three samples remain
	 * a small bounded cost at the stable target and let subsequent ordinary
	 * traversal cost govern the retained budget. */
	const unsigned int calibrationProbeCount =
	    this->d->lodBudgetPolicy.probeCount();
	BObolLodBudgetPolicy::CalibrationInputs calibrationInputs;
	calibrationInputs.activeCost = sceneActiveCost;
	calibrationInputs.targetNanoseconds = stableTargetNanoseconds;
	calibrationInputs.observedNanoseconds = observedStableNanoseconds;
	calibrationInputs.calibratedBudget = calibratedStableCostBudget;
	calibrationInputs.interactive = this->d->lodInteractive != FALSE;
	calibrationInputs.stablePresentationHandoff =
		    this->d->lodPresentationPolicy.handoffPending();
	calibrationInputs.passAdmittedWork =
	    this->d->lodPassAdmittedWork != FALSE;
	calibrationInputs.retainedAllocation =
	    completedPassRetainedAllocation;
	const BObolLodBudgetPolicy::CalibrationDecision calibration =
	    this->d->lodBudgetPolicy.finishBlockedPass(calibrationInputs);
	/* A structural repair is already ordered by projected importance.  If its
	 * finite scene allowance cannot replace the entire exact box frontier,
	 * preserve the admitted prominent meshes and aggregate the remaining
	 * small-object tail before another pass.  This avoids tens of thousands of
	 * individually affordable frames while retaining immutable resident data
	 * for later views and faster renderers. */
	BObolLodPointProxyCalibrationPolicy::Decision structuralAggregation;
	const size_t structuralTailBlockedCount =
	    completedStructuralRepair ?
		completedMissingMeshBudgetBlockedCount : 0;
	const float structuralPointThresholdBefore =
	    this->d->lodPresentationPointProxyPixelThreshold;
	if (structuralTailBlockedCount > 0 &&
	    !this->d->lodInteractive &&
	    this->d->pointProxyAggregationApplicable()) {
	    structuralAggregation = this->d->lodPointProxyCalibrationPolicy.
		structuralCoverageBlocked(
		    this->d->lodPresentationPointProxyPixelThreshold,
		    structuralTailBlockedCount);
	    if (structuralAggregation.changed) {
		this->d->lodPresentationPointProxyPixelThreshold =
		    structuralAggregation.threshold;
		BObolViewLodState *presentationState =
		    this->d->viewAttachment->getViewLodState();
		if (presentationState)
		    presentationState->
			setCadPresentationPointProxyPixelThreshold(
			    structuralAggregation.threshold);
		this->d->lodStablePointProxyCalibrationPending = TRUE;
		this->requestRender("lod-structural-tail-aggregation");
	    }
	}
	if (getenv("BOBOL_LOD_TRACE_STRUCTURAL_BUDGET") &&
	    structuralTailBlockedCount > 0)
	    bu_log("BObol LoD structural allocation tail=%zu "
		   "interactive=%d applicable=%d threshold=%.6g "
		   "next=%.6g changed=%d\n",
		   structuralTailBlockedCount,
		   this->d->lodInteractive ? 1 : 0,
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
		       "active_faces=%zu active_cost=%zu cost_budget=%zu admitted=%d "
		       "probe_candidate=%d probe_eligible=%d "
		       "probe_count=%u probe_next=%d "
		       "observed_ms=%.3f target_ms=%.3f "
		       "calibrated_mcost_s=%.3f\n",
		       sceneActiveFaces, sceneActiveCost,
		       this->d->lodBudgetPolicy.currentBudget(),
		       this->d->lodPassAdmittedWork ? 1 : 0,
		       calibration.probeCandidate ? 1 : 0,
		       calibration.probeEligible ? 1 : 0,
		       calibrationProbeCount,
		       calibration.calibrationFrame ? 1 : 0,
		       observedStableNanoseconds / 1000000.0,
		       stableTargetNanoseconds / 1000000.0,
		       static_cast<double>(
			   this->d->lodStableCalibratedRenderCostPerSecond /
			       1000000.0L));
	}
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = FALSE;
	if (calibration.requestFrame) {
	    const char *frameReason = calibration.calibrationFrame ?
		"lod-budget-calibration" : "lod-budget-admission";
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
	if (!structuralAggregation.changed &&
	    !this->d->lodBudgetPolicy.rescanAfterFrame() &&
	    this->d->lodBudgetPolicy.stableBudgetLimited() &&
	    !this->d->lodInteractive &&
	    this->d->pointProxyAggregationApplicable() &&
	    sceneLodState && sceneLodState->hasCadPresentationAssemblies() &&
	    sceneActiveCost <= (sceneLodState ?
		sceneLodState->minimumActiveRenderCost() : sceneActiveCost) &&
	    sceneActiveCost > this->d->lodBudgetPolicy.currentBudget() &&
	    stableTargetNanoseconds > 0 &&
	    observedStableNanoseconds >
		static_cast<long double>(stableTargetNanoseconds) * 1.20L) {
		const BObolLodPointProxyCalibrationPolicy::Decision pressure =
		    this->d->lodPointProxyCalibrationPolicy.interrupted(
			this->d->lodPresentationPointProxyPixelThreshold,
			static_cast<uint64_t>(observedStableNanoseconds),
			this->d->quietAllocationTargetFps());
	    if (pressure.changed) {
		this->d->lodPresentationPointProxyPixelThreshold =
		    pressure.threshold;
		BObolViewLodState *presentationState =
		    this->d->viewAttachment->getViewLodState();
		if (presentationState)
		    presentationState->
			setCadPresentationPointProxyPixelThreshold(
			    pressure.threshold);
		this->d->lodStablePointProxyCalibrationPending = TRUE;
		/* The multiplicative pressure correction intentionally lands on
		 * the safe side of the target.  Its next unchanged frame continues
		 * the bounded bracket search, so terminal quality is the finest cut
		 * which meets the stable FPS contract. */
		this->requestRender("lod-stable-point-calibration");
	    }
	}
	append_controller_lod_diagnostic(this->d->lastLodDiagnostics,
	    SbString(""),
	    this->d->lodStablePointProxyCalibrationPending ?
		"scene LoD calibrating small-part aggregation" :
	    this->d->lodBudgetPolicy.rescanAfterFrame() ?
		(calibration.calibrationFrame ?
		    "scene LoD probing stable calibrated capacity" :
		    "scene LoD admission awaiting calibrated frame") :
		"scene LoD reached its calibrated face budget");
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodBudgetPolicy.resetPass();
    } else {
	this->d->lodSubmissionPending = completedPass ? FALSE : TRUE;
	if (completedPass)
	    this->d->lodPassAdmittedWork = FALSE;
	if (completedPass)
	    this->d->lodBudgetPolicy.clearBudgetLimit();
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
    if (completedPass && this->d->lodBudgetPolicy.retainedAdmission() &&
	completedPassChangedCut) {
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = FALSE;
    }

    if (completedPass && this->d->lodRetainedResidencyPending &&
	!completedPassChangedCut && !this->d->lodSubmissionPending &&
	!this->d->lodFrameObligation.pending()) {
	/* The minimax allocation did not need to change the currently drawable
	 * prefix, but one or more allocated levels still need a resident suffix.
	 * There is no presentation barrier to wait for, so start the ordinary
	 * provider pass now.  Treating this as a mere quality observation made a
	 * quiet Hubble view report ready until an unrelated erase or selection
	 * event happened to wake the missing cache requests. */
	this->d->lodRetainedResidencyPending = FALSE;
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = sources.empty() ? FALSE : TRUE;
	this->d->lodBudgetPolicy.resetPass();
	if (this->d->lodSubmissionPending)
	    this->markProgressiveWorkPending();
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
    BObolLodPresentationPolicy::CompletedPassInputs handoffInputs;
    handoffInputs.completed = completedPass != FALSE;
    handoffInputs.submissionPending =
	this->d->lodSubmissionPending != FALSE;
    handoffInputs.rescanAfterFrame =
	this->d->lodBudgetPolicy.rescanAfterFrame();
    handoffInputs.changedCut = completedPassChangedCut != FALSE;

    /* A retained-allocation commit proves the owner-thread CAD revision it
     * observed, but it cannot predict a provider result which is still in
     * flight.  Releasing the global deadline ceiling while that stream was
     * growing made each newly installed batch invalidate the measured
     * population, producing an exact constrained -> release -> deadline
     * cycle on 50k warm-cache scenes.  Finish handoff only after both the
     * immutable-result stream and its presentation/coalescing edges are
     * quiet.  The resident-growth policy below supplies the eventual
     * scene-wide retry rather than busy-scanning while workers run. */
    const bool handoffServiceQuiescent = !service ||
	(service->activeTaskCountForGeneration(generation) == 0 &&
	 service->queuedResultCountForGeneration(generation) == 0 &&
	 this->d->lodResultsPending.load() == 0);
    const bool handoffPopulationQuiescent =
	handoffServiceQuiescent &&
	retainedPopulationSettled &&
	!this->d->lodFrameObligation.pending() &&
	!this->d->lodPublicationPolicy.pending() &&
	!this->d->lodResidentGrowthPolicy.pending();
    handoffInputs.retainedAllocationCompleted =
	completedPassRetainedAllocation;
    handoffInputs.populationQuiescent =
	handoffPopulationQuiescent;
    const size_t handoffReconciliationBudget =
	this->d->lodPresentationPolicy.handoffReconciliationBudget();
    const BObolRetainedAllocationResult &allocationCertificate =
	this->d->lodRetainedAllocationCertificate;
    const bool allocationCertificateCurrent = sceneLodState &&
	allocationCertificate.allocationPlanSerial != 0 &&
	allocationCertificate.allocationPlanSerial ==
	    sceneLodState->activeCadAllocationPlan() &&
	allocationCertificate.viewRevision ==
	    this->d->lodViewRevision.value() &&
	allocationCertificate.policyRevision ==
	    this->d->lodPolicyRevision.value() &&
	std::isfinite(allocationCertificate.pointProxyPixelThreshold) &&
	std::fabs(allocationCertificate.pointProxyPixelThreshold -
	    this->d->lodPresentationPointProxyPixelThreshold) <= 1.0e-6f &&
	sceneLodState->cadAllocationPlanCoversCurrentPopulation(
	    allocationCertificate.allocationPlanSerial,
	    allocationCertificate.viewRevision,
	    allocationCertificate.policyRevision,
	    allocationCertificate.fixedCadPresentationCost);
    handoffInputs.retainedAllocationCertified =
	allocationCertificateCurrent;
    handoffInputs.presentationLimitsReconciled =
	BObolLodPresentationPolicy::presentationLimitsReconciled(
	    completedPassRetainedAllocation,
	    allocationCertificateCurrent,
	    allocationCertificate.selectedPresentationCost,
	    allocationCertificate.certifiedPresentationBudget,
	    handoffReconciliationBudget);
    handoffInputs.retainedRefinementPending =
	this->d->lodRetainedRefinementPending != FALSE;
    handoffInputs.retainedRefinementBudgetBlocked =
	completedPassBudgetBlocked;
    const BObolLodPresentationPolicy::CompletedPassDecision handoff =
	this->d->lodPresentationPolicy.completePass(handoffInputs);
    if (completedPass && controller_lod_trace_enabled(
	    "BOBOL_LOD_TRACE_PASS", this->d->lodViewRevision.value()))
	bu_log("BObol LoD handoff pass pending=%d rescan=%d changed=%d "
	       "reconciled=%d "
	       "allocation=%d certified=%d selected=%zu budget=%zu "
	       "quiescent=%d "
	       "retained=%d blocked=%d presentation_pending=%d "
	       "finish=%d request_allocation=%d request_rescan=%d "
	       "retire=%d preserve=%d\n",
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
	       handoff.retireRetainedObservation ? 1 : 0,
	       handoff.preservePresentationLimits ? 1 : 0);

    const bool handoffResidentPopulationPending =
	!handoffServiceQuiescent ||
	this->d->lodResultsPending.load() != 0 ||
	this->d->lodPublicationPolicy.pending();
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
	this->d->lodResidentGrowthPolicy.noteRicherPrefixAvailable();
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
	    this->d->lodBudgetPolicy.requestPresentationReconciliation(
		reconciliationBudget);
	else
	    this->d->lodBudgetPolicy.requestRetainedReallocation();
	/* This is another pass in the same camera/policy/source epoch.  The
	 * explicit pending cursor is sufficient to bypass the completed-pass fast
	 * path.  Clearing the epoch witness as well makes the wrapper classify the
	 * already-pending cursor as a view change during submission and append an
	 * unnecessary full rescan after every allocation. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = sources.empty() ? FALSE : TRUE;
	/* These latches summarize one complete bounded transaction.  The
	 * allocation request below starts a new transaction in the same epoch;
	 * carrying a preceding cut-advanced bit made an unchanged pass report
	 * changed and prevented handoff completion indefinitely. */
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedResidencyPending = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	if (this->d->lodSubmissionPending)
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
	const BObolLodPointProxyCalibrationPolicy::Decision reduction =
	    this->d->lodPointProxyCalibrationPolicy.
		structuralCoverageBlocked(
		    this->d->lodPresentationPointProxyPixelThreshold,
		    unresolved);
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
	    if (retryBudget > 0)
		this->d->lodBudgetPolicy.requestPresentationReconciliation(
		    retryBudget);
	    else
		this->d->lodBudgetPolicy.requestRetainedReallocation();
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodSubmissionPending = sources.empty() ? FALSE : TRUE;
	    this->d->lodPassAdmittedWork = FALSE;
	    this->d->lodRetainedRefinementPending = FALSE;
	    this->d->lodRetainedResidencyPending = FALSE;
	    this->d->lodRetainedRefinementCutAdvanced = FALSE;
	    this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	    this->d->lodBudgetPolicy.resetPass();
	    if (this->d->lodSubmissionPending)
		this->markProgressiveWorkPending();
	} else {
	    /* Reaching this state means every aggregatable occurrence is already
	     * collapsed and the protected minimum itself exceeds the measured hard
	     * budget.  Keep the completed framebuffer responsive and wait for an
	     * explicit capacity/resource/view edge; silently reinstating a global
	     * terminal PoP ceiling would violate the quality contract. */
	    append_controller_lod_diagnostic(
		this->d->lastLodDiagnostics, SbString("retained-allocation"),
		"protected visible minimum exceeds the measured static frame "
		"budget; retaining the last completed frame");
	    this->d->lodPresentationPolicy.cancelHandoff();
	    this->d->lodStaticOverscanActive = FALSE;
	    this->d->lodStaticOverscanLeapAvailable = FALSE;
	    this->d->lodStaticOverscanRejected = TRUE;
	}
    }
    /* The policy distinguishes an unsatisfied quality observation from an
     * actual progress witness.  Retiring only the former prevents a
     * performance-limited terminal cut from becoming an endless compaction
     * defer/BACKGROUND loop.  A later result, camera/policy epoch, or capacity
     * edge starts a fresh pass from the retained cut. */
    if (handoff.retireRetainedObservation)
	this->d->lodRetainedRefinementPending = FALSE;
    if (handoff.finishHandoff) {
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	/* Retained occurrence cuts and their immutable PoP suffix are the
	 * continuity baseline across camera epochs.  The former normalization
	 * path rewrote every cut to its minimum merely because a renderer-only
	 * ceiling had protected one handoff frame.  It then walked the same
	 * resident data upward again, producing the conspicuous coarse flash after
	 * every zoom and invalidating otherwise reusable command records.
	 *
	 * Remove only the reversible ceiling and measure the existing cut when the
	 * completed pass actually reconciled the retained occurrences.  A hard
	 * deadline handoff deliberately sets preservePresentationLimits: its first
	 * completed constrained frame is proof that the richer cut was too costly,
	 * not permission to retry that cut immediately.  Clearing the ceiling in
	 * that case produced an exact cycle (rich abort -> coarse completion -> rich
	 * abort), visible as alternating boxes/coarse meshes and a progress HUD that
	 * never settled.
	 *
	 * The endpoint render remains deadline-bounded and double buffered.  An
	 * explicit later capacity/view edge can relax a preserved presentation cut;
	 * memory-pressure compaction remains the sole authority for shrinking the
	 * retained occurrence prefixes. */
	if (!handoff.preservePresentationLimits) {
	    this->d->lodInteractiveProgressiveCeiling = -1;
	    if (presentationState)
		presentationState->setCadPresentationProgressiveCutCeiling(-1);
	}
	/* The full static overscan phase is complete.  Keep its rejection witness:
	 * a later bounded marginal allocation may still use the static deadline,
	 * but must not reopen the rejected all-or-nothing protected-floor trial. */
	if (this->d->lodStaticOverscanActive &&
	    this->d->lodStaticOverscanRejected) {
	    this->d->lodStaticOverscanActive = FALSE;
	    this->d->lodStaticOverscanLeapAvailable = FALSE;
	}
	if (presentationState)
	    presentationState->setCadPresentationCameraMotionFrameReuse(FALSE);
	this->d->lodStablePointProxyCalibrationPending =
	    presentationState &&
	    presentationState->hasCadPresentationAssemblies() &&
	    this->d->lodPresentationPointProxyPixelThreshold > 1.01f ?
		TRUE : FALSE;
	/* The motion-time renderer ceiling may have prevented an otherwise
	 * resident occurrence from advancing at all.  That leaves a real richer
	 * demand but no changed cut, worker, or result which could schedule the
	 * next pass.  Preserve the bounded first ceiling-free presentation, then
	 * turn that demand into the existing post-frame rescan witness.  Leaving
	 * the raw retained-pending annotation armed here creates an unwitnessed
	 * progressive/compaction loop; submitting immediately would skip the
	 * handoff frame and could expose an unbounded population. */
	if (handoff.requestRetainedRescan) {
	    this->d->lodRetainedRefinementPending = FALSE;
	    this->d->lodBudgetPolicy.requestRescanAfterFrame(true);
	}
	this->requestRender("lod-stable-handoff");
    }

    /* A retained-prefix recovery can discover that the current occurrence
     * population is already at its minimum cut.  No geometry changes in that
     * case, so no presentation barrier is requested and completeRenderTiming()
     * will not be called merely to retire this planning latch.  Complete the
     * no-op recovery here; recoveries which did change a cut still pass
     * through the render barrier below and are retired by the measured-frame
     * path. */
    if (this->d->lodPointProxyTriangleRecoveryPending && completedPass &&
	!completedPassChangedCut && !this->d->lodSubmissionPending &&
	!this->d->lodSubmissionRescanPending &&
	!this->d->lodRetainedRefinementPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodPresentationPolicy.handoffPending()) {
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	const BObolRetainedAllocationResult &recoveryAllocation =
	    this->d->lodRetainedAllocationCertificate;
	if (presentationState && recoveryAllocation.allocationPlanSerial != 0 &&
	    recoveryAllocation.viewRevision ==
		this->d->lodViewRevision.value() &&
	    recoveryAllocation.policyRevision ==
		this->d->lodPolicyRevision.value() &&
	    presentationState->cadAllocationPlanCoversCurrentPopulation(
		recoveryAllocation.allocationPlanSerial,
		recoveryAllocation.viewRevision,
		recoveryAllocation.policyRevision,
		recoveryAllocation.fixedCadPresentationCost))
	    this->d->lodPointProxyTriangleRecoverySaturatedPlanSerial =
		recoveryAllocation.allocationPlanSerial;
	this->d->lodPointProxyTriangleRecoveryPending = FALSE;
	BObolViewLodState::CadStructuralProjectionHistogram
	    structuralProjection;
	const SbBool exactStructuralProjection = presentationState &&
	    presentationState->lastCadStructuralProjectionHistogram(
		structuralProjection);
	const float recoveredPointThreshold =
	    BObolLodPointProxyCalibrationPolicy::
		triangleRecoveryPointThreshold(
		    this->d->lodPresentationPointProxyPixelThreshold,
		    exactStructuralProjection != FALSE,
		    exactStructuralProjection ?
			structuralProjection.visibleCount : SIZE_MAX);
	const SbBool pointCutChanged = std::fabs(
	    recoveredPointThreshold -
	    this->d->lodPresentationPointProxyPixelThreshold) > 0.01f ?
		TRUE : FALSE;
	this->d->lodInteractiveProgressiveCeiling = -1;
	this->d->lodPresentationPointProxyPixelThreshold =
	    recoveredPointThreshold;
	this->d->lodPointProxyCalibrationPolicy.reset();
	if (presentationState) {
	    presentationState->setCadPresentationProgressiveCutCeiling(-1);
	    presentationState->setCadPresentationPointProxyPixelThreshold(
		recoveredPointThreshold);
	}
	this->d->lodStablePointProxyCalibrationPending =
	    pointCutChanged && presentationState &&
	    presentationState->hasCadPresentationAssemblies() ? TRUE : FALSE;
	if (this->d->lodStablePointProxyCalibrationPending) {
	    this->requestRender("lod-stable-point-restore");
	} else {
	    this->resumeLodAfterRetainedRecovery();
	}
    }
    if (!this->d->lodSubmissionPending &&
	this->d->lodRetainedRefinementCutAdvanced &&
	(boundedScenePass || this->d->lodRetainedRefinementPending ||
	 this->d->lodPresentationPolicy.handoffPending())) {
	/* Richer prefixes selected by the completed pass are already in memory;
	 * expose that coherent budgeted wave in one completed frame.  This
	 * presentation barrier is also required when
	 * the just-completed pass reached every requested cut: in that case
	 * lodRetainedRefinementPending is false, but the motion-to-stable
	 * handoff still needs one frame with the newly selected cuts before a
	 * second unchanged pass may remove its renderer-wide ceiling.  Gating
	 * this on "pending" stranded the handoff forever at large-scene scale.
	 * Retargeting is metadata/draw-count only; no provider task, cache read,
	 * or geometry rebuild is involved. */
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	this->scheduleLodRefinementFrame("lod-cut");
    }
    /* The bounded calibration series may end before its first unchanged,
     * reusable CAD replay.  Request one explicit terminal probe now, while
     * this pass still owns the progress edge.  Its exact view/policy/budget
     * witness is one-shot; no later selection or HUD repaint may substitute
     * for this frame or reopen a view which already reported stable. */
    if (completedPass)
	this->armStableLodHeadroomProbeIfReady();
    if (completedPass && !this->d->lodSubmissionPending &&
	!this->d->lodRetainedRefinementPending) {
	this->d->lodBudgetPolicy.resetPass();
    }
    if (completedPass && !this->d->lodSubmissionPending)
	this->d->lodDiscretePopulationTrialAvailable = FALSE;
    if (completedPass)
	this->d->lodSubmissionRetainedAdmissionMode = FALSE;
    if (this->d->lodSubmissionPending)
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
    if (completedPass && this->d->lodSubmissionDeltaActive &&
	!this->d->lodSubmissionPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodFrameObligation.pending()) {
	this->d->lodSubmissionDeltaActive = FALSE;
	this->d->lodSubmissionDeltaSources.clear();
	this->d->lodSubmissionDeltaPlans.clear();
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
    this->d->lodResultsPending.store(
	service->queuedResultCountForGeneration(
	    drainGeneration) > 0 ? 1 : 0);
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
	if (extendsCumulativeCadAsset &&
	    drained[i].geometry.activeCut < residentCadBefore->activeCut &&
	    (residentCadBefore->presentedChunks.empty() ?
		drained[i].progressiveMesh->canDrawCut(
		    residentCadBefore->activeCut) :
		drained[i].progressiveMesh->canDrawChunksAtCut(
		    residentCadBefore->presentedChunks,
		    residentCadBefore->activeCut))) {
	    drained[i].geometry.activeCut = residentCadBefore->activeCut;
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
	    const bool residentSupersedes =
		(cadMesh &&
		 (cadPayload->policyRevision >
		      drained[i].request.policyRevision ||
		  (cadPayload->policyRevision ==
		       drained[i].request.policyRevision &&
		   cadPayload->viewRevision >=
		       drained[i].request.viewRevision))) ||
		(shapeMesh &&
		 (meshPayload->policyRevision >
		      drained[i].request.policyRevision ||
		  (meshPayload->policyRevision ==
		       drained[i].request.policyRevision &&
		   meshPayload->viewRevision >=
		       drained[i].request.viewRevision)));
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
	    BObolLodResidentGrowthPolicy::canRetainPresentation(
		this->d->lodResidentGrowthResidencyDrainActive != FALSE,
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
    this->d->lodResultsPending.store(
	service->queuedResultCountForGeneration(
	    this->d->lodActiveGeneration) > 0 ? 1 : 0);

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
	    this->d->lodResidentGrowthPolicy.noteRicherPrefixAvailable();
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		    this->d->lodViewRevision.value()))
		bu_log("BObol LoD resident growth noted applied=%zu "
		       "interactive=%d scale_changed=%d\n",
		       applied, this->d->lodInteractive ? 1 : 0,
		       this->d->lodViewDemandPolicy.interactionScaleChanged() ?
			   1 : 0);
	    /* The missing suffix may complete after this scale epoch has already
	     * spent its one quality probe.  Residency alone is not a visible
	     * improvement, so rearm exactly one bounded probe for the newly
	     * available population.  An ordinary interactive/first-useful result
	     * publishes below; a quiet residency-only wave retains the coherent
	     * frame until its aggregate reallocation consumes this latch. */
	    (void)this->d->lodViewDemandPolicy.rearmAfterResidentGrowth(
		this->d->lodInteractive != FALSE);
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
		this->d->lodResultsFirstReadyMicroseconds.load();
	    this->d->lodPublicationPolicy.noteApplied(
		applied, firstReady > 0 ? firstReady : now);
	}
	const bool firstUsefulMesh = !retainedPresentationBatch &&
	    this->getActiveLodMeshPayloadCount() <= applied;
	const bool serviceProducer =
	    service->queuedResultCountForGeneration(
		this->d->lodActiveGeneration) > 0 ||
	    service->activeTaskCountForGeneration(
		this->d->lodActiveGeneration) > 0;
    const bool submissionPausedByPresentation =
	BObolLodPointProxyCalibrationPolicy::blocksSourceAdmission(
	    this->d->lodAdmissionPointProxyFramePending != FALSE,
	    this->d->lodStablePointProxyCalibrationPending != FALSE);
	const bool streamIdle =
	    !BObolLodProducerPolicy::canProduceGeometry(
		this->d->lodSubmissionPending != FALSE,
		submissionPausedByPresentation,
		this->d->progressiveProviderPendingCount > 0,
		serviceProducer);
	BObolLodPublicationPolicy::Inputs publicationInputs;
	publicationInputs.nowMicroseconds = now;
	publicationInputs.observedRenderNanoseconds = std::max(
	    this->d->lastSceneRenderTimeNanoseconds,
	    this->d->lastRenderTimeNanoseconds);
	publicationInputs.interactive = this->d->lodInteractive != FALSE;
	publicationInputs.firstUseful = firstUsefulMesh;
	publicationInputs.streamIdle = streamIdle;
	const BObolLodPublicationPolicy::Decision publicationDecision =
	    this->d->lodPublicationPolicy.decide(publicationInputs);
	const bool publishNow =
	    this->d->lodPublicationPolicy.framePending();
	if (partialRefinementCandidate &&
	    !this->d->lodFrameObligation.pending()) {
	    /*
	     * Hold the next prefix immediately, but allow independent coverage
	     * results to accumulate until the adaptive publication deadline.
	     * scheduleLodRefinementFrame() below supplies the host wakeup once
	     * that deadline is reached.
	     */
	    (void)this->d->lodFrameObligation.arm(
		BObolLodFrameObligation::REASON_RESULT_PUBLICATION,
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
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = TRUE;
	this->d->lodRetainedRefinementPending = TRUE;
	this->d->lodPassAdmittedWork = FALSE;
	this->markProgressiveWorkPending();
	this->requestRender("lod-current-demand-replay");
    }

    if (this->d->lastLodResultCount > 0)
	this->d->reconcilePhase(
	    BObolLodStateMachine::Event::RESULT_PUBLISHED);

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
    this->d->lodPolicyRevision.set(newRevision);
    this->d->lodStaticOverscanActive = FALSE;
    this->d->lodStaticOverscanRejected = FALSE;
    struct bv_lod_policy lodPolicy;
    bv_lod_policy_init(&lodPolicy);
    this->d->viewAttachment->getLodPolicy(&lodPolicy);
    this->d->lodCoveragePolicy.setRequired(
	this->d->lodAutoSubmit && this->d->lodService &&
	lodPolicy.policy != BV_LOD_OFF && lodPolicy.mesh_enabled);
    this->d->lodHeadroomPolicy.cancelRetry();
    /*
     * A quality-policy revision changes the requested PoP cut, not camera
     * visibility.  Retain the last complete current-view denominator until a
     * camera or source-inventory revision invalidates it.  Clearing it here
     * let the quiet transition report a 50k scene converged with a zero
     * visible target while its unsatisfied-frontier fast path deliberately
     * skipped another full visibility scan.
     */
    this->d->resetLodConvergenceFraction();
    this->d->lodBudgetPolicy.clearBudgetLimit();
    this->d->lodBudgetPolicy.clearRetainedRecoveryCeiling();
    this->d->lodBudgetPolicy.clearDeadlineCapacityCeiling();
    this->d->lodBudgetPolicy.clearRetainedQualityFloorBudget();
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodPointProxyTriangleRecoveryPending = FALSE;
    this->d->lodPointProxyCalibrationPolicy.reset();
    this->requestRender("lod-policy");
    this->d->reconcilePhase(BObolLodStateMachine::Event::POLICY_CHANGED);
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
    this->d->lodViewRevision.advance();
    this->d->lodInterruptedPresentationReplayPending = FALSE;
    this->d->resetRetainedAdmissionQualityProof();
    this->d->lodStaticOverscanActive = FALSE;
    this->d->lodStaticOverscanRejected = FALSE;
    this->d->lodHeadroomPolicy.cancelRetry();
    this->d->clearLodConvergenceCandidates();
    this->d->resetLodConvergenceFraction();
    this->d->lodBudgetPolicy.clearBudgetLimit();
    this->d->lodBudgetPolicy.clearRetainedQualityFloorBudget();
    this->d->lodBudgetPolicy.clearRetainedRecoveryCeiling();
    this->d->lodBudgetPolicy.clearDeadlineCapacityCeiling();
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodPointProxyTriangleRecoveryPending = FALSE;
    this->d->lodPointProxyCalibrationPolicy.reset();
    this->d->lodPresentationPolicy.viewInvalidated();
    this->d->lodRetainPoseOccurrenceCuts = FALSE;
    this->d->lodRetainedImportanceCensusPending = FALSE;
    this->d->resetDeadlineSafePresentation();
    this->d->lodCoveragePolicy.activate(true);
    this->d->lodViewDemandPolicy.refreshForViewRevision(
	this->d->lodInteractive != FALSE);
    this->d->reconcilePhase(BObolLodStateMachine::Event::VIEW_INVALIDATED);
}

void
BObolViewController::advanceLodPolicyRevision(
    SbBool preserveScaleDemandRefresh)
{
    this->d->lodPolicyRevision.advance();
    this->d->lodInterruptedPresentationReplayPending = FALSE;
    this->d->resetRetainedAdmissionQualityProof();
    /* Every quality-policy epoch first converges at the preferred quiet
     * cadence.  Static overscan is armed later from an exact terminal frame
     * and does not itself mutate the semantic policy epoch. */
    this->d->lodStaticOverscanActive = FALSE;
    struct bv_lod_policy lodPolicy;
    bv_lod_policy_init(&lodPolicy);
    this->d->viewAttachment->getLodPolicy(&lodPolicy);
    this->d->lodCoveragePolicy.setRequired(
	this->d->lodAutoSubmit && this->d->lodService &&
	lodPolicy.policy != BV_LOD_OFF && lodPolicy.mesh_enabled);
    this->d->lodHeadroomPolicy.cancelRetry();
    /* A policy epoch may change the target cadence (most importantly the
	 * 60 Hz motion -> 20 Hz quiet transition), refinement authority, and
	 * presentation handoff.  A partially consumed pass is meaningful only for
	 * the inputs which initialized it.  Carrying its retained-admission flag
	 * into the new epoch made a 0.03% budget rounding difference normalize all
	 * 2,500 Hubble occurrences to their minimum prefixes after rotation.  Keep
	 * the calibrated total budget, but recompute the pass decision from the new
	 * epoch's exact inputs. */
    this->d->lodBudgetPolicy.resetPass();
    /* Quality changes preserve the current view's proven visibility
     * denominator.  Source and camera revisions clear it explicitly. */
    this->d->resetLodConvergenceFraction();
    this->d->lodBudgetPolicy.clearBudgetLimit();
    this->d->lodBudgetPolicy.clearRetainedQualityFloorBudget();
    this->d->lodBudgetPolicy.clearDeadlineCapacityCeiling();
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodPointProxyTriangleRecoveryPending = FALSE;
    this->d->lodPointProxyCalibrationPolicy.reset();
    this->d->lodCoveragePolicy.activate(false);
    this->d->lodViewDemandPolicy.refreshForPolicyRevision(
	preserveScaleDemandRefresh != FALSE,
	this->d->lodInteractive != FALSE);
    this->d->reconcilePhase(BObolLodStateMachine::Event::POLICY_CHANGED);
}

void
BObolViewController::beginLodInteraction(void)
{
    if (!this->d->lodAutoSubmit || this->d->lodGestureActive)
	return;

    this->d->lodHeadroomPolicy.cancelRetry();
    this->d->lodStaticOverscanActive = FALSE;
    const float previousPixelError = this->d->lodTargetPixelError;
    const bool newInteractionEpoch = !this->d->lodInteractive;
    int initialQualityFloor = -1;
    if (newInteractionEpoch) {
	BObolLodConvergenceStatus priorStatus;
	this->getLodConvergenceStatus(priorStatus);
	this->d->lodInteractionStartScaleSignature =
	    this->d->lodViewScaleSignature;
	this->d->lodInteractionStartScaleSignatureValid =
	    this->d->lodViewSignatureValid;
	this->d->lodInteractionStartedFromReadyView =
	    priorStatus.viewReady;
	this->d->lodStableBudgetBeforeInteraction =
	    this->d->lodBudgetPolicy.currentBudget();
	this->d->lodStableBudgetBeforeInteractionValid = TRUE;
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
	    this->d->pointProxyAggregationApplicable() ?
		this->d->lodPresentationPointProxyPixelThreshold : 1.0f,
	    controller_lod_presentation_population(snapshotState,
		this->d->lodViewQualityDomainRevision),
	    this->d->lodViewRevision);
	this->d->seedInteractiveCalibrationFromStable();
    }
    this->d->lodRetainPoseOccurrenceCuts = FALSE;
    this->d->lodRetainedImportanceCensusPending = FALSE;
    this->d->lodGestureActive = TRUE;
    this->d->lodViewDemandPolicy.beginGesture(newInteractionEpoch);
    if (newInteractionEpoch)
	this->d->lodViewDemandPolicy.seedQualityFloor(initialQualityFloor);
    this->d->lodDiscretePopulationTrialAvailable = FALSE;
    this->d->lodReleaseCutFloorActive = FALSE;
    this->d->lodPresentationPolicy.cancelHandoff();
    this->d->lodBudgetPolicy.resetProbeSeries();
    this->d->lodBudgetPolicy.resetOverloadRecovery();
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodPointProxyTriangleRecoveryPending = FALSE;
    this->d->lodInteractive = TRUE;
    this->d->lodSettleAfterRenderSerial = 0;
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodBudgetPolicy.resetPass();
    this->d->lodLastViewChangeMicroseconds = bu_gettime();
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
    if (this->d->pointProxyAggregationApplicable()) {
	this->d->lodPresentationPointProxyPixelThreshold =
	    BObolLodQualityPolicy::pointProxyThreshold(
		this->d->lodPresentationPointProxyPixelThreshold,
		interactionTimingSample,
		this->d->lodInteractiveTargetFps);
    } else {
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->lodPointProxyCalibrationPolicy.reset();
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
    this->d->reconcilePhase(
	BObolLodStateMachine::Event::INTERACTION_STARTED);
}

void
BObolViewController::endLodInteraction(void)
{
    if (!this->d->lodGestureActive)
	return;

    this->d->lodGestureActive = FALSE;
    this->d->lodViewDemandPolicy.endGesture();
    this->d->lodDiscretePopulationTrialAvailable = FALSE;
    {
	std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
	this->d->lastInteractivePresentationTimestampNanoseconds = 0;
    }
    BObolViewLodState *viewState =
	this->d->viewAttachment->getViewLodState();
    if (viewState)
	viewState->setCadPresentationCameraMotionFrameReuse(FALSE);
    /* Pose-only orthographic interaction never changed the occurrence cuts;
     * only renderer-side motion limits hid part of the resident presentation.
     * Restore the last measured stable limits on release, while retaining the
     * normal debounce for new projection/refinement work.  A zoom/projection
     * change must instead keep the responsive cut until it is re-budgeted. */
    BObolLodPresentationPolicy::RestoreInputs restoreInputs;
    restoreInputs.orthographic = this->d->activeCamera &&
	this->d->activeCamera->isOfType(
	    SoOrthographicCamera::getClassTypeId());
    restoreInputs.scaleChanged =
	this->d->lodViewDemandPolicy.interactionScaleChanged();
    restoreInputs.viewEpoch = this->d->lodViewRevision;
    restoreInputs.population =
	controller_lod_presentation_population(viewState,
	    this->d->lodViewQualityDomainRevision);
    restoreInputs.currentTargetPixelError = this->d->lodTargetPixelError;
    restoreInputs.currentProgressiveCeiling =
	this->d->lodInteractiveProgressiveCeiling;
    restoreInputs.currentPointProxyPixelThreshold =
	this->d->lodPresentationPointProxyPixelThreshold;
    const BObolLodPresentationPolicy::RestoreDecision restore =
	this->d->lodPresentationPolicy.restorePrior(restoreInputs);
    if (restore.apply) {
	/* Restore the semantic demand together with the renderer-only limits.
	 * The occurrence cuts remain untouched until the quiet exact-importance
	 * census, but leaving the motion pixel error active made diagnostics and
	 * any intervening producer result describe a coarse target after the rich
	 * presentation had already been reinstated. */
	this->d->lodTargetPixelError = restore.targetPixelError;
	this->d->lodInteractiveProgressiveCeiling =
	    restore.progressiveCeiling;
	this->d->lodPresentationPointProxyPixelThreshold =
	    restore.pointProxyPixelThreshold;
	this->d->lodReleaseCutFloorActive = FALSE;
	if (viewState) {
	    viewState->setCadPresentationProgressiveCutCeiling(
		this->d->lodInteractiveProgressiveCeiling);
	    viewState->setCadPresentationPointProxyPixelThreshold(
		this->d->lodPresentationPointProxyPixelThreshold);
	}
    }
    /*
     * Do not promote the face budget to activeFaceCount() here.  The active
     * retained population may be orders of magnitude richer than the
     * responsive cut actually shown through the renderer-wide PoP ceiling.
     * Zoom and changed-population releases retain that ceiling through the
     * quiet handoff; a validated pose-only snapshot was restored above.
     */
    if (!restore.restoredPriorStable)
	this->d->lodReleaseCutFloorActive = TRUE;
    /* A gesture pass may have exhausted the old, smaller allowance.  Force
     * the release pass to derive its admission state from the coherent floor
     * above rather than reusing that stale remainder. */
    this->d->lodBudgetPolicy.resetPass();
    /* Keep the interaction epoch through the normal quiet-view debounce.
     * Pose-only presentation may already be restored, but new projection and
     * refinement work still waits so release cannot block on a full software
     * planning frame. */
    this->d->lodLastViewChangeMicroseconds = bu_gettime();
    this->d->lodSettleAfterRenderSerial = 0;
    this->markProgressiveWorkPending();
    this->requestRender("lod-interaction-end");
    this->d->reconcilePhase(BObolLodStateMachine::Event::INTERACTION_ENDED);
}

SbBool
BObolViewController::isLodInteractionActive(void) const
{
    return this->d->lodInteractive;
}

SbBool
BObolViewController::isLodGestureActive(void) const
{
    return this->d->lodGestureActive;
}

SbBool
BObolViewController::isLodScaleChangingInteraction(void) const
{
    return this->d->lodViewDemandPolicy.scaleChangingInteraction(
	this->d->lodInteractive != FALSE) ? TRUE : FALSE;
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
    if (fabsf(this->d->lodInteractiveTargetFps - interactiveFps) <=
	    std::numeric_limits<float>::epsilon() &&
	fabsf(this->d->lodStableTargetFps - stableFps) <=
	    std::numeric_limits<float>::epsilon())
	return;
    this->d->resetLodViewQualityHistory();
    this->d->lodInteractiveTargetFps = interactiveFps;
    this->d->lodStableTargetFps = stableFps;
    this->d->lodStaticOverscanActive = FALSE;
    this->d->lodStaticOverscanRejected = FALSE;
    this->d->lodBudgetPolicy.resetPass();
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
    return this->d->lodBudgetPolicy.currentBudget();
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
    static_assert(
	static_cast<int>(BObolLodStateMachine::Phase::FALLBACK) ==
	    BOBOL_LOD_COORDINATOR_FALLBACK &&
	static_cast<int>(BObolLodStateMachine::Phase::COVERAGE) ==
	    BOBOL_LOD_COORDINATOR_COVERAGE &&
	static_cast<int>(BObolLodStateMachine::Phase::INTERACTIVE) ==
	    BOBOL_LOD_COORDINATOR_INTERACTIVE &&
	static_cast<int>(BObolLodStateMachine::Phase::SETTLING) ==
	    BOBOL_LOD_COORDINATOR_SETTLING &&
	static_cast<int>(BObolLodStateMachine::Phase::STABLE) ==
	    BOBOL_LOD_COORDINATOR_STABLE &&
	static_cast<int>(BObolLodStateMachine::Phase::COMPACTING) ==
	    BOBOL_LOD_COORDINATOR_COMPACTING,
	"public and private LoD coordinator phases must agree");
    static_assert(
	static_cast<int>(BObolLodStateMachine::Event::INITIALIZE) ==
	    BOBOL_LOD_EVENT_INITIALIZE &&
	static_cast<int>(BObolLodStateMachine::Event::FRAME_COMPLETED) ==
	    BOBOL_LOD_EVENT_FRAME_COMPLETED &&
	static_cast<int>(BObolLodStateMachine::Event::WORK_SCHEDULED) ==
	    BOBOL_LOD_EVENT_WORK_SCHEDULED &&
	static_cast<int>(BObolLodStateMachine::Event::WORK_PUMPED) ==
	    BOBOL_LOD_EVENT_WORK_PUMPED &&
	static_cast<int>(BObolLodStateMachine::Event::RESULT_PUBLISHED) ==
	    BOBOL_LOD_EVENT_RESULT_PUBLISHED &&
	static_cast<int>(BObolLodStateMachine::Event::SERVICE_CHANGED) ==
	    BOBOL_LOD_EVENT_SERVICE_CHANGED &&
	static_cast<int>(BObolLodStateMachine::Event::GENERATION_CANCELLED) ==
	    BOBOL_LOD_EVENT_GENERATION_CANCELLED &&
	static_cast<int>(BObolLodStateMachine::Event::AUTO_SUBMIT_CHANGED) ==
	    BOBOL_LOD_EVENT_AUTO_SUBMIT_CHANGED &&
	static_cast<int>(BObolLodStateMachine::Event::VIEW_INVALIDATED) ==
	    BOBOL_LOD_EVENT_VIEW_INVALIDATED &&
	static_cast<int>(BObolLodStateMachine::Event::POLICY_CHANGED) ==
	    BOBOL_LOD_EVENT_POLICY_CHANGED &&
	static_cast<int>(BObolLodStateMachine::Event::INTERACTION_STARTED) ==
	    BOBOL_LOD_EVENT_INTERACTION_STARTED &&
	static_cast<int>(BObolLodStateMachine::Event::INTERACTION_ENDED) ==
	    BOBOL_LOD_EVENT_INTERACTION_ENDED &&
	static_cast<int>(BObolLodStateMachine::Event::VIEW_OBSERVED) ==
	    BOBOL_LOD_EVENT_VIEW_OBSERVED,
	"public and private LoD coordinator events must agree");
    const BObolLodStateMachine::Phase coordinatorPhase =
	this->d->currentPhase();
    status.coordinatorPhase = static_cast<int>(coordinatorPhase);
    status.coordinatorEvent =
	static_cast<int>(this->d->lastDispatchedEvent());
    status.coordinatorTransitionSerial =
	this->d->phaseTransitionSerial();
    status.coordinatorProgressSequence =
	this->d->phaseWitness(coordinatorPhase).sequence;
    status.coordinatorDispatchSerial = this->d->dispatchSerial();
    status.coordinatorStagnantDispatchCount =
	this->d->stagnantDispatchCount();
    status.coordinatorInvariantViolationCount =
	this->d->invariantViolationCount();
    status.coordinatorInvariantMask = this->d->lastInvariantMask();
    status.coordinatorInvariantHistoryMask =
	this->d->invariantHistoryMask();
    status.viewQualityHistoryEntryCount =
	this->d->lodViewQualityHistory.size();
    status.viewQualityHistoryRememberCount =
	this->d->lodViewQualityHistory.rememberCount();
    status.viewQualityHistoryRecallCount =
	this->d->lodViewQualityHistory.recallCount();
    status.viewRevision = this->d->lodViewRevision.value();
    status.failedSourceCount = this->getLastFailedSourceCount();
    status.visibleTargetCount =
	this->d->lodConvergenceCandidateCount();
    status.activeFaces = this->getActiveLodFaceCount();
    status.activeRenderCost = this->getActiveLodRenderCost();
    status.renderCostBudget = this->getCurrentLodRenderCostBudget();

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
    }

    const SbBool structuralPending =
	status.expectedLeafCount > status.availableLeafCount;
    const SbBool sourcePreparationPending =
	this->d->progressiveProviderPendingCount > 0 ||
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
	this->d->lodPublicationPolicy.pending();
    const SbBool calibrationPending =
	this->hasPendingLodRefinementFrame();
    status.refinementFramePending =
	this->d->lodFrameObligation.pending() ? TRUE : FALSE;
    status.activeGeneration = this->d->lodActiveGeneration;
    status.submissionSourceIndex = this->d->lodSubmissionSourceIndex;
    status.submissionEntryOffset = this->d->lodSubmissionEntryOffset;
    status.budgetCalibrationPending =
	this->d->lodBudgetPolicy.rescanAfterFrame() ||
	this->d->lodHeadroomPolicy.retryPending();
    status.stablePresentationHandoffPending =
	this->d->lodPresentationPolicy.handoffPending() ? TRUE : FALSE;
    status.pointProxyCalibrationPending =
	(this->d->lodAdmissionPointProxyFramePending ||
	 this->d->lodStablePointProxyCalibrationPending ||
	 this->d->lodPointProxyTriangleRecoveryPending) ? TRUE : FALSE;
    status.pointProxyAdmissionFramePending =
	this->d->lodAdmissionPointProxyFramePending;
    status.stablePointProxyCalibrationPending =
	this->d->lodStablePointProxyCalibrationPending;
    status.pointProxyTriangleRecoveryPending =
	this->d->lodPointProxyTriangleRecoveryPending;
    status.residentGrowthReallocationPending =
	this->d->lodResidentGrowthPolicy.pending() ? TRUE : FALSE;
    status.publicationFramePending = publicationPending;
    status.sourcePreparationPending = sourcePreparationPending;
    status.sourcePreparationProviderCount =
	this->d->progressiveProviderPendingCount;

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
	this->d->lodSubmissionPending != FALSE;
    convergenceInputs.resultPending = resultPending != FALSE;
    convergenceInputs.publicationPending = publicationPending != FALSE;
    convergenceInputs.calibrationPending = calibrationPending != FALSE;
    convergenceInputs.interactive = this->d->lodInteractive != FALSE;
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
	this->d->lodBudgetPolicy.stableBudgetLimited();
    convergenceInputs.presentationLimited =
	!this->d->forceTerminalLodRefinement &&
	(this->d->lodInteractiveProgressiveCeiling >= 0 ||
	 this->d->lodPresentationPointProxyPixelThreshold > 1.01f ||
	 /* A rejected static-quality trial is completed-frame proof that the
	  * next richer population missed this view epoch's hard presentation
	  * deadline.  Its temporary renderer ceiling may already have been
	  * reconciled and removed, but that must not make the accepted terminal
	  * cut appear unconstrained to the HUD or qualification harness. */
	 this->d->lodStaticOverscanRejected);
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
    status.phase = static_cast<int>(convergence.phase);
    status.fraction = convergence.fraction;
    status.viewReady = convergence.viewReady ? TRUE : FALSE;
    status.backgroundPending =
	convergence.backgroundPending ? TRUE : FALSE;
    status.memoryLimited = convergence.memoryLimited ? TRUE : FALSE;
    status.performanceLimited =
	convergence.performanceLimited ? TRUE : FALSE;
}

double
BObolViewController::getCalibratedLodRenderCostPerSecond(void) const
{
    const long double value = this->d->lodInteractive ?
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
	this->d->lodViewSignatureValid;
    const BObolLodViewScaleSnapshot priorScaleSignature =
	this->d->lodViewScaleSignature;
    /* advanceLodViewRevision() clears the view-local visibility census.  Its
     * transient zero must not discard an already proven aggregate point cut;
     * sample this presentation fact before invalidating the census. */
    const SbBool priorPointProxyAggregationApplicable =
	BObolLodPointProxyCalibrationPolicy::
	    applicableAcrossCameraInvalidation(
		this->d->pointProxyAggregationApplicable(),
		this->d->lodPresentationPointProxyPixelThreshold) ?
	    TRUE : FALSE;
    SbBool priorViewReady = FALSE;
    if (this->d->lodAutoSubmit && !this->d->lodInteractive) {
	BObolLodConvergenceStatus priorStatus;
	this->getLodConvergenceStatus(priorStatus);
	priorViewReady = priorStatus.viewReady;
    }

    if (this->d->lodViewSignatureValid &&
	this->d->lodViewSignature.same(signatures.view))
	return;

    const SbBool scaleChanged =
	this->d->lodViewSignatureValid &&
	!this->d->lodViewScaleSignature.same(signatures.scale) ?
	    TRUE : FALSE;
    this->d->lodViewSignature = signatures.view;
    this->d->lodViewScaleSignature = signatures.scale;
    this->d->lodViewSignatureValid = TRUE;
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
	this->d->lodDiscretePopulationTrialAvailable = FALSE;
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
	     * setting lodInteractive below.
	     */
	    const SbBool continuingInteractive =
		this->d->lodInteractive;
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
		this->d->lodInteractionStartScaleSignature =
		    priorScaleSignature;
		this->d->lodInteractionStartScaleSignatureValid =
		    priorViewSignatureValid;
		this->d->lodInteractionStartedFromReadyView =
		    priorViewReady;
		this->d->lodStableBudgetBeforeInteraction =
		    this->d->lodBudgetPolicy.currentBudget();
		this->d->lodStableBudgetBeforeInteractionValid = TRUE;
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
		    pointProxyAggregationForCameraChange ?
			this->d->lodPresentationPointProxyPixelThreshold : 1.0f,
		    controller_lod_presentation_population(snapshotState,
			this->d->lodViewQualityDomainRevision),
		    this->d->lodViewRevision);
		this->d->seedInteractiveCalibrationFromStable();
	    }
	    this->d->lodRetainPoseOccurrenceCuts = FALSE;
	    this->d->lodRetainedImportanceCensusPending = FALSE;
	    this->d->lodViewDemandPolicy.beginCameraInteraction(
		!continuingInteractive, scaleChanged != FALSE);
	    if (!continuingInteractive)
		this->d->lodViewDemandPolicy.seedQualityFloor(
		    initialQualityFloor);
	    const int64_t now = bu_gettime();
	    this->d->lodLastViewChangeMicroseconds = now;
	    this->d->lodInteractive = TRUE;
	    this->d->lodReleaseCutFloorActive = FALSE;
	    this->d->lodPresentationPolicy.cancelHandoff();
	    this->d->lodBudgetPolicy.resetProbeSeries();
	    this->d->lodBudgetPolicy.resetOverloadRecovery();
	    this->d->lodRefinementNotBeforeMicroseconds = 0;
	    this->d->lodBudgetPolicy.resetPass();
	    this->d->lodSettleAfterRenderSerial =
		this->d->renderCompletionSerial + 1;
	    if (this->d->lodSettleAfterRenderSerial == 0)
		this->d->lodSettleAfterRenderSerial = 1;
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
		this->d->lodPointProxyCalibrationPolicy.reset();
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
    this->d->reconcilePhase(BObolLodStateMachine::Event::VIEW_OBSERVED);
}
