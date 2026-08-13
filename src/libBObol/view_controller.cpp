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
#include "raytrace.h"
#include "rt/view.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
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
static std::vector<SoBRLDatabaseSource *> controller_render_database_sources(
    const BObolViewController *controller);
static std::vector<SoBRLDatabaseSource *> controller_render_database_source_roots(
    const BObolViewController *controller);
static std::vector<SoBRLMeshShape *> controller_render_mesh_shapes(
    const BObolViewController *controller);

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

/* Enter recognizable-geometry protection before the public diagnostic's
 * sixteen-pixel "prominent" boundary.  The hysteresis prevents a discrete
 * PoP or point-proxy transition from moving an occurrence directly into a
 * reported deficit. */
static constexpr double controller_lod_protected_footprint_pixels = 12.0;

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

/* SoCADAssembly replaces an unselected micro-occurrence with one point in its
 * retained aggregate point channel.  Scene admission must price the same
 * presentation.  Charging the occurrence's minimum triangle prefix after the
 * renderer has collapsed it creates a fictitious per-leaf floor: on a 150k
 * scene the allocator can report that the minimum population exceeds the
 * measured budget while the completed frame still has ample room for richer
 * wheels, blades, panels, and other visible meshes.
 *
 * The renderer applies hysteresis around its threshold.  Use its conservative
 * new-collapse boundary here.  This may under-count a little reusable point
 * capacity, but can never reserve triangle headroom for an occurrence which
 * the renderer has not yet proved collapsible. */
static bool
controller_lod_point_proxy_candidate(
    SoBRLDatabaseSource *source,
    const BObolViewLodState::CadPayload *payload,
    float pointProxyPixelThreshold,
    BObolCompactLodPlanningSummary *planningSummary = NULL)
{
    if (!source || !payload || !payload->progressiveMesh ||
	!payload->progressiveMesh->isValid() ||
	!std::isfinite(pointProxyPixelThreshold) ||
	pointProxyPixelThreshold <= 0.0f ||
	!std::isfinite(payload->projectedPixelDiameter) ||
	payload->projectedPixelDiameter <= 0.0f)
	return false;

    if (!payload->projectedBoundsContained)
	return false;

    BObolCompactLodPlanningSummary summary;
    if (!source->getCompactLodPlanningSummaryForKey(
	    payload->sourceInstanceKey.getString(), summary) ||
	!summary.valid || !summary.visible || summary.selected ||
	summary.highlighted)
	return false;
    if (planningSummary)
	*planningSummary = summary;

    /* Adaptive point aggregation bounds draw-call cost; it is never a
     * license to replace a recognizable wheel, blade, panel, or silhouette
     * part with one pixel.  The renderer receives the same per-instance
     * protection below, so allocator accounting and presentation agree even
     * when the global pressure threshold reaches tens of pixels. */
    if (controller_lod_visual_footprint(payload) >=
	    controller_lod_protected_footprint_pixels)
	return false;

    return payload->projectedPixelDiameter <=
	pointProxyPixelThreshold * 0.75f;
}

static size_t
controller_lod_point_proxy_render_cost(int drawMode)
{
    BObolLodCounts point;
    point.pointCount = 1;
    /* Point proxies share one aggregate draw and therefore do not carry the
     * retained-mesh per-occurrence command cost. */
    return bobol_lod_render_cost_units(point, drawMode, 0);
}

static void
controller_lod_effective_population_cost(
    const std::vector<SoBRLDatabaseSource *> &sources,
    const BObolViewLodState *viewState, float pointProxyPixelThreshold,
    const SbMatrix *viewProjection, size_t &activeCost, size_t &minimumCost)
{
    activeCost = 0;
    minimumCost = 0;
    if (!viewState)
	return;

    const auto add = [](size_t total, size_t value) {
	return value > SIZE_MAX - total ? SIZE_MAX : total + value;
    };
    for (SoBRLDatabaseSource *source : sources) {
	if (!source)
	    continue;
	std::vector<const BObolViewLodState::CadPayload *> payloads;
	viewState->findCadPayloadsUnordered(source, payloads);
	for (const BObolViewLodState::CadPayload *payload : payloads) {
	    if (!payload || !payload->isValid())
		continue;
	    const int drawMode = source->getEffectiveLodDrawMode();
	    if (controller_lod_point_proxy_candidate(source, payload,
		    pointProxyPixelThreshold)) {
		const size_t cost = controller_lod_point_proxy_render_cost(
		    drawMode);
		activeCost = add(activeCost, cost);
		minimumCost = add(minimumCost, cost);
		continue;
	    }

	    std::array<BObolLodCounts,
		BOBOL_MESH_LOD_CUT_COUNT_MAX> visibleCounts;
	    bool haveVisibleCounts = false;
	    if (viewProjection && payload->progressiveMesh &&
		payload->progressiveMesh->hasSpatialClusters()) {
		BObolCompactLodPlanningSummary summary;
		if (source->getCompactLodPlanningSummaryForKey(
			payload->sourceInstanceKey.getString(), summary) &&
		    summary.valid && summary.visible)
		    haveVisibleCounts = payload->progressiveMesh->
			visibleCountsAtCuts(summary.localToSource,
			    *viewProjection, payload->hasNormals,
			    visibleCounts.data(), visibleCounts.size());
	    }
	    const BObolLodCounts &activeCounts = haveVisibleCounts &&
		payload->activeCut >= 0 ?
		visibleCounts[static_cast<size_t>(payload->activeCut)] :
		payload->counts;
	    activeCost = add(activeCost, bobol_lod_render_cost_units(
		activeCounts, drawMode, 1));
	    BObolLodCounts minimumCounts = payload->counts;
	    if (payload->progressiveMesh &&
		payload->progressiveMesh->isValid()) {
		const int minimumCut =
		    payload->progressiveMesh->minimumCut();
		if (minimumCut >= 0 &&
		    payload->progressiveMesh->canDrawCut(minimumCut))
		    minimumCounts = haveVisibleCounts ?
			visibleCounts[static_cast<size_t>(minimumCut)] :
			bobol_lod_progressive_counts(
			    payload->progressiveMesh, minimumCut,
			    payload->hasNormals);
	    }
	    minimumCost = add(minimumCost, bobol_lod_render_cost_units(
		minimumCounts, drawMode, 1));
	}
    }
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
	std::vector<size_t> &entryIndices)
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
	    viewState->findCadForOccurrence(source, key);
	const int emphasis = haveSummary && summary.selected ? 2 :
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
	std::vector<size_t> &entryIndices)
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
	if (!source->getCompactInstanceIndex(
		payload->sourceInstanceKey.getString(), entryIndex))
	    continue;
	BObolCompactLodPlanningSummary summary;
	const bool haveSummary = source->getCompactLodPlanningSummary(
	    static_cast<int>(entryIndex), summary) && summary.valid &&
	    summary.visible;
	if (!haveSummary)
	    continue;
	const int emphasis = summary.selected ? 2 :
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

/* Compute one scene-wide minimax screen-error target for a retained
 * population.  All PoP prefixes are already immutable and resident; this is
 * a discrete draw-count allocation only.  A logarithmic binary search is
 * O(N) per iteration, independent of hierarchy depth, and prevents an early
 * rich object from consuming the allowance before another prominent object
 * receives a recognizable cut. */
static double
controller_lod_retained_error_ceiling(
    const std::vector<SoBRLDatabaseSource *> &sources,
    BObolViewLodState *viewState, size_t sceneBudget,
    bool allowProtectedFloor, size_t maximumProtectedBudget,
    size_t *allocationBudget, uint64_t *allocationSignature,
    uint64_t viewRevision, uint64_t policyRevision,
    float pointProxyPixelThreshold, const SbMatrix *viewProjection)
{
    struct Candidate {
	const BObolViewLodState::CadPayload *payload;
	BObolLodProgressiveMeshPtr mesh;
	int minimumCut;
	int maximumCut;
	int drawMode;
	SbBool hasNormals;
	double projectedPixelDiameter;
	double errorWeight;
	int protectedCut;
	std::unique_ptr<std::array<BObolLodCounts,
	    BOBOL_MESH_LOD_CUT_COUNT_MAX>> visibleCounts;
    };
    if (allocationBudget)
	*allocationBudget = 0;
    if (allocationSignature)
	*allocationSignature = 0;
    if (!viewState || sceneBudget == SIZE_MAX)
	return std::numeric_limits<double>::infinity();

    std::vector<Candidate> candidates;
    candidates.reserve(viewState->cadMeshPayloadCount());
    std::unordered_map<SoCADAssembly *, std::vector<Obol::InstanceId>>
	pointProxyProtectedInstances;
    size_t fixedCost = 0;
    double maximumMinimumError = 1.0;
    for (SoBRLDatabaseSource *source : sources) {
	if (!source)
	    continue;
	SoCADAssembly *cadAssembly = viewState->findCadPresentation(source);
	if (cadAssembly)
	    (void)pointProxyProtectedInstances[cadAssembly];
	std::vector<const BObolViewLodState::CadPayload *> payloads;
	viewState->findCadPayloadsUnordered(source, payloads);
	for (const BObolViewLodState::CadPayload *payload : payloads) {
	    if (!payload || !payload->isValid())
		continue;
	    if (!payload->progressiveMesh ||
		!payload->progressiveMesh->isValid()) {
		const size_t cost = bobol_lod_render_cost_units(
		    payload->counts, source->getEffectiveLodDrawMode(), 1);
		fixedCost = cost > SIZE_MAX - fixedCost ? SIZE_MAX :
		    fixedCost + cost;
		continue;
	    }
	    const int minimumCut = payload->progressiveMesh->minimumCut();
	    /* Budget the complete view-requested hierarchy, not only the prefix
	     * which happens to be resident at this instant.  Immutable hierarchy
	     * counts are available before its suffix is loaded.  Capping this at
	     * residentCut let worker completions make later bounded windows
	     * eligible for unreserved detail, exhausting the allowance before
	     * other prominent occurrences were reached.  Presentation still
	     * activates only canDrawCut() cuts; this is capacity reservation. */
	    int maximumCut = payload->requestedCut;
	    maximumCut = std::max(minimumCut,
		std::min(maximumCut,
		    payload->progressiveMesh->maximumCut()));
	    BObolCompactLodPlanningSummary summary;
	    const bool haveSummary = source->getCompactLodPlanningSummaryForKey(
		payload->sourceInstanceKey.getString(), summary) &&
		summary.valid && summary.visible;
	    if (haveSummary && controller_lod_point_proxy_candidate(
		    source, payload, pointProxyPixelThreshold)) {
		/* Preserve one drawable mesh prefix behind the reversible point
		 * presentation, but charge only the aggregate point.  The allocated
		 * minimum prevents a stale richer allocation from consuming GPU atlas
		 * space; the immutable resident suffix remains available when the
		 * threshold relaxes or the view zooms in. */
		const size_t pointCost =
		    controller_lod_point_proxy_render_cost(
			source->getEffectiveLodDrawMode());
		fixedCost = pointCost > SIZE_MAX - fixedCost ? SIZE_MAX :
		    fixedCost + pointCost;
		(void)viewState->setCadAllocatedCut(payload, minimumCut,
		    viewRevision, policyRevision,
		    source->getEffectiveLodDrawMode());
		continue;
	    }
	    const double emphasis = haveSummary && summary.selected ? 4.0 :
		(haveSummary && summary.highlighted ? 2.0 : 1.0);
	    const double target = std::max(1.0,
		static_cast<double>(payload->targetPixelError));
	    /* Equal normalized error is not equal visual harm.  A large wheel,
	     * blade, panel, or silhouette segment remains conspicuous at an error
	     * which is harmless on a four-pixel fastener.  Weight the minimax
	     * demand continuously by projected footprint: two pixels and below
	     * retain weight one, and the diagnostic "prominent" boundary at
	     * sixteen pixels reaches the bounded maximum weight of eight.  That
	     * stronger fallback matters only when the explicit protected floor
	     * below is mathematically infeasible inside the measured allowance.
	     *
	     * This is intentionally a bounded screen-space multiplier, not object
	     * identity or model-space size.  It survives arbitrary vehicle scale,
	     * changes correctly with the view, and cannot let one huge hull consume
	     * an unbounded share before other visible parts receive a useful cut. */
	    const double area = std::max(0.0,
		static_cast<double>(payload->projectedPixelArea));
	    const double perimeter = std::max(0.0,
		static_cast<double>(payload->projectedPixelPerimeter));
	    const double diameter = std::max(0.0,
		static_cast<double>(payload->projectedPixelDiameter));
	    const double footprint = std::max(std::sqrt(area),
		std::max(perimeter * 0.25, diameter * 0.25));
	    const double significance = std::max(1.0,
		std::min(8.0, footprint * 0.50));
	    double errorWeight = emphasis * significance / target;
	    if (!std::isfinite(errorWeight) || errorWeight <= 0.0)
		errorWeight = 1.0;
	    Candidate candidate;
	    candidate.payload = payload;
	    candidate.mesh = payload->progressiveMesh;
	    candidate.minimumCut = minimumCut;
	    candidate.maximumCut = maximumCut;
	    /* A retained progressive generation is shared by shaded and wire
	     * presentations.  The payload records the mode which last published
	     * it, but the admission action charges the source's current effective
	     * mode.  Price the same current channel here or the minimax population
	     * can be mathematically feasible in one currency and over budget when
	     * its bounded windows apply it in another. */
	    candidate.drawMode = source->getEffectiveLodDrawMode();
	    candidate.hasNormals = payload->hasNormals;
	    candidate.projectedPixelDiameter = diameter;
	    candidate.errorWeight = errorWeight;
	    if (haveSummary && viewProjection &&
		candidate.mesh->hasSpatialClusters()) {
		candidate.visibleCounts.reset(new std::array<BObolLodCounts,
		    BOBOL_MESH_LOD_CUT_COUNT_MAX>);
		if (!candidate.mesh->visibleCountsAtCuts(
		    summary.localToSource, *viewProjection,
		    candidate.hasNormals, candidate.visibleCounts->data(),
		    candidate.visibleCounts->size()))
		    candidate.visibleCounts.reset();
	    }
	    /* Reserve recognizable quality before the soft minimax objective.
	     * Selection receives pixel-exact priority, highlight a two-pixel floor,
	     * and ordinary prominent parts the three-pixel visual contract.  The
	     * previous two-pixel reservation commonly selected a whole additional
	     * PoP cut for Hubble's hundreds of visible parts, pushing an otherwise
	     * useful software frame from roughly 50 ms to the 100 ms hard boundary.
	     * Enter protection at
	     * twelve pixels, below the diagnostic's sixteen-pixel boundary.  This
	     * deliberate hysteresis absorbs the small metric changes introduced by
	     * discrete PoP retargeting and bounded current-view scans; an object may
	     * not cross into the reported prominent set immediately after an
	     * allocation with no reserved cut.  Non-prominent parts retain their
	     * minimum mesh until the second-stage allocator finds them affordable. */
	    const double protectedError = haveSummary && summary.selected ? 1.0 :
		(haveSummary && summary.highlighted ? 2.0 :
		    (footprint >= controller_lod_protected_footprint_pixels ?
			3.0 :
			std::numeric_limits<double>::infinity()));
	    candidate.protectedCut = minimumCut;
	    if (std::isfinite(protectedError) && diameter > 0.0) {
		const int protectedCut = candidate.mesh->cutForScreenError(
		    diameter, protectedError);
		if (protectedCut >= 0)
		    candidate.protectedCut = std::max(minimumCut,
			std::min(maximumCut, protectedCut));
	    }
	    if (std::isfinite(protectedError) && cadAssembly) {
		size_t entryIndex = 0;
		BObolCompactInstanceHandle handle;
		if (source->getCompactInstanceIndex(
			payload->sourceInstanceKey.getString(), entryIndex) &&
		    entryIndex <= static_cast<size_t>(
			std::numeric_limits<int>::max()) &&
		    source->getCompactInstanceHandle(
			static_cast<int>(entryIndex), handle) &&
		    handle.isValid()) {
		    Obol::InstanceId instance;
		    instance.w0 = handle.instanceWord0;
		    instance.w1 = handle.instanceWord1;
		    pointProxyProtectedInstances[cadAssembly].push_back(instance);
		}
	    }
	    candidates.push_back(std::move(candidate));
	    const double minimumError =
		payload->progressiveMesh->projectedErrorAtCut(
		    minimumCut, diameter) * errorWeight;
	    if (std::isfinite(minimumError))
		maximumMinimumError = std::max(
		    maximumMinimumError, minimumError);
	}
    }

    /* SoCADAssembly owns the actual aggregate-point classification.  Publish
     * the allocator's complete view-local protection set atomically per
     * assembly so a high scene-pressure threshold cannot hide protected
     * instances behind the allocator's retained mesh accounting. */
    for (auto &protectedSet : pointProxyProtectedInstances)
	protectedSet.first->setPointProxyProtectedInstances(protectedSet.second);
    if (candidates.empty())
	return std::numeric_limits<double>::infinity();

    size_t protectedFloorCost = fixedCost;
    uint64_t protectedFloorSignature =
	static_cast<uint64_t>(candidates.size()) * 0x9e3779b97f4a7c15ULL;
    const auto countsAtCut = [](const Candidate &candidate, int cut) {
	return candidate.visibleCounts ?
	    (*candidate.visibleCounts)[static_cast<size_t>(cut)] :
	    bobol_lod_progressive_counts(
		candidate.mesh, cut, candidate.hasNormals);
    };
    for (const Candidate &candidate : candidates) {
	const BObolLodCounts counts = countsAtCut(
	    candidate, candidate.protectedCut);
	const size_t candidateCost = bobol_lod_render_cost_units(
	    counts, candidate.drawMode, 1);
	protectedFloorCost = candidateCost > SIZE_MAX - protectedFloorCost ?
	    SIZE_MAX : protectedFloorCost + candidateCost;
	/* Runtime identity is sufficient here: this witness never leaves the
	 * controller or survives a process.  Combine an order-independent mixed
	 * occurrence/cut token because payload discovery deliberately has no
	 * ordering contract. */
	uint64_t token = static_cast<uint64_t>(
	    reinterpret_cast<uintptr_t>(candidate.payload));
	token ^= static_cast<uint64_t>(
	    static_cast<uint32_t>(candidate.protectedCut)) << 32;
	token ^= static_cast<uint64_t>(
	    static_cast<uint32_t>(candidate.drawMode)) << 16;
	token ^= candidate.hasNormals ? 0xd6e8feb86659fd93ULL : 0;
	token += 0x9e3779b97f4a7c15ULL;
	token = (token ^ (token >> 30)) * 0xbf58476d1ce4e5b9ULL;
	token = (token ^ (token >> 27)) * 0x94d049bb133111ebULL;
	token ^= token >> 31;
	protectedFloorSignature ^= token;
    }
    protectedFloorSignature ^=
	static_cast<uint64_t>(protectedFloorCost) * 0x517cc1b727220a95ULL;
    if (allocationSignature)
	*allocationSignature = protectedFloorSignature;
    /* The ordinary stable allowance targets 20 Hz, but falling just below
     * that target must not make a large wheel, blade, panel, or silhouette
     * unrecognizable.  Permit the protected floor to spend calibrated
     * headroom down to the existing 10 Hz hard quality deadline.  If even
     * that bounded allowance cannot hold the complete floor, fall back to
     * the soft footprint-weighted minimax solution. */
    const size_t protectedBudget = std::max(
	sceneBudget, maximumProtectedBudget);
    /* PoP cuts are discrete and the throughput estimator is deliberately
     * safety-scaled.  A complete recognizable population can therefore sit
     * just above its predicted allowance even when the current frame has
     * ample room below the hard quality deadline.  Permit one bounded direct
     * trial up to 25 percent beyond that estimate.  The installed floor is a
     * fixed point; repeated ceiling-free deadline misses explicitly reject it
     * for this view, so this trial cannot become an admission oscillation. */
    const size_t protectedTrialBudget =
	allowProtectedFloor && maximumProtectedBudget > 0 ?
	(protectedBudget > SIZE_MAX - protectedBudget / 4 ? SIZE_MAX :
	    protectedBudget + protectedBudget / 4) : protectedBudget;
    const bool enforceProtectedFloor =
	allowProtectedFloor && protectedFloorCost <= protectedTrialBudget;
    const size_t effectiveBudget = enforceProtectedFloor ?
	std::max(sceneBudget, protectedFloorCost) : sceneBudget;
    if (allocationBudget)
	*allocationBudget = enforceProtectedFloor ? protectedFloorCost : 0;
    const auto allocatedCut = [enforceProtectedFloor](
	const Candidate &candidate, double ceiling) {
	const double targetError = candidate.errorWeight > 0.0 ?
	    ceiling / candidate.errorWeight : ceiling;
	int minimaxCut = candidate.mesh->cutForScreenError(
	    candidate.projectedPixelDiameter, targetError);
	if (minimaxCut < 0)
	    minimaxCut = candidate.minimumCut;
	minimaxCut = std::max(candidate.minimumCut,
	    std::min(candidate.maximumCut, minimaxCut));
	return enforceProtectedFloor ?
	    std::max(minimaxCut, candidate.protectedCut) : minimaxCut;
    };
    const auto populationCost = [&candidates, fixedCost,
	&allocatedCut, &countsAtCut](double ceiling) {
	size_t cost = fixedCost;
	for (const Candidate &candidate : candidates) {
	    const int cut = allocatedCut(candidate, ceiling);
	    const BObolLodCounts counts = countsAtCut(candidate, cut);
	    const size_t candidateCost = bobol_lod_render_cost_units(
		counts, candidate.drawMode, 1);
	    cost = candidateCost > SIZE_MAX - cost ? SIZE_MAX :
		cost + candidateCost;
	}
	return cost;
    };

    if (getenv("BOBOL_LOD_TRACE_BUDGET"))
	bu_log("BObol LoD retained allocator candidates=%zu payloads=%zu "
	       "fixed_cost=%zu calculated_minimum=%zu view_minimum=%zu "
	       "protected_floor=%zu protected=%d budget=%zu "
	       "maximum_minimum_error=%.9g\n",
	       candidates.size(), viewState->cadMeshPayloadCount(), fixedCost,
	       populationCost(maximumMinimumError),
	       viewState->minimumActiveRenderCost(), protectedFloorCost,
	       enforceProtectedFloor ? 1 : 0, effectiveBudget,
	       maximumMinimumError);

    double ceiling = 1.0;
    if (populationCost(ceiling) > effectiveBudget) {
	ceiling = maximumMinimumError;
	if (populationCost(ceiling) <= effectiveBudget) {
	    double lowLog = 0.0;
	    double highLog = std::log2(maximumMinimumError);
	    for (int iteration = 0; iteration < 24; ++iteration) {
		const double middleLog = (lowLog + highLog) * 0.5;
		const double middle = std::exp2(middleLog);
		if (populationCost(middle) <= effectiveBudget)
		    highLog = middleLog;
		else
		    lowLog = middleLog;
	    }
	    ceiling = std::exp2(highLog) * (1.0 + 1.0e-9);
	}
    }

    /* A common scalar ceiling is fair, but discrete PoP cuts can make it
     * badly under-fill the allowance: one infinitesimal threshold change may
     * advance hundreds of symmetric parts at once.  Preserve the minimax
     * result as the baseline, then spend otherwise stranded headroom on
     * individual marginal upgrades.  Quality-floor violations sort first;
     * the remaining order is projected error reduction per render-cost unit.
     * No geometry is inspected or copied here, and each accepted transition
     * advances one producer-authored cut, so the work is bounded by the number
     * of transitions the calibrated scene budget can actually afford. */
    std::vector<int> finalCuts(candidates.size(), -1);
    std::vector<size_t> finalCosts(candidates.size(), 0);
    size_t finalCost = fixedCost;
    for (size_t i = 0; i < candidates.size(); ++i) {
	finalCuts[i] = allocatedCut(candidates[i], ceiling);
	const BObolLodCounts counts = countsAtCut(
	    candidates[i], finalCuts[i]);
	finalCosts[i] = bobol_lod_render_cost_units(
	    counts, candidates[i].drawMode, 1);
	finalCost = finalCosts[i] > SIZE_MAX - finalCost ? SIZE_MAX :
	    finalCost + finalCosts[i];
    }

    struct MarginalUpgrade {
	size_t candidateIndex;
	int nextCut;
	size_t nextCost;
	size_t addedCost;
	bool qualityFloorViolation;
	double weightedError;
	double valuePerCost;
    };
    struct MarginalUpgradeLess {
	bool operator()(const MarginalUpgrade &a,
		const MarginalUpgrade &b) const
	{
	    if (a.qualityFloorViolation != b.qualityFloorViolation)
		return !a.qualityFloorViolation;
	    if (a.qualityFloorViolation &&
		(a.weightedError < b.weightedError ||
		 a.weightedError > b.weightedError))
		return a.weightedError < b.weightedError;
	    if (a.valuePerCost < b.valuePerCost ||
		a.valuePerCost > b.valuePerCost)
		return a.valuePerCost < b.valuePerCost;
	    if (a.weightedError < b.weightedError ||
		a.weightedError > b.weightedError)
		return a.weightedError < b.weightedError;
	    return a.candidateIndex > b.candidateIndex;
	}
    };
    std::priority_queue<MarginalUpgrade,
	std::vector<MarginalUpgrade>, MarginalUpgradeLess> upgrades;
    const auto queueNextUpgrade = [&candidates, &finalCuts, &finalCosts,
	&upgrades, &countsAtCut](size_t candidateIndex) {
	const Candidate &candidate = candidates[candidateIndex];
	const int currentCut = finalCuts[candidateIndex];
	if (currentCut < candidate.minimumCut ||
	    currentCut >= candidate.maximumCut)
	    return;
	const double currentError = candidate.mesh->projectedErrorAtCut(
	    currentCut, candidate.projectedPixelDiameter);
	for (int nextCut = currentCut + 1;
	    nextCut <= candidate.maximumCut; ++nextCut) {
	    const BObolLodCounts counts = countsAtCut(candidate, nextCut);
	    const size_t nextCost = bobol_lod_render_cost_units(
		counts, candidate.drawMode, 1);
	    const double nextError = candidate.mesh->projectedErrorAtCut(
		nextCut, candidate.projectedPixelDiameter);
	    if (nextCost <= finalCosts[candidateIndex] &&
		!(nextError < currentError))
		continue;
	    const size_t addedCost = nextCost > finalCosts[candidateIndex] ?
		nextCost - finalCosts[candidateIndex] : 0;
	    const double weightedError = currentError * candidate.errorWeight;
	    const double weightedBenefit = std::max(0.0,
		(currentError - nextError) * candidate.errorWeight);
	    MarginalUpgrade upgrade;
	    upgrade.candidateIndex = candidateIndex;
	    upgrade.nextCut = nextCut;
	    upgrade.nextCost = nextCost;
	    upgrade.addedCost = addedCost;
	    upgrade.qualityFloorViolation =
		candidate.protectedCut > currentCut;
	    upgrade.weightedError = std::isfinite(weightedError) ?
		weightedError : std::numeric_limits<double>::max();
	    upgrade.valuePerCost = addedCost > 0 ?
		weightedBenefit / static_cast<double>(addedCost) :
		weightedBenefit + 1.0;
	    upgrades.push(upgrade);
	    return;
	}
    };

    if (finalCost < effectiveBudget) {
	for (size_t i = 0; i < candidates.size(); ++i)
	    queueNextUpgrade(i);
	while (!upgrades.empty() && finalCost < effectiveBudget) {
	    const MarginalUpgrade upgrade = upgrades.top();
	    upgrades.pop();
	    if (upgrade.addedCost > effectiveBudget - finalCost)
		continue;
	    finalCost += upgrade.addedCost;
	    finalCuts[upgrade.candidateIndex] = upgrade.nextCut;
	    finalCosts[upgrade.candidateIndex] = upgrade.nextCost;
	    queueNextUpgrade(upgrade.candidateIndex);
	}
    }

    double realizedCeiling = 1.0;
    for (size_t i = 0; i < candidates.size(); ++i) {
	const Candidate &candidate = candidates[i];
	const double error = candidate.mesh->projectedErrorAtCut(
	    finalCuts[i], candidate.projectedPixelDiameter) *
	    candidate.errorWeight;
	if (std::isfinite(error))
	    realizedCeiling = std::max(realizedCeiling, error);
	(void)viewState->setCadAllocatedCut(candidate.payload, finalCuts[i],
	    viewRevision, policyRevision, candidate.drawMode);
    }
    return realizedCeiling;
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
    maxProviderMicroseconds(4000),
    forceTerminalLodRefinement(FALSE)
{
}

BObolProgressiveStatus::BObolProgressiveStatus(void)
{
    this->clear();
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
    this->residentGrowthReallocationPending = FALSE;
    this->publicationFramePending = FALSE;
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

/* Camera-rig directions are expressed in eye space (camera looks down -Z,
 * +X right, +Y up).  Directional-light vectors describe photon travel, so a
 * positive X component places the source to the viewer's left.  The studio
 * rig is intentionally asymmetric: an equal ring behaves like ambient light
 * and erases the shape contrast this policy is meant to recover. */
static SbVec3f
bobol_headlight_default_offset(void)
{
    SbVec3f direction(0.35f, -0.25f, -1.0f);
    direction.normalize();
    return direction;
}

static SbVec3f
bobol_mged_headlight_offset(void)
{
    return SbVec3f(0.0f, 0.0f, -1.0f);
}

static SbVec3f
bobol_studio_fill_offset(void)
{
    SbVec3f direction(-0.45f, 0.15f, -1.0f);
    direction.normalize();
    return direction;
}

static SbVec3f
bobol_studio_rim_offset(void)
{
    SbVec3f direction(-0.25f, -0.35f, 1.0f);
    direction.normalize();
    return direction;
}

static double
controller_aspect_from_region(const SbViewportRegion &region)
{
    SbVec2s window = region.getWindowSize();
    if (window[0] <= 0 || window[1] <= 0)
	return 0.0;

    return static_cast<double>(window[0]) / static_cast<double>(window[1]);
}

static const char *
controller_render_environment_name(void)
{
    return "BObolRenderEnvironment";
}

static const char *
controller_headlight_name(void)
{
    return "BObolHeadlight";
}

static const char *
controller_studio_fill_name(void)
{
    return "BObolStudioFill";
}

static const char *
controller_studio_rim_name(void)
{
    return "BObolStudioRim";
}

static const char *
controller_clip_plane_name(SbBool minimum)
{
    return minimum ? "BObolClipMinimum" : "BObolClipMaximum";
}

static SoGroup *
controller_find_render_environment(SoSeparator *root)
{
    if (!root)
	return NULL;

    const char *name = controller_render_environment_name();
    for (int i = 0; i < root->getNumChildren(); i++) {
	SoNode *child = root->getChild(i);
	if (child &&
	    child->isOfType(SoGroup::getClassTypeId()) &&
	    bu_strcmp(child->getName().getString(), name) == 0)
	    return static_cast<SoGroup *>(child);
    }

    return NULL;
}

static SoDirectionalLight *
controller_find_camera_light(SoSeparator *root, const char *name)
{
    if (!root || !name)
	return NULL;
    for (int i = 0; i < root->getNumChildren(); i++) {
	SoNode *child = root->getChild(i);
	if (child && child->isOfType(SoDirectionalLight::getClassTypeId()) &&
	    bu_strcmp(child->getName().getString(), name) == 0)
	    return static_cast<SoDirectionalLight *>(child);
    }
    return NULL;
}

/* A fixed-function OpenGL light is transformed into eye space when its
 * direction is submitted, not when geometry is drawn.  Keep all camera-rig
 * fields in world space for shader/RT consumers and traverse them immediately
 * after the camera in deterministic key/fill/rim order. */
static void
controller_place_camera_rig(SoViewport *viewport,
			    SoDirectionalLight *key,
			    SoDirectionalLight *fill,
			    SoDirectionalLight *rim)
{
    if (!viewport || !viewport->getRoot() || !key || !fill || !rim)
	return;

    SoSeparator *root = viewport->getRoot();
    SoDirectionalLight *lights[3] = {key, fill, rim};
    SbBool retained[3] = {FALSE, FALSE, FALSE};
    for (size_t i = 0; i < 3; i++) {
	const int oldIndex = root->findChild(lights[i]);
	if (oldIndex >= 0) {
	    lights[i]->ref();
	    retained[i] = TRUE;
	    root->removeChild(oldIndex);
	}
    }

    int cameraIndex = -1;
    for (int i = 0; i < root->getNumChildren(); i++) {
	SoNode *child = root->getChild(i);
	if (child && child->isOfType(SoCamera::getClassTypeId())) {
	    cameraIndex = i;
	    break;
	}
    }
    int insertAt = cameraIndex >= 0 ? cameraIndex + 1 :
	root->getNumChildren();
    for (size_t i = 0; i < 3; i++) {
	root->insertChild(lights[i], insertAt++);
	if (retained[i])
	    lights[i]->unref();
    }
}

static void
controller_ensure_camera_rig(SoViewport *viewport,
			     SoDirectionalLight *migratedKey)
{
    if (!viewport || !viewport->getRoot())
	return;
    SoSeparator *root = viewport->getRoot();
    SoDirectionalLight *key = controller_find_camera_light(root,
	controller_headlight_name());
    if (!key) {
	key = migratedKey ? migratedKey : new SoDirectionalLight;
	key->color = SbColor(1.0f, 1.0f, 1.0f);
	key->intensity = 0.68f;
	key->direction = bobol_headlight_default_offset();
    }
    key->setName(SbName(controller_headlight_name()));

    SoDirectionalLight *fill = controller_find_camera_light(root,
	controller_studio_fill_name());
    if (!fill) {
	fill = new SoDirectionalLight;
	fill->setName(SbName(controller_studio_fill_name()));
	fill->color = SbColor(1.0f, 1.0f, 1.0f);
	fill->intensity = 0.22f;
	fill->direction = bobol_studio_fill_offset();
    }
    SoDirectionalLight *rim = controller_find_camera_light(root,
	controller_studio_rim_name());
    if (!rim) {
	rim = new SoDirectionalLight;
	rim->setName(SbName(controller_studio_rim_name()));
	rim->color = SbColor(1.0f, 1.0f, 1.0f);
	rim->intensity = 0.18f;
	rim->direction = bobol_studio_rim_offset();
    }
    controller_place_camera_rig(viewport, key, fill, rim);
}

static void
controller_configure_render_environment(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return;

    SoSeparator *root = viewport->getRoot();
    SoGroup *renderEnvironment = controller_find_render_environment(root);
    if (renderEnvironment) {
	int hasDepthBuffer = 0;
	int hasEnvironment = 0;
	int hasLightModel = 0;
	int hasClipMinimum = 0;
	int hasClipMaximum = 0;
	SoDirectionalLight *legacyHeadlight = NULL;
	for (int i = 0; i < renderEnvironment->getNumChildren(); i++) {
	    SoNode *child = renderEnvironment->getChild(i);
	    if (child && child->isOfType(SoDepthBuffer::getClassTypeId())) {
		hasDepthBuffer = 1;
	    }
	    if (child && child->isOfType(SoEnvironment::getClassTypeId()))
		hasEnvironment = 1;
	    if (child && child->isOfType(SoLightModel::getClassTypeId()))
		hasLightModel = 1;
	    if (child && child->isOfType(SoClipPlane::getClassTypeId())) {
		const char *name = child->getName().getString();
		hasClipMinimum |= bu_strcmp(name,
		    controller_clip_plane_name(TRUE)) == 0;
		hasClipMaximum |= bu_strcmp(name,
		    controller_clip_plane_name(FALSE)) == 0;
	    }
	    if (child &&
		child->isOfType(SoDirectionalLight::getClassTypeId()))
		legacyHeadlight = static_cast<SoDirectionalLight *>(child);
	}
	if (!hasDepthBuffer)
	    renderEnvironment->insertChild(new SoDepthBuffer, 0);
	/* A retained scene created by an older controller may have only the depth
	 * and clip nodes.  Complete it in place so named lighting profiles cannot
	 * silently lose their ambient or PHONG policy after attachment migration. */
	if (!hasEnvironment) {
	    SoEnvironment *environment = new SoEnvironment;
	    environment->ambientIntensity = 0.18f;
	    environment->ambientColor = SbColor(1.0f, 1.0f, 1.0f);
	    renderEnvironment->addChild(environment);
	}
	if (!hasLightModel) {
	    SoLightModel *lightModel = new SoLightModel;
	    lightModel->model = SoLightModel::PHONG;
	    renderEnvironment->addChild(lightModel);
	}
	if (!hasClipMinimum) {
	    SoClipPlane *plane = new SoClipPlane;
	    plane->setName(SbName(controller_clip_plane_name(TRUE)));
	    plane->on = FALSE;
	    renderEnvironment->addChild(plane);
	}
	if (!hasClipMaximum) {
	    SoClipPlane *plane = new SoClipPlane;
	    plane->setName(SbName(controller_clip_plane_name(FALSE)));
	    plane->on = FALSE;
	    renderEnvironment->addChild(plane);
	}
	if (legacyHeadlight) {
	    legacyHeadlight->ref();
	    renderEnvironment->removeChild(legacyHeadlight);
	    legacyHeadlight->setName(SbName(controller_headlight_name()));
	}
	const int index = root->findChild(renderEnvironment);
	if (index > 0) {
	    renderEnvironment->ref();
	    root->removeChild(renderEnvironment);
	    root->insertChild(renderEnvironment, 0);
	    renderEnvironment->unref();
	}
	controller_ensure_camera_rig(viewport, legacyHeadlight);
	if (legacyHeadlight)
	    legacyHeadlight->unref();
	return;
    }

    renderEnvironment = new SoGroup;
    renderEnvironment->setName(SbName(controller_render_environment_name()));

    SoDepthBuffer *depthBuffer = new SoDepthBuffer;
    depthBuffer->test = TRUE;
    depthBuffer->write = TRUE;
    renderEnvironment->addChild(depthBuffer);

    SoEnvironment *environment = new SoEnvironment;
    environment->ambientIntensity = 0.18f;
    environment->ambientColor = SbColor(1.0f, 1.0f, 1.0f);
    renderEnvironment->addChild(environment);

    SoLightModel *lightModel = new SoLightModel;
    lightModel->model = SoLightModel::PHONG;
    renderEnvironment->addChild(lightModel);

    SoDirectionalLight *headlight = new SoDirectionalLight;
    headlight->setName(SbName(controller_headlight_name()));
    headlight->color = SbColor(1.0f, 1.0f, 1.0f);
    headlight->intensity = 0.68f;
    headlight->direction = bobol_headlight_default_offset();

    SoClipPlane *clipMinimum = new SoClipPlane;
    clipMinimum->setName(SbName(controller_clip_plane_name(TRUE)));
    clipMinimum->on = FALSE;
    renderEnvironment->addChild(clipMinimum);

    SoClipPlane *clipMaximum = new SoClipPlane;
    clipMaximum->setName(SbName(controller_clip_plane_name(FALSE)));
    clipMaximum->on = FALSE;
    renderEnvironment->addChild(clipMaximum);

    root->insertChild(renderEnvironment, 0);
    controller_ensure_camera_rig(viewport, headlight);
}

static SoClipPlane *
controller_clip_plane(SoViewport *viewport, SbBool minimum)
{
    if (!viewport || !viewport->getRoot())
	return NULL;
    controller_configure_render_environment(viewport);
    SoGroup *environment =
	controller_find_render_environment(viewport->getRoot());
    if (!environment)
	return NULL;
    const char *wanted = controller_clip_plane_name(minimum);
    for (int i = 0; i < environment->getNumChildren(); i++) {
	SoNode *child = environment->getChild(i);
	if (child && child->isOfType(SoClipPlane::getClassTypeId()) &&
	    bu_strcmp(child->getName().getString(), wanted) == 0)
	    return static_cast<SoClipPlane *>(child);
    }
    return NULL;
}

static SoDepthBuffer *
controller_depth_buffer(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return NULL;
    controller_configure_render_environment(viewport);
    SoGroup *environment =
	controller_find_render_environment(viewport->getRoot());
    if (!environment)
	return NULL;
    for (int i = 0; i < environment->getNumChildren(); i++) {
	SoNode *child = environment->getChild(i);
	if (child && child->isOfType(SoDepthBuffer::getClassTypeId()))
	    return static_cast<SoDepthBuffer *>(child);
    }
    return NULL;
}

static SoLightModel *
controller_light_model(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return NULL;
    controller_configure_render_environment(viewport);
    SoGroup *environment =
	controller_find_render_environment(viewport->getRoot());
    if (!environment)
	return NULL;
    for (int i = 0; i < environment->getNumChildren(); i++) {
	SoNode *child = environment->getChild(i);
	if (child && child->isOfType(SoLightModel::getClassTypeId()))
	    return static_cast<SoLightModel *>(child);
    }
    return NULL;
}

static SoEnvironment *
controller_environment(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return NULL;
    controller_configure_render_environment(viewport);
    SoGroup *renderEnvironment =
	controller_find_render_environment(viewport->getRoot());
    if (!renderEnvironment)
	return NULL;
    for (int i = 0; i < renderEnvironment->getNumChildren(); i++) {
	SoNode *child = renderEnvironment->getChild(i);
	if (child && child->isOfType(SoEnvironment::getClassTypeId()))
	    return static_cast<SoEnvironment *>(child);
    }
    return NULL;
}

static SoDirectionalLight *
controller_headlight(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return NULL;
    controller_configure_render_environment(viewport);
    return controller_find_camera_light(viewport->getRoot(),
	controller_headlight_name());
}

static SoDirectionalLight *
controller_studio_fill(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return NULL;
    controller_configure_render_environment(viewport);
    return controller_find_camera_light(viewport->getRoot(),
	controller_studio_fill_name());
}

static SoDirectionalLight *
controller_studio_rim(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return NULL;
    controller_configure_render_environment(viewport);
    return controller_find_camera_light(viewport->getRoot(),
	controller_studio_rim_name());
}

static const char *
controller_scene_lights_group_name(void)
{
    return "BObolSceneLights";
}

/* Locate (creating if needed) the in-scene lights group, always positioned
 * after the camera rig in the viewport root so fixed-function GL
 * transforms their world-space positions/directions into eye space when the
 * light nodes are traversed. */
static SoGroup *
controller_scene_lights_group(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return NULL;
    controller_configure_render_environment(viewport);
    SoSeparator *root = viewport->getRoot();

    SoGroup *group = NULL;
    for (int i = 0; i < root->getNumChildren(); i++) {
	SoNode *child = root->getChild(i);
	if (child && child->isOfType(SoGroup::getClassTypeId()) &&
	    bu_strcmp(child->getName().getString(),
		controller_scene_lights_group_name()) == 0) {
	    group = static_cast<SoGroup *>(child);
	    break;
	}
    }

    const SbBool existed = (group != NULL);
    if (!group) {
	group = new SoGroup;
	group->setName(SbName(controller_scene_lights_group_name()));
    } else {
	group->ref();
	root->removeChild(group);
    }

    int cameraIndex = -1;
    for (int i = 0; i < root->getNumChildren(); i++) {
	if (root->getChild(i) &&
	    root->getChild(i)->isOfType(SoCamera::getClassTypeId())) {
	    cameraIndex = i;
	    break;
	}
    }
    int insertAt = (cameraIndex >= 0) ? cameraIndex + 1 :
	root->getNumChildren();
    SoDirectionalLight *cameraLights[3] = {
	controller_find_camera_light(root, controller_headlight_name()),
	controller_find_camera_light(root, controller_studio_fill_name()),
	controller_find_camera_light(root, controller_studio_rim_name())
    };
    for (size_t i = 0; i < 3; i++) {
	const int lightIndex = cameraLights[i] ?
	    root->findChild(cameraLights[i]) : -1;
	if (lightIndex >= insertAt)
	    insertAt = lightIndex + 1;
    }
    root->insertChild(group, insertAt);
    if (existed)
	group->unref();
    return group;
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

static std::vector<SoBRLDatabaseSource *>
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

static std::vector<SoBRLDatabaseSource *>
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

/**
 * Owner-thread progressive-display state.  This first extraction preserves
 * the exact field layout and algorithms formerly embedded in Impl; it adds
 * no allocation, virtual dispatch, or per-occurrence indirection.
 */
struct BObolLodViewSnapshot {
    bool same(const BObolLodViewSnapshot &other) const
    {
	return this->haveCamera == other.haveCamera &&
	    this->width == other.width &&
	    this->height == other.height &&
	    memcmp(&this->size, &other.size, sizeof(this->size)) == 0 &&
	    memcmp(&this->lodScale, &other.lodScale,
		sizeof(this->lodScale)) == 0 &&
	    memcmp(&this->curveScale, &other.curveScale,
		sizeof(this->curveScale)) == 0 &&
	    memcmp(&this->pointScale, &other.pointScale,
		sizeof(this->pointScale)) == 0 &&
	    this->botThreshold == other.botThreshold &&
	    memcmp(this->viewVolumeMatrix, other.viewVolumeMatrix,
		sizeof(this->viewVolumeMatrix)) == 0;
    }

    uint8_t haveCamera = 0;
    int32_t width = 0;
    int32_t height = 0;
    double size = 0.0;
    double lodScale = 0.0;
    double curveScale = 0.0;
    double pointScale = 0.0;
    uint32_t botThreshold = 0;
    float viewVolumeMatrix[16] = {};
};

struct BObolLodViewScaleSnapshot {
    bool same(const BObolLodViewScaleSnapshot &other) const
    {
	return this->haveCamera == other.haveCamera &&
	    this->width == other.width &&
	    this->height == other.height &&
	    memcmp(&this->size, &other.size, sizeof(this->size)) == 0 &&
	    memcmp(&this->lodScale, &other.lodScale,
		sizeof(this->lodScale)) == 0 &&
	    memcmp(&this->curveScale, &other.curveScale,
		sizeof(this->curveScale)) == 0 &&
	    memcmp(&this->pointScale, &other.pointScale,
		sizeof(this->pointScale)) == 0 &&
	    this->botThreshold == other.botThreshold &&
	    this->viewportWidth == other.viewportWidth &&
	    this->viewportHeight == other.viewportHeight &&
	    this->cameraTypeKey == other.cameraTypeKey &&
	    memcmp(&this->aspectRatio, &other.aspectRatio,
		sizeof(this->aspectRatio)) == 0 &&
	    memcmp(&this->focalDistance, &other.focalDistance,
		sizeof(this->focalDistance)) == 0 &&
	    memcmp(&this->projectionScale, &other.projectionScale,
		sizeof(this->projectionScale)) == 0;
    }

    uint8_t haveCamera = 0;
    int32_t width = 0;
    int32_t height = 0;
    double size = 0.0;
    double lodScale = 0.0;
    double curveScale = 0.0;
    double pointScale = 0.0;
    uint32_t botThreshold = 0;
    int16_t viewportWidth = 0;
    int16_t viewportHeight = 0;
    uint64_t cameraTypeKey = 0;
    float aspectRatio = 0.0f;
    float focalDistance = 0.0f;
    float projectionScale = 0.0f;
};

static_assert(std::is_trivially_copyable<BObolLodViewSnapshot>::value,
    "view signatures must remain allocation-free values");
static_assert(std::is_trivially_copyable<BObolLodViewScaleSnapshot>::value,
    "view scale signatures must remain allocation-free values");

struct BObolLodCoordinator : BObolLodStateMachine {
    float quietAllocationTargetFps(void) const
    {
	/* The ordinary quiet cadence keeps cold streaming and multi-frame
	 * convergence responsive.  Once that phase has produced a terminal exact
	 * frame, an event-driven view may spend up to the separate hard static
	 * presentation deadline on a richer framebuffer which is then retained
	 * without redraw. */
	if (!lodStaticOverscanActive ||
	    !stablePresentationFrameDeadlineNanoseconds)
	    return lodStableTargetFps;
	const long double deadlineFps = 1000000000.0L /
	    static_cast<long double>(
		stablePresentationFrameDeadlineNanoseconds);
	if (!std::isfinite(deadlineFps) || deadlineFps <= 0.0L)
	    return lodStableTargetFps;
	const float bounded = static_cast<float>(deadlineFps);
	return lodStableTargetFps > 0.0f ?
	    std::min(lodStableTargetFps, bounded) : bounded;
    }

    void reconcilePhase(Event event, size_t pendingCount = 0)
    {
	Inputs inputs;
	inputs.interactive = lodInteractive;
	inputs.gestureActive = lodGestureActive;
	/* Coverage is a prerequisite only while automatic mesh LoD is an
	 * active consumer.  A disabled/removed service has no mesh inventory to
	 * prove, and a planning cursor from the preceding source epoch must not
	 * be projected as concurrently compacting an invalidated inventory. */
	inputs.coverageComplete = lodCoveragePolicy.effectiveComplete();
	inputs.coverageActive = lodCoveragePolicy.effectiveActive();
	inputs.compacting = lodCompactionPolicy.planning() &&
	    lodCoveragePolicy.compactionAllowed();
	inputs.cpuMemoryPressure = lodResourcePolicy.cpuPressure();
	inputs.gpuMemoryPressure = lodResourcePolicy.gpuPressure();
	inputs.resourceRecoveryPending =
	    lodResourcePolicy.recoveryPending();
	inputs.generationActive = lodActiveGeneration != 0;
	inputs.settlingWork =
	    lodSubmissionPending || lodSubmissionRescanPending ||
	    lodRetainedRefinementPending ||
	    lodRetainedResidencyPending ||
	    lodBudgetPolicy.rescanAfterFrame() ||
	    lodPresentationPolicy.handoffPending() ||
	    lodRefinementAwaitingFrame ||
	    lodPublicationPolicy.pending() || pendingCount != 0 ||
	    lodBudgetPolicy.calibrationFramesRemaining() != 0 ||
	    lodStablePointProxyCalibrationPending ||
	    lodPointProxyTriangleRecoveryPending ||
	    lodHeadroomPolicy.retryPending();

	Witness witness;
	witness.viewEpoch = lodViewRevision.value();
	witness.policyEpoch = lodPolicyRevision.value();
	witness.renderSerial = renderCompletionSerial;
	witness.activeGeneration = lodActiveGeneration;
	witness.residentDemandRevision =
	    lodCompactionPolicy.demandRevision();
	witness.resourcePressureRevision = lodResourcePolicy.revision();
	witness.visibleCount = lodCoveragePolicy.visibleCount();
	witness.completedCount = lodCoveragePolicy.coveredCount();
	witness.pendingCount = pendingCount;
	(void)this->dispatch(event, inputs, witness);
    }

    uint64_t lastRenderTimeNanoseconds = 0;
    uint64_t smoothedRenderTimeNanoseconds = 0;
    uint64_t lastBackgroundRenderTimeNanoseconds = 0;
    uint64_t lastSceneRenderTimeNanoseconds = 0;
    /* Full CAD plan/atlas construction is a real latency problem, but not a
     * measurement of steady triangle throughput.  Only an unchanged prepared
     * replay may drive quiet-view capacity cuts. */
    SbBool lodLastRenderWasPreparedCadReplay = TRUE;
    SbBool lodLastRenderWasReusableCadPresentation = TRUE;
    uint64_t lodLastCadGpuSampleSerial = 0;
    uint64_t lodLastCadGpuTimeNanoseconds = 0;
    uint64_t lodLastCadGpuGeometryUploadBytes = 0;
    SbBool lodLastCadGpuGeometryUploadBytesValid = FALSE;
    uint64_t lastProgressiveAdvanceTimeNanoseconds = 0;
    uint64_t lastLodResultProcessingTimeNanoseconds = 0;
    uint64_t lastProgressiveProviderTimeNanoseconds = 0;
    uint64_t lastLodSubmissionTimeNanoseconds = 0;
    uint64_t lastPresentationSyncTimeNanoseconds = 0;
    uint64_t interactivePresentationFrameDeadlineNanoseconds = 40000000ULL;
    uint64_t stablePresentationFrameDeadlineNanoseconds = 100000000ULL;
    uint64_t interruptedPresentationFrameCount = 0;
    uint64_t consecutiveInterruptedPresentationFrames = 0;
    uint64_t lastInterruptedPresentationTimeNanoseconds = 0;
    uint64_t renderCompletionSerial = 0;
    mutable std::mutex presentationTimingMutex;
    uint64_t presentedFrameSerial = 0;
    uint64_t lastPresentationTimestampNanoseconds = 0;
    uint64_t smoothedPresentationIntervalNanoseconds = 0;
    uint64_t displayedPresentationIntervalNanoseconds = 0;
    uint64_t lastInteractivePresentationTimestampNanoseconds = 0;
    uint64_t smoothedInteractivePresentationIntervalNanoseconds = 0;
    std::vector<BObolProgressiveProviderRecord> progressiveProviders;
    uint64_t progressiveProviderNextToken = 1;
    std::atomic<int> progressiveWorkPending {0};
    BObolProgressiveOptions defaultProgressiveOptions;
    BObolLodService *lodService = NULL;
    std::shared_ptr<BObolLodService> managedLodService;
    size_t managedLodWorkerCount = 0;
    uint64_t lodResultSubscriberId = 0;
    std::atomic<int> lodResultsPending {0};
    /* Timestamp the empty->non-empty worker-result edge.  Bulk warm-cache
     * streams use it to coalesce tiny timer-tick dribbles into bounded-latency
     * scene updates. */
    std::atomic<int64_t> lodResultsFirstReadyMicroseconds {0};
    SbBool lodAutoSubmit = FALSE;
    uint64_t lodActiveGeneration = 0;
    size_t lodSubmissionSourceIndex = 0;
    size_t lodSubmissionEntryOffset = 0;
    SoBRLDatabaseSource *lodSubmissionPlanSource = NULL;
    std::vector<size_t> lodSubmissionPlanEntries;
    SbBool lodSubmissionPlanValid = FALSE;
    SbBool lodSubmissionPlanRetainedAdmission = FALSE;
    SbBool lodSubmissionRetainedAdmissionMode = FALSE;
    double lodRetainedAdmissionMaximumNormalizedError =
	std::numeric_limits<double>::infinity();
    SoBRLDatabaseSource *lodSubmissionVisibleCountSource = NULL;
    size_t lodSubmissionVisibleCount = 0;
    /* Large compact scenes are consumed in bounded GUI-thread windows.  A
     * coverage pass proves that every projected leaf has a useful structural
     * presentation before any of those leaves opens a PoP hierarchy.  The
     * following quality pass promotes view-significant leaves under the
     * aggregate budget, preventing early entries from starving later ones. */
    BObolLodCoveragePolicy lodCoveragePolicy;
    /*
     * A camera-scale refresh is not cold coverage.  Existing occurrence cuts
     * remain useful and a zoom must be allowed to retarget them, consume a
     * richer resident prefix, or request the missing cumulative suffix while
     * the bounded visibility scan continues.  New/source inventory coverage
     * keeps the stricter minimum-mesh-first rule below.
     */
    /*
     * The newest completed camera pass proved that every projected mesh
     * occurrence already has a retained mesh presentation.  While input is
     * active, richer worker results may then remain queued without exposing
     * boxes or missing newly visible geometry.
     */
    std::vector<std::pair<SoBRLDatabaseSource *, size_t>>
	lodConvergenceCandidateCounts;
    /*
     * Authoritative visible-mesh denominator from the newest completed
     * all-source coverage pass.  Per-source candidate counts are useful while
     * a bounded pass is still being assembled, but its final window used to
     * clear the accumulated coverage counters before convergence retained
     * their total.  A quality-only policy revision preserves this view proof;
     * camera and source-inventory revisions invalidate it.
    */
    mutable BObolLodConvergencePolicy lodConvergencePolicy;
    void clearLodSubmissionPlan(void)
    {
	lodSubmissionPlanSource = NULL;
	lodSubmissionPlanEntries.clear();
	lodSubmissionPlanValid = FALSE;
	lodSubmissionPlanRetainedAdmission = FALSE;
	lodSubmissionVisibleCountSource = NULL;
	lodSubmissionVisibleCount = 0;
    }
    void clearLodConvergenceCandidates(void)
    {
	lodConvergenceCandidateCounts.clear();
	lodCoveragePolicy.clearCompleteVisibleCount();
    }
    void resetLodConvergenceFraction(void)
    {
	lodConvergencePolicy.resetFraction();
    }
    void setLodConvergenceCandidateCount(
	SoBRLDatabaseSource *source, size_t count)
    {
	if (!source)
	    return;
	for (std::pair<SoBRLDatabaseSource *, size_t> &entry :
	    lodConvergenceCandidateCounts) {
	    if (entry.first != source)
		continue;
	    entry.second = count;
	    return;
	}
	lodConvergenceCandidateCounts.push_back(
	    std::make_pair(source, count));
    }
    size_t lodConvergenceCandidateCount(void) const
    {
	if (lodCoveragePolicy.hasCompleteVisibleCount())
	    return lodCoveragePolicy.completeVisibleCount();
	size_t total = 0;
	for (const std::pair<SoBRLDatabaseSource *, size_t> &entry :
	    lodConvergenceCandidateCounts)
	    total = entry.second > SIZE_MAX - total ?
		SIZE_MAX : total + entry.second;
	return total;
    }
    SbBool lodSubmissionPending = FALSE;
    SbBool lodSubmissionRescanPending = FALSE;
    SbBool lodRetainedRefinementPending = FALSE;
    /* A retained minimax pass selected a level whose PoP suffix is not yet
     * resident.  This is provider work, not a performance-limited quality
     * observation, and must survive the coherent cut presentation barrier. */
    SbBool lodRetainedResidencyPending = FALSE;
    SbBool lodRetainedRefinementCutAdvanced = FALSE;
    SbBool lodRetainedRefinementBudgetBlocked = FALSE;
    SbBool lodPassAdmittedWork = FALSE;
    SbBool forceTerminalLodRefinement = FALSE;
    SbBool lodSubmissionRefreshMissing = TRUE;
    int lodSubmissionReset = 0;
    BObolLodViewEpoch lodLastSubmittedViewRevision;
    BObolLodPolicyEpoch lodLastSubmittedPolicyRevision;
    struct LodSourceSnapshot {
	SoBRLDatabaseSource *source = NULL;
	struct db_i *database = NULL;
	BObolLodSourceRoutingId routingId;
	BObolLodInventoryEpoch inventoryRevision;
	SbString databaseId;
	SbString path;
	int drawMode = 0;
	int representationMode = 0;
	SbBool visible = FALSE;
	int lodBotThreshold = 0;
	uint32_t sourceRevision = 0;
	uint32_t inputsRevision = 0;

	bool sameIdentity(const LodSourceSnapshot &other) const
	{
	    return this->database == other.database &&
		this->routingId == other.routingId &&
		this->databaseId == other.databaseId &&
		this->path == other.path &&
		this->drawMode == other.drawMode &&
		this->representationMode == other.representationMode &&
		this->visible == other.visible &&
		this->lodBotThreshold == other.lodBotThreshold &&
		this->sourceRevision == other.sourceRevision &&
		this->inputsRevision == other.inputsRevision;
	}
    };
    std::vector<LodSourceSnapshot> lodLastSubmittedSources;
    SbBool lodSubmissionDeltaActive = FALSE;
    std::vector<SoBRLDatabaseSource *> lodSubmissionDeltaSources;
    std::vector<std::pair<SoBRLDatabaseSource *, std::vector<size_t>>>
	lodSubmissionDeltaPlans;
    BObolLodViewSnapshot lodViewSignature;
    BObolLodViewScaleSnapshot lodViewScaleSignature;
    SbBool lodViewSignatureValid = FALSE;

    /* Net scale, not the mere presence of a wheel event, decides whether a
     * quiet orthographic view must retarget every occurrence.  Keep the
     * signature and ready-view proof from the start of the continuous input
     * epoch.  A bracketed pose gesture may begin before an unbracketed wheel
     * epoch's debounce expires; in that case these values deliberately keep
     * describing the original stable view. */
    BObolLodViewScaleSnapshot lodInteractionStartScaleSignature;
    SbBool lodInteractionStartScaleSignatureValid = FALSE;
    SbBool lodInteractionStartedFromReadyView = FALSE;
    size_t lodStableBudgetBeforeInteraction = 0;
    SbBool lodStableBudgetBeforeInteractionValid = FALSE;
    /* A completed, fully covered orthographic pose may verify visibility,
     * selection, and resource pressure, but may not rewrite retained PoP
     * cuts merely because the interaction debounce ended. */
    SbBool lodRetainPoseOccurrenceCuts = FALSE;
    /* A settled retained camera epoch first records exact projected demand
     * for every visible occurrence, then performs one scene-budgeted
     * importance reallocation.  Keeping this as a census-completion latch
     * prevents stale previous-view metrics and makes the redistribution
     * explicitly one-shot. */
    SbBool lodRetainedImportanceCensusPending = FALSE;
    /* Worker completions extend immutable residency, but active occurrence
     * cuts remain an aggregate scene-budget decision.  Coalesce a whole
     * completion wave before invoking that allocator. */
    BObolLodResidentGrowthPolicy lodResidentGrowthPolicy;
    BObolLodViewEpoch lodViewRevision {1};
    BObolLodPolicyEpoch lodPolicyRevision {1};
    BObolLodViewDemandPolicy lodViewDemandPolicy;
    /* One scale-quality or exact stable-headroom probe may measure one
     * otherwise-unaffordable populated PoP transition across the entire
     * scene.  Individual submit actions are time-sliced and source-local, so
     * uniqueness must live here. */
    SbBool lodDiscretePopulationTrialAvailable = FALSE;
    int64_t lodLastViewChangeMicroseconds = 0;
    SbBool lodInteractive = FALSE;
    SbBool lodGestureActive = FALSE;
    /* A scale gesture whose camera has paused is still interactive, but it
     * should not remain a magnified copy of the motion cut until button-up.
     * One stable-budget quality probe may expose already resident detail;
     * the next camera event immediately returns to motion calibration. */
    /* A scale-changing frame may be followed by one bounded quality frame.
     * This is distinct from quiet/stable refinement: continuous wheel or
     * trackpad input must be able to expose a newly resident PoP suffix
     * without first ending the interaction. */
    /* Keep the complete held-button presentation through button-release
     * debounce.  This preserves the measured aggregate presentation, not the
     * richer retained per-occurrence cuts which may be hidden by the
     * renderer-wide motion ceiling. */
    SbBool lodReleaseCutFloorActive = FALSE;
    /* Renderer-only presentation continuity and motion-to-stable handoff.
     * The policy owns its snapshot/latch proof; this coordinator only applies
     * returned ceiling and point-proxy values to the retained view state. */
    BObolLodPresentationPolicy lodPresentationPolicy;
    uint64_t lodSettleAfterRenderSerial = 0;
    BObolLodCompactionPolicy lodCompactionPolicy;
    /* Complete-frame resource pressure is retained independently of HUD
     * queries.  A pressure edge requests one safe compaction pass; if the
     * current visible working set itself exceeds a renderer limit, that pass
     * may finish in a stable memory-limited presentation rather than loop. */
    BObolLodResourcePolicy lodResourcePolicy;
    uint64_t lodResidentAdmissionRevision = 0;
    SbBool lodRefinementAwaitingFrame = FALSE;
    uint64_t lodRefinementResumeAfterRenderSerial = 0;
    /* A result-publication frame may include one-time CPU/GPU buffer
     * construction.  It controls refinement pacing, but must not teach the
     * steady-state face-capacity estimator that every later frame costs the
     * same amount. */
    BObolLodPublicationPolicy lodPublicationPolicy;
    int64_t lodRefinementNotBeforeMicroseconds = 0;
    float lodTargetPixelError = 1.0f;
    float lodInteractiveTargetFps = 60.0f;
    float lodStableTargetFps = 20.0f;
    /* A terminal ordinary quiet frame may prove capacity for one or more
     * bounded static-image overscan allocations.  This never changes motion
     * calibration and is reset by every camera, policy, service, or input
     * epoch. */
    SbBool lodStaticOverscanActive = FALSE;
    BObolLodBudgetPolicy lodBudgetPolicy;
    long double lodInteractiveCalibratedRenderCostPerSecond = 0.0L;
    long double lodStableCalibratedRenderCostPerSecond = 0.0L;
    void seedInteractiveCalibrationFromStable(void)
    {
	if (!(lodStableCalibratedRenderCostPerSecond > 0.0L))
	    return;
	/* Stable retained rendering is a measured geometry-throughput witness,
	 * but motion also pays camera/cut/update overhead.  Carry half of that
	 * known capacity into a new gesture.  This is deliberately a floor, not
	 * a reset: an established faster interaction estimate is retained, while
	 * a 118-triangle underloaded frame can no longer become a permanent
	 * self-fulfilling capacity ceiling. */
	const long double seed =
	    lodStableCalibratedRenderCostPerSecond * 0.5L;
	lodInteractiveCalibratedRenderCostPerSecond = std::max(
	    lodInteractiveCalibratedRenderCostPerSecond, seed);
    }
    /* A late reusable frame may demonstrate headroom after the ordinary
     * bounded probe series has ended.  Remember the exact camera/policy and
     * budget already retried so each newly enlarged allowance gets one
     * admission pass without permitting an unchanged discrete-cut loop. */
    BObolLodHeadroomPolicy lodHeadroomPolicy;
    int lodInteractiveProgressiveCeiling = -1;
    /*
     * Current renderer-side small-occurrence aggregation cut.  Interaction
     * raises it immediately when measured frames miss their target.  A quiet
     * performance-limited view may also raise it when every retained mesh is
     * already at its minimum prefix and that irreducible population still
     * exceeds the calibrated scene budget.
     */
    float lodPresentationPointProxyPixelThreshold = 1.0f;
    BObolLodPointProxyCalibrationPolicy lodPointProxyCalibrationPolicy;
    /*
     * A changed quiet aggregation threshold needs one measured presentation
     * before convergence is authoritative.  Keep this distinct from PoP
     * admission: no provider scan or retained-array mutation is required.
     */
    SbBool lodStablePointProxyCalibrationPending = FALSE;
    /* A temporary multi-pixel point cut may be relaxed only after reducible
     * PoP triangle detail has been compacted to the measured scene capacity.
     * This latch keeps convergence behind that retained-prefix pass. */
    SbBool lodPointProxyTriangleRecoveryPending = FALSE;
    uint64_t lodInteractiveCeilingFeedbackRenderSerial = 0;
    SbBool lodUseForcedCut = FALSE;
    int lodForcedCut = 0;
    uint64_t maxExactFullDetailFaceCount = 0;
    uint64_t maxExactFullDetailPointCount = 0;
    std::vector<BObolRtPickCache *> rtPickCaches;
    std::vector<SbString> rtPickCachePaths;
    std::vector<struct db_i *> rtPickCacheDatabases;
    std::vector<uint32_t> rtPickCacheSourceRevisions;
    SbBool meshResidencyBudgetEnabled = FALSE;
    size_t maxResidentMeshBytes = 0;
    SbBool meshResidencyEvictDisplayPayloads = TRUE;
    size_t lastMeshBudgetInitialResidentBytes = 0;
    size_t lastMeshBudgetFinalResidentBytes = 0;
    size_t lastMeshBudgetFreedFullDetailBytes = 0;
    size_t lastMeshBudgetFreedDisplayBytes = 0;
    unsigned int lastMeshBudgetVisitedMeshCount = 0;
    unsigned int lastMeshBudgetEvictedFullDetailMeshCount = 0;
    unsigned int lastMeshBudgetEvictedDisplayMeshCount = 0;
    unsigned int lastLodVisitedMeshCount = 0;
    unsigned int lastLodSubmittedTaskCount = 0;
    unsigned int lastLodUpdatedCutCount = 0;
    unsigned int lastLodSkippedMeshCount = 0;
    size_t lastLodResultCount = 0;
    unsigned int lastLodMatchedResultCount = 0;
    unsigned int lastLodAppliedResultCount = 0;
    unsigned int lastLodRejectedResultCount = 0;
    unsigned int lastLodUnmatchedResultCount = 0;
    SbString lastLodDiagnostics = SbString("");
};

static BObolLodPresentationPolicy::Population
controller_lod_presentation_population(const BObolViewLodState *state)
{
    BObolLodPresentationPolicy::Population population;
    population.available = state != NULL;
    population.activeFaces = state ? state->activeFaceCount() : 0;
    population.cadRevision = state ? state->cadRevision() : 0;
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

struct BObolViewController::Impl : BObolLodCoordinator {
    explicit Impl(BObolViewController *owner) :
	featureStore(new BObolFeatureStore(owner)),
	polygonStore(new BObolPolygonStore(owner)),
	selectionStore(new BObolSelectionStore)
    {
    }

    BObolSceneController sceneController;
    SoViewport *viewport = new SoViewport;
    SoBRLViewLodGroup *renderLodRoot = NULL;
    SoNode *renderBatchRoot = NULL;
    SoNode *renderPresentationRoot = NULL;
    SoGroup *framebufferUnderlayRoot = NULL;
    SoGroup *framebufferInterlayRoot = NULL;
    SoGroup *framebufferOverlayRoot = NULL;
    BObolViewAttachment *viewAttachment = new BObolViewAttachment;
    SoRenderManager *renderManager = new SoRenderManager;
    SoOffscreenRenderer *imageRenderer = NULL;
    SoDB::ContextManager *imageRendererManager = NULL;
    SoCamera *activeCamera = NULL;
    SbViewportRegion viewportRegion = SbViewportRegion(1, 1);
    SbColor backgroundBottom = SbColor(0.0f, 0.0f, 0.0f);
    SbColor backgroundTop = SbColor(0.0f, 0.0f, 0.0f);
    SoftwareWireMode softwareWireMode = SOFTWARE_WIRE_AUTO;
    std::atomic<int> endpointGraphicalRenderingEnabled {1};
    SbBool transparencyEnabled = TRUE;
    SbBool antialiasingEnabled = FALSE;
    /* Camera-driven lighting state.  The selected profile owns ambient level,
     * relative rig intensities, and fill/rim directions; headlightOffsetEye is
     * the independently adjustable primary/key direction. */
    BObolViewController::LightingProfile lightingProfile =
	BObolViewController::LIGHTING_STUDIO;
    SbBool headlightEnabled = TRUE;
    SbBool headlightCameraTracked = TRUE;
    SbBool sceneLightsEnabled = FALSE;
    SbVec3f headlightOffsetEye = bobol_headlight_default_offset();
    SbRotation lastCameraOrientation = SbRotation::identity();
    /* In-scene lights, supplied by the GED layer (cache-immune) rather than the
     * geometry realize walk (which is skipped on a warm LoD cache). */
    std::vector<BObolSceneLightRealization> sceneLights;
    double clipMinimum = BV_VIEW_MIN;
    double clipMaximum = BV_VIEW_MAX;
    mutable std::mutex renderRequestMutex;
    SbBool renderRequested = FALSE;
    /* Coalesced requests are capacity relevant if any constituent request is.
     * This prevents a selection repaint from masking a simultaneous PoP cut
     * publication, while keeping an isolated style patch out of LoD feedback. */
    SbBool renderLodCapacityRelevant = FALSE;
    SbString renderReason = SbString("");
    uint64_t renderRequestSerial = 0;
    std::mutex frameRequestMutex;
    BObolFrameRequestCallback frameRequestCallback = NULL;
    void *frameRequestUserData = NULL;
    mutable std::mutex presentationSyncMutex;
    BObolPresentationSyncCallback presentationSyncCallback = NULL;
    void *presentationSyncUserData = NULL;
    BObolFeatureStore *featureStore;
    BObolPolygonStore *polygonStore;
    BObolSelectionStore *selectionStore;
};

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
BObolViewController::setCamera(SoCamera *camera)
{
    if (camera)
	camera->ref();
    if (this->d->activeCamera)
	this->d->activeCamera->unref();
    this->d->activeCamera = camera;
    this->d->viewport->setCamera(camera);
    controller_configure_render_environment(this->d->viewport);
    /* Keep the in-scene light group positioned immediately after the (possibly
     * new) camera so its world-space light positions stay correct. */
    (void)controller_scene_lights_group(this->d->viewport);
    this->syncRenderManager();
    this->syncLodViewSignature(TRUE);
    this->requestRender("camera");
}

SoCamera *
BObolViewController::getCamera(void) const
{
    return this->d->activeCamera;
}

void
BObolViewController::setViewportRegion(const SbViewportRegion &region)
{
    this->d->viewportRegion = region;
    this->d->viewport->setViewportRegion(region);
    this->syncRenderManager();
    this->syncLodViewSignature(TRUE);
    this->requestRender("viewport");
}

const SbViewportRegion &
BObolViewController::getViewportRegion(void) const
{
    return this->d->viewportRegion;
}

void
BObolViewController::setViewportSize(unsigned int width, unsigned int height)
{
    if (width > 32767)
	width = 32767;
    if (height > 32767)
	height = 32767;
    if (width == 0)
	width = 1;
    if (height == 0)
	height = 1;
    const SbVec2s current = this->d->viewportRegion.getWindowSize();
    const SbVec2s currentOrigin =
	this->d->viewportRegion.getViewportOriginPixels();
    const SbVec2s currentViewport =
	this->d->viewportRegion.getViewportSizePixels();
    if (current[0] == static_cast<short>(width) &&
	current[1] == static_cast<short>(height) &&
	currentOrigin[0] == 0 && currentOrigin[1] == 0 &&
	currentViewport[0] == static_cast<short>(width) &&
	currentViewport[1] == static_cast<short>(height))
	return;
    this->d->viewportRegion.setWindowSize((short)width, (short)height);
    this->d->viewportRegion.setViewportPixels(0, 0, (short)width,
	(short)height);
    this->d->viewport->setViewportRegion(this->d->viewportRegion);
    this->syncRenderManager();
    this->syncLodViewSignature(TRUE);
    this->requestRender("viewport-size");
}

void
BObolViewController::setBackgroundColors(const SbColor &bottom,
	const SbColor &top)
{
    if (this->d->backgroundBottom == bottom && this->d->backgroundTop == top)
	return;
    this->d->backgroundBottom = bottom;
    this->d->backgroundTop = top;
    SoEnvironment *environment = controller_environment(this->d->viewport);
    if (environment)
	environment->fogColor = top;
    this->requestRender("background");
}

const SbColor &
BObolViewController::getBackgroundBottomColor(void) const
{
    return this->d->backgroundBottom;
}

const SbColor &
BObolViewController::getBackgroundTopColor(void) const
{
    return this->d->backgroundTop;
}

void
BObolViewController::setDepthTestEnabled(SbBool enabled)
{
    SoDepthBuffer *depth = controller_depth_buffer(this->d->viewport);
    if (!depth || (depth->test.getValue() == enabled &&
	depth->write.getValue() == enabled))
	return;
    depth->test = enabled;
    depth->write = enabled;
    this->requestRender("depth-test");
}

SbBool
BObolViewController::isDepthTestEnabled(void) const
{
    SoDepthBuffer *depth = controller_depth_buffer(this->d->viewport);
    return depth ? depth->test.getValue() : TRUE;
}

void
BObolViewController::setTransparencyEnabled(SbBool enabled)
{
    enabled = enabled ? TRUE : FALSE;
    if (this->d->transparencyEnabled == enabled)
	return;
    this->d->transparencyEnabled = enabled;
    const SoGLRenderAction::TransparencyType type = enabled ?
	SoGLRenderAction::BLEND : SoGLRenderAction::NONE;
    this->d->renderManager->getGLRenderAction()->setTransparencyType(type);
    if (this->d->imageRenderer)
	this->d->imageRenderer->getGLRenderAction()->setTransparencyType(type);
    this->requestRender("transparency");
}

SbBool
BObolViewController::isTransparencyEnabled(void) const
{
    return this->d->transparencyEnabled;
}

void
BObolViewController::setAntialiasingEnabled(SbBool enabled)
{
    enabled = enabled ? TRUE : FALSE;
    if (this->d->antialiasingEnabled == enabled)
	return;
    this->d->antialiasingEnabled = enabled;
    this->d->renderManager->setAntialiasing(enabled, 1);
    if (this->d->imageRenderer) {
	SoGLRenderAction *action = this->d->imageRenderer->getGLRenderAction();
	if (action) {
	    action->setSmoothing(enabled);
	    action->setNumPasses(1);
	}
    }
    this->requestRender("antialiasing");
}

SbBool
BObolViewController::isAntialiasingEnabled(void) const
{
    return this->d->antialiasingEnabled;
}

SbBool
BObolViewController::setClipBounds(double minimum, double maximum)
{
    if (!std::isfinite(minimum) || !std::isfinite(maximum) ||
	minimum > maximum)
	return FALSE;
    this->d->clipMinimum = minimum;
    this->d->clipMaximum = maximum;
    this->requestRender("clip-bounds");
    return TRUE;
}

void
BObolViewController::getClipBounds(double &minimum, double &maximum) const
{
    minimum = this->d->clipMinimum;
    maximum = this->d->clipMaximum;
}

size_t
BObolViewController::getActiveClipPlanes(SbPlane planes[2]) const
{
    if (!planes)
	return 0;
    size_t count = 0;
    SoClipPlane *minimum = controller_clip_plane(this->d->viewport, TRUE);
    SoClipPlane *maximum = controller_clip_plane(this->d->viewport, FALSE);
    if (minimum && minimum->on.getValue())
	planes[count++] = minimum->plane.getValue();
    if (maximum && maximum->on.getValue())
	planes[count++] = maximum->plane.getValue();
    return count;
}

/* Rewrite the camera rig's world-space directions from the stored camera
 * orientation.  A forced update is used when changing profiles while tracking
 * is disabled: the newly selected rig is aimed once, then remains scene-fixed. */
void
BObolViewController::applyTrackedHeadlight(SbBool force)
{
    if (!this->d->headlightEnabled ||
	(!force && !this->d->headlightCameraTracked))
	return;
    SoDirectionalLight *lights[3] = {
	controller_headlight(this->d->viewport),
	controller_studio_fill(this->d->viewport),
	controller_studio_rim(this->d->viewport)
    };
    if (!lights[0] || !lights[1] || !lights[2])
	return;
    SbVec3f eyeDirections[3] = {
	this->d->headlightOffsetEye,
	bobol_studio_fill_offset(),
	bobol_studio_rim_offset()
    };
    for (size_t i = 0; i < 3; i++) {
	SbVec3f worldDir;
	this->d->lastCameraOrientation.multVec(eyeDirections[i], worldDir);
	if (worldDir.normalize() > 0.0f &&
	    lights[i]->direction.getValue() != worldDir)
	    lights[i]->direction = worldDir;
    }
}

void
BObolViewController::setLightingEnabled(SbBool enabled)
{
    SoLightModel *model = controller_light_model(this->d->viewport);
    SoDirectionalLight *lights[3] = {
	controller_headlight(this->d->viewport),
	controller_studio_fill(this->d->viewport),
	controller_studio_rim(this->d->viewport)
    };
    if (!model || !lights[0] || !lights[1] || !lights[2])
	return;
    const int requested = enabled ? SoLightModel::PHONG :
	SoLightModel::BASE_COLOR;
    /* Master lighting picks PHONG vs flat shading.  The headlight only actually
     * contributes when both master lighting and the per-view headlight toggle
     * are on, so restore the headlight to its finer-grained state rather than
     * unconditionally forcing it on. */
    const SbBool lightOn = enabled && this->d->headlightEnabled;
    const SbBool studio = this->d->lightingProfile == LIGHTING_STUDIO;
    SbBool changed = model->model.getValue() != requested;
    model->model = requested;
    const SbBool desired[3] = {lightOn, lightOn && studio,
	lightOn && studio};
    for (size_t i = 0; i < 3; i++) {
	changed |= lights[i]->on.getValue() != desired[i];
	lights[i]->on = desired[i];
    }
    if (changed)
	this->requestRender("lighting");
}

SbBool
BObolViewController::isLightingEnabled(void) const
{
    SoLightModel *model = controller_light_model(this->d->viewport);
    return model && model->model.getValue() == SoLightModel::PHONG;
}

void
BObolViewController::setLightingProfile(LightingProfile profile)
{
    if (profile != LIGHTING_STUDIO && profile != LIGHTING_MGED)
	return;
    if (this->d->lightingProfile == profile)
	return;

    SoDirectionalLight *lights[3] = {
	controller_headlight(this->d->viewport),
	controller_studio_fill(this->d->viewport),
	controller_studio_rim(this->d->viewport)
    };
    SoEnvironment *environment = controller_environment(this->d->viewport);
    if (!lights[0] || !lights[1] || !lights[2] || !environment)
	return;

    this->d->lightingProfile = profile;
    const SbBool studio = profile == LIGHTING_STUDIO;
    this->d->headlightOffsetEye = studio ?
	bobol_headlight_default_offset() : bobol_mged_headlight_offset();
    this->d->headlightOffsetEye.normalize();

    environment->ambientColor = SbColor(1.0f, 1.0f, 1.0f);
    environment->ambientIntensity = studio ? 0.18f : 0.30f;
    for (size_t i = 0; i < 3; i++)
	lights[i]->color = SbColor(1.0f, 1.0f, 1.0f);
    lights[0]->intensity = studio ? 0.68f : 1.0f;
    lights[1]->intensity = 0.22f;
    lights[2]->intensity = 0.18f;

    const SbBool lightOn = this->d->headlightEnabled &&
	this->isLightingEnabled();
    lights[0]->on = lightOn;
    lights[1]->on = lightOn && studio;
    lights[2]->on = lightOn && studio;
    this->applyTrackedHeadlight(TRUE);
    this->requestRender("lighting-profile");
}

BObolViewController::LightingProfile
BObolViewController::getLightingProfile(void) const
{
    return this->d->lightingProfile;
}

float
BObolViewController::getLightingAmbientIntensity(void) const
{
    SoEnvironment *environment = controller_environment(this->d->viewport);
    return environment ? environment->ambientIntensity.getValue() : 0.0f;
}

void
BObolViewController::setNormalStyle(BObolViewLodState::NormalStyle style,
	float creaseAngleDegrees)
{
    BObolViewLodState *viewState = this->getViewLodState();
    if (!viewState)
	return;
    const BObolViewLodState::NormalStyle beforeStyle =
	viewState->getNormalStyle();
    const float beforeAngle = viewState->getNormalCreaseAngle();
    viewState->setNormalStyle(style, creaseAngleDegrees);
    if (viewState->getNormalStyle() != beforeStyle ||
	std::fabs(viewState->getNormalCreaseAngle() - beforeAngle) > 1.0e-6f)
	this->requestRender("normal-style");
}

BObolViewLodState::NormalStyle
BObolViewController::getNormalStyle(void) const
{
    const BObolViewLodState *viewState = this->getViewLodState();
    return viewState ? viewState->getNormalStyle() :
	BObolViewLodState::NORMAL_AUTHORED;
}

float
BObolViewController::getNormalCreaseAngle(void) const
{
    const BObolViewLodState *viewState = this->getViewLodState();
    return viewState ? viewState->getNormalCreaseAngle() : 60.0f;
}

void
BObolViewController::setHeadlightEnabled(SbBool enabled)
{
    if (this->d->headlightEnabled == enabled)
	return;
    this->d->headlightEnabled = enabled;
    SoDirectionalLight *lights[3] = {
	controller_headlight(this->d->viewport),
	controller_studio_fill(this->d->viewport),
	controller_studio_rim(this->d->viewport)
    };
    if (lights[0] && lights[1] && lights[2]) {
	/* Only illuminate when master lighting (PHONG) is also on. */
	const SbBool lightOn = enabled && this->isLightingEnabled();
	lights[0]->on = lightOn;
	const SbBool studio = this->d->lightingProfile == LIGHTING_STUDIO;
	lights[1]->on = lightOn && studio;
	lights[2]->on = lightOn && studio;
	if (lightOn)
	    this->applyTrackedHeadlight(TRUE);
    }
    this->requestRender("lighting");
}

SbBool
BObolViewController::isHeadlightEnabled(void) const
{
    return this->d->headlightEnabled;
}

void
BObolViewController::setHeadlightCameraTracked(SbBool tracked)
{
    if (this->d->headlightCameraTracked == tracked)
	return;
    this->d->headlightCameraTracked = tracked;
    /* Re-enabling: recompute direction now from the last camera orientation so
     * the change is visible without waiting for the next camera motion.
     * Disabling: leave the current direction in place (scene-fixed). */
    if (tracked)
	this->applyTrackedHeadlight();
    this->requestRender("lighting");
}

SbBool
BObolViewController::isHeadlightCameraTracked(void) const
{
    return this->d->headlightCameraTracked;
}

void
BObolViewController::setHeadlightOffset(const SbVec3f &eyeDir)
{
    SbVec3f dir = eyeDir;
    if (dir.normalize() <= 0.0f)
	return; /* ignore a degenerate direction */
    if (this->d->headlightOffsetEye == dir)
	return;
    this->d->headlightOffsetEye = dir;
    this->applyTrackedHeadlight(TRUE);
    this->requestRender("lighting");
}

SbVec3f
BObolViewController::getHeadlightOffset(void) const
{
    return this->d->headlightOffsetEye;
}

SbVec3f
BObolViewController::getHeadlightDirection(void) const
{
    /* Current world-space travel direction of the headlight node (updated each
     * view sync by applyTrackedHeadlight when camera-tracked). */
    SoDirectionalLight *light = controller_headlight(this->d->viewport);
    return light ? light->direction.getValue() : SbVec3f(0.0f, 0.0f, -1.0f);
}

void
BObolViewController::getCameraLights(
	std::vector<BObolSceneLightRealization> &lights) const
{
    lights.clear();
    if (!this->isLightingEnabled() || !this->d->headlightEnabled)
	return;
    SoDirectionalLight *nodes[3] = {
	controller_headlight(this->d->viewport),
	controller_studio_fill(this->d->viewport),
	controller_studio_rim(this->d->viewport)
    };
    const char *names[3] = {"camera-key", "camera-fill", "camera-rim"};
    for (size_t i = 0; i < 3; i++) {
	if (!nodes[i] || !nodes[i]->on.getValue() ||
	    nodes[i]->intensity.getValue() <= 0.0f)
	    continue;
	BObolSceneLightRealization light;
	light.kind = BOBOL_SCENE_LIGHT_DIRECTIONAL;
	light.name = names[i];
	light.direction = nodes[i]->direction.getValue();
	if (light.direction.length() > 0.0f)
	    light.direction.normalize();
	light.color = nodes[i]->color.getValue();
	light.intensity = nodes[i]->intensity.getValue();
	lights.push_back(light);
    }
}

void
BObolViewController::setSceneLightsEnabled(SbBool enabled)
{
    this->d->sceneLightsEnabled = enabled;
    /* Apply to the realized in-scene light group (may be empty until the next
     * realize populates it; rebuildSceneLights() reapplies this flag then). */
    SoGroup *group = controller_scene_lights_group(this->d->viewport);
    if (group) {
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SoNode *child = group->getChild(i);
	    if (child && child->isOfType(SoLight::getClassTypeId()))
		static_cast<SoLight *>(child)->on = enabled;
	}
    }
    this->requestRender("lighting");
}

SbBool
BObolViewController::isSceneLightsEnabled(void) const
{
    return this->d->sceneLightsEnabled;
}

SoNode *
BObolViewController::getSceneLightsRoot(void) const
{
    if (!this->d->viewport || !this->d->viewport->getRoot())
	return NULL;
    SoSeparator *root = this->d->viewport->getRoot();
    for (int i = 0; i < root->getNumChildren(); i++) {
	SoNode *child = root->getChild(i);
	if (child && child->isOfType(SoGroup::getClassTypeId()) &&
	    bu_strcmp(child->getName().getString(),
		controller_scene_lights_group_name()) == 0)
	    return child;
    }
    return NULL;
}

/* Rebuild the in-scene light group from every database source's realized light
 * snapshots.  Called after realization; positions/directions are already
 * world-space so nodes are added directly. */
void
BObolViewController::setSceneLights(
	const std::vector<BObolSceneLightRealization> &lights)
{
    this->d->sceneLights = lights;
    this->rebuildSceneLights();
}

void
BObolViewController::rebuildSceneLights(void)
{
    SoGroup *group = controller_scene_lights_group(this->d->viewport);
    if (!group)
	return;
    group->removeAllChildren();

    for (size_t li = 0; li < this->d->sceneLights.size(); li++) {
	const BObolSceneLightRealization &L = this->d->sceneLights[li];
	SoLight *node = NULL;
	if (L.kind == BOBOL_SCENE_LIGHT_DIRECTIONAL) {
	    SoDirectionalLight *dl = new SoDirectionalLight;
	    dl->direction = L.direction;
	    node = dl;
	} else if (L.kind == BOBOL_SCENE_LIGHT_SPOT) {
	    SoSpotLight *sl = new SoSpotLight;
	    sl->location = L.position;
	    sl->direction = L.direction;
	    /* DB angle is full beam dispersion; Coin cutOffAngle is the
	     * half-angle from the axis.  Clamp to <= 90 degrees. */
	    float cutoff = static_cast<float>(L.coneAngleDeg * 0.5 * DEG2RAD);
	    const float maxCutoff = static_cast<float>(M_PI_2);
	    if (cutoff > maxCutoff)
		cutoff = maxCutoff;
	    sl->cutOffAngle = cutoff;
	    node = sl;
	} else {
	    SoPointLight *pl = new SoPointLight;
	    pl->location = L.position;
	    node = pl;
	}
	node->color = L.color;
	node->intensity = L.intensity;
	node->on = this->d->sceneLightsEnabled;
	group->addChild(node);
    }
}

void
BObolViewController::setHeadlightColor(const SbColor &color)
{
    SoDirectionalLight *light = controller_headlight(this->d->viewport);
    if (!light)
	return;
    const SbColor clamped(
	std::max(0.0f, std::min(1.0f, color[0])),
	std::max(0.0f, std::min(1.0f, color[1])),
	std::max(0.0f, std::min(1.0f, color[2])));
    if (light->color.getValue() == clamped)
	return;
    light->color = clamped;
    this->requestRender("lighting");
}

SbColor
BObolViewController::getHeadlightColor(void) const
{
    SoDirectionalLight *light = controller_headlight(this->d->viewport);
    return light ? light->color.getValue() : SbColor(1.0f, 1.0f, 1.0f);
}

void
BObolViewController::setHeadlightIntensity(float intensity)
{
    SoDirectionalLight *light = controller_headlight(this->d->viewport);
    if (!light || !std::isfinite(intensity))
	return;
    const float clamped = std::max(0.0f, std::min(1.0f, intensity));
    if (std::fabs(light->intensity.getValue() - clamped) <= 1.0e-6f)
	return;
    light->intensity = clamped;
    this->requestRender("lighting");
}

float
BObolViewController::getHeadlightIntensity(void) const
{
    SoDirectionalLight *light = controller_headlight(this->d->viewport);
    return light ? light->intensity.getValue() : 1.0f;
}

void
BObolViewController::setDepthCueEnabled(SbBool enabled)
{
    SoEnvironment *environment = controller_environment(this->d->viewport);
    if (!environment)
	return;
    const int requested = enabled ? SoEnvironment::HAZE :
	SoEnvironment::NONE;
    if (environment->fogType.getValue() == requested)
	return;
    environment->fogType = requested;
    environment->fogColor = this->d->backgroundTop;
    /* Zero delegates visibility distance to the active camera volume. */
    environment->fogVisibility = 0.0f;
    this->requestRender("depth-cue");
}

SbBool
BObolViewController::isDepthCueEnabled(void) const
{
    SoEnvironment *environment = controller_environment(this->d->viewport);
    return environment && environment->fogType.getValue() !=
	SoEnvironment::NONE;
}

void
BObolViewController::setSoftwareWireMode(SoftwareWireMode mode)
{
    if (mode < SOFTWARE_WIRE_AUTO || mode > SOFTWARE_WIRE_FAST)
	mode = SOFTWARE_WIRE_AUTO;
    if (this->d->softwareWireMode == mode)
	return;
    this->d->softwareWireMode = mode;
    if (this->d->renderLodRoot)
	this->d->renderLodRoot->setSoftwareWireMode(mode);
    SoBRLCadRenderBatch *batch =
	dynamic_cast<SoBRLCadRenderBatch *>(this->d->renderBatchRoot);
    if (batch)
	batch->setSoftwareWireMode(mode);
    this->requestRender("software-wire-mode");
}

BObolViewController::SoftwareWireMode
BObolViewController::getSoftwareWireMode(void) const
{
    return this->d->softwareWireMode;
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
BObolViewController::requestRender(const char *reason)
{
    SbBool wakeEndpoint = FALSE;
    {
	std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
	wakeEndpoint = this->d->renderRequested ? FALSE : TRUE;
	this->d->renderRequested = TRUE;
	this->d->renderLodCapacityRelevant = TRUE;
	this->d->renderReason = reason ? reason : "";
	this->d->renderRequestSerial++;
    }
    /*
     * A render request is itself the complete host contract.  Requiring every
     * caller to remember a separate notifyFrameRequest() made camera sync and
     * final calibration replays easy to strand after the progressive timer
     * became idle.  Wake only on the empty-to-pending edge; repeated reasons
     * still update diagnostics without flooding the Qt event queue.
     */
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
	this->d->renderRequestSerial++;
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
    if (callback && (this->hasProgressiveWorkPending() ||
	this->isRenderRequested()))
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

void
BObolViewController::clearRenderRequest(void)
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    this->d->renderRequested = FALSE;
    this->d->renderLodCapacityRelevant = FALSE;
    this->d->renderReason = "";
    this->d->renderRequestSerial++;
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
    this->d->renderRequestSerial++;
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
    this->d->renderRequestSerial++;
}

uint64_t
BObolViewController::renderRequestSerialGet(void) const
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    return this->d->renderRequestSerial;
}

struct ControllerRenderDeadlineContext {
    BObolViewController *controller = NULL;
    uint64_t deadlineNanoseconds = 0;
    SoGLRenderAction::SoGLRenderAbortCB *previous = NULL;
    void *previousData = NULL;
};

static SoGLRenderAction::AbortCode
controller_render_deadline_cb(void *userData)
{
    ControllerRenderDeadlineContext *context =
	static_cast<ControllerRenderDeadlineContext *>(userData);
    if (!context)
	return SoGLRenderAction::CONTINUE;
    if (context->previous) {
	const SoGLRenderAction::AbortCode prior =
	    (*context->previous)(context->previousData);
	if (prior != SoGLRenderAction::CONTINUE)
	    return prior;
    }
    if (!context->controller || !context->deadlineNanoseconds)
	return SoGLRenderAction::CONTINUE;
    return context->controller->beginRenderTiming() >=
	context->deadlineNanoseconds ?
	SoGLRenderAction::ABORT : SoGLRenderAction::CONTINUE;
}

SbBool
BObolViewController::renderPending(SbBool clearWindow,
				     SbBool clearZBuffer,
				     SbString *reason)
{
    /* LoD off -> classic behavior: realize the whole scene before presenting.
     * LoD on -> stay on the progressive coarse-first path (no whole-tree
     * realize here; geometry streams in via advanceProgressiveWork). */
    if (this->isForceRealizeDisplay())
	(void)this->realizePending();
    (void)this->advanceProgressiveWork(NULL, NULL);
	this->synchronizePresentation();


    if (!this->d->renderManager || !this->getRenderContextManager() ||
	!this->d->activeCamera || !this->getRenderRoot())
	return FALSE;

    SbString renderReasonCopy;
    SbBool lodCapacityRelevant = TRUE;
    if (!this->consumeRenderRequest(&renderReasonCopy,
	    &lodCapacityRelevant))
	return FALSE;

    if (reason)
	*reason = renderReasonCopy;

    const uint64_t started = this->beginRenderTiming();
    SoGLRenderAction *renderAction =
	this->d->renderManager->getGLRenderAction();
    ControllerRenderDeadlineContext deadlineContext;
    const uint64_t deadlineDuration = lodCapacityRelevant ?
	this->getCurrentPresentationFrameDeadline() : 0;
    if (renderAction && deadlineDuration) {
	deadlineContext.controller = this;
	deadlineContext.deadlineNanoseconds =
	    started > UINT64_MAX - deadlineDuration ? UINT64_MAX :
	    started + deadlineDuration;
	renderAction->getAbortCallback(
	    deadlineContext.previous, deadlineContext.previousData);
	renderAction->setAbortCallback(
	    controller_render_deadline_cb, &deadlineContext);
    }
    if (clearWindow)
	this->renderBackground();
    else if (clearZBuffer)
	glClear(GL_DEPTH_BUFFER_BIT);
    const uint64_t backgroundCompleted = this->beginRenderTiming();
    const BObolViewLodState *presentationState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    const uint64_t cadExecutionBefore = presentationState ?
	presentationState->cadPresentationExecutionSerial() : 0;
    const uint64_t cadPreparationBefore = presentationState ?
	presentationState->cadPresentationPreparationSerial() : 0;
    if (presentationState)
	presentationState->beginCadPresentationFrame();
    this->d->renderManager->render(static_cast<SbBool>(FALSE), static_cast<SbBool>(FALSE));
    const uint64_t sceneCompleted = this->beginRenderTiming();
    const uint64_t cadExecutionAfter = presentationState ?
	presentationState->cadPresentationExecutionSerial() : 0;
    const uint64_t cadPreparationAfter = presentationState ?
	presentationState->cadPresentationPreparationSerial() : 0;
    if (presentationState)
	presentationState->refreshCadPresentationFrameStatus();
    const SbBool cadFrameIncomplete = presentationState &&
	!presentationState->lastCadPresentationFrameExact();
    const SbBool interrupted =
	(renderAction && renderAction->hasTerminated()) || cadFrameIncomplete;
    if (renderAction && deadlineDuration)
	renderAction->setAbortCallback(
	    deadlineContext.previous, deadlineContext.previousData);
    this->d->lastBackgroundRenderTimeNanoseconds =
	backgroundCompleted >= started ?
	    backgroundCompleted - started : 0;
    this->d->lastSceneRenderTimeNanoseconds =
	sceneCompleted >= backgroundCompleted ?
	    sceneCompleted - backgroundCompleted : 0;
    if (interrupted) {
	this->notePresentationRenderInterrupted(
	    sceneCompleted > started ? sceneCompleted - started : 1,
	    cadExecutionAfter != cadExecutionBefore ? TRUE : FALSE,
	    cadPreparationAfter != cadPreparationBefore ? TRUE : FALSE,
	    lodCapacityRelevant);
	return FALSE;
    }
    this->completeRenderTiming(started, lodCapacityRelevant);
    return TRUE;
}

uint64_t
BObolViewController::beginRenderTiming(void) const
{
    return static_cast<uint64_t>(
	std::chrono::duration_cast<std::chrono::nanoseconds>(
	    std::chrono::steady_clock::now().time_since_epoch()).count());
}

void
BObolViewController::setPresentationFrameDeadlines(
    uint64_t interactiveNanoseconds, uint64_t stableNanoseconds)
{
    this->d->interactivePresentationFrameDeadlineNanoseconds =
	interactiveNanoseconds;
    this->d->stablePresentationFrameDeadlineNanoseconds =
	stableNanoseconds;
    this->d->consecutiveInterruptedPresentationFrames = 0;
}

uint64_t
BObolViewController::getInteractivePresentationFrameDeadline(void) const
{
    return this->d->interactivePresentationFrameDeadlineNanoseconds;
}

uint64_t
BObolViewController::getStablePresentationFrameDeadline(void) const
{
    return this->d->stablePresentationFrameDeadlineNanoseconds;
}

uint64_t
BObolViewController::getCurrentPresentationFrameDeadline(void) const
{
    if (!this->d->lodInteractive && !this->d->lodGestureActive)
	return this->d->stablePresentationFrameDeadlineNanoseconds;

    uint64_t base =
	this->d->interactivePresentationFrameDeadlineNanoseconds;
    /* An ordinary motion frame targets the configured interactive cadence.
     * A deliberate zoom-quality probe has a separate 10 Hz hard floor;
     * aborting it at the ordinary 40 ms software deadline prevents the frame
     * which would calibrate and retain a useful richer cut from completing. */
    if (this->d->lodViewDemandPolicy.qualityBudgetActive())
	base = std::max<uint64_t>(base,
	    BObolLodViewDemandPolicy::qualityFrameDurationNanoseconds());
    if (!base || !this->d->consecutiveInterruptedPresentationFrames)
	return base;

    /* A hard deadline without forward-progress backoff can starve the first
     * useful software frame: every richer coherent retry may cost 45-60 ms,
     * so a 40 ms limit keeps presenting the old background forever.  Permit
     * bounded 50% steps after consecutive aborts, capped by the quiet-frame
     * deadline.  The first completed frame resets this allowance. */
    const uint64_t steps = std::min<uint64_t>(
	this->d->consecutiveInterruptedPresentationFrames, 4u);
    const uint64_t increment = base / 2u;
    const uint64_t candidate =
	increment && steps <= (UINT64_MAX - base) / increment ?
	base + steps * increment : UINT64_MAX;
    const uint64_t stableCap =
	this->d->stablePresentationFrameDeadlineNanoseconds;
    return stableCap ? std::min(candidate, std::max(base, stableCap)) :
	candidate;
}

void
BObolViewController::notePresentationRenderInterrupted(
    uint64_t elapsedNanoseconds, SbBool cadDrawAttempted,
    SbBool cadPreparationChanged, SbBool lodCapacityRelevant)
{
    if (!elapsedNanoseconds)
	return;
    this->d->interruptedPresentationFrameCount++;
    this->d->lastInterruptedPresentationTimeNanoseconds =
	elapsedNanoseconds;
    if (!lodCapacityRelevant) {
	/* A presentation patch may still encounter SoCADAssembly's resumable
	 * command preparation.  Finish that patch, but do not let it cancel a
	 * geometry headroom witness, lower a PoP ceiling, or alter a scene budget. */
	this->requestPresentationRender("render-presentation-replay");
	return;
    }
    this->d->consecutiveInterruptedPresentationFrames++;
    /* A pending headroom witness belongs to the cut which just failed its
     * hard presentation contract.  It may not reopen admission after the
     * constrained retry completes. */
    this->d->lodHeadroomPolicy.cancelRetry();

    const BObolViewLodState *state = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    const SbBool interactive =
	this->d->lodInteractive || this->d->lodGestureActive;
    const int interruptedActiveMaximum = state ?
	state->maximumActiveProgressiveCut() : -1;
    const int interruptedPresentedMaximum =
	interruptedActiveMaximum >= 0 &&
	this->d->lodInteractiveProgressiveCeiling >= 0 ?
	    std::min(interruptedActiveMaximum,
		this->d->lodInteractiveProgressiveCeiling) :
	    interruptedActiveMaximum;
    /* Any abort which advanced retained preparation is forward progress, not
	 * draw-capacity evidence.  Large command records are deliberately built in
	 * resumable deadline slices and may need more than one slice at 150k+
	 * occurrences.  Coarsening the second slice invalidates the work just
	 * prepared and creates a normalize/refine loop.  A genuinely stuck retry
	 * changes no preparation token and reaches the ordinary pressure path. */
    const SbBool preparationOnlyRetry =
	cadPreparationChanged ? TRUE : FALSE;
    if (getenv("BOBOL_LOD_TRACE_DEADLINE"))
	bu_log("BObol LoD render deadline elapsed_ms=%.3f "
	       "draw=%d preparation=%d interactive=%d "
	       "active_max=%d presented_max=%d ceiling=%d "
	       "active_cost=%zu minimum_cost=%zu handoff=%d\n",
	       elapsedNanoseconds / 1000000.0,
	       cadDrawAttempted ? 1 : 0,
	       cadPreparationChanged ? 1 : 0,
	       interactive ? 1 : 0,
	       interruptedActiveMaximum,
	       interruptedPresentedMaximum,
	       this->d->lodInteractiveProgressiveCeiling,
	       state ? state->activeRenderCost() : 0,
	       state ? state->minimumActiveRenderCost() : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0);
    if (!interactive && cadDrawAttempted && !preparationOnlyRetry &&
	this->d->lodInteractiveProgressiveCeiling < 0 &&
	this->d->lodBudgetPolicy.retainedQualityFloorActive()) {
	const size_t floorBudget =
	    this->d->lodBudgetPolicy.retainedQualityFloorBudget();
	const uint64_t floorSignature =
	    this->d->lodBudgetPolicy.retainedQualityFloorSignature();
	const unsigned int missesBefore =
	    this->d->lodBudgetPolicy.retainedQualityFloorMissCount();
	const bool rejected =
	    this->d->lodBudgetPolicy.noteRetainedQualityFloorMiss();
	if (getenv("BOBOL_LOD_TRACE_BUDGET"))
	    bu_log("BObol LoD quality-floor miss budget=%zu "
		   "signature=%llu miss=%u elapsed_ms=%.3f "
		   "active_cost=%zu active_max=%d rejected=%d\n",
		   floorBudget, (unsigned long long)floorSignature,
		   missesBefore + 1, elapsedNanoseconds / 1000000.0,
		   state ? state->activeRenderCost() : 0,
		   interruptedActiveMaximum, rejected ? 1 : 0);
    }

    /* A Coin abort bounds total endpoint-thread traversal, not a completed
     * geometry presentation.  In particular, a newly rotated many-part CAD
     * scene may spend the whole deadline rebuilding command records before
     * the GPU cut can be measured.  activeRenderCost() is retained occurrence
     * demand (often richer than the renderer ceiling), so pairing it with
     * this elapsed time produced a fictitious low triangle rate and
     * destructively compacted otherwise proven resident PoP prefixes.
     *
     * Keep deadline recovery presentation-only below: lower the reversible
     * renderer ceiling and/or aggregate sub-pixel occurrences.  Only a
     * completed frame publishes exact presented work and may update the
     * persistent scene budget or throughput calibration. */

    /* A deadline abort is already sufficient evidence that the submitted
     * progressive cut is too expensive; waiting for a completed frame to
     * correct it creates a retry deadlock.  This applies both to a scale
     * interaction and to a quiet admission whose first coherent frame
     * exceeded the stable deadline.  Lower only the renderer-wide prefix
     * ceiling by one PoP bit.  The occurrence cut and resident suffix remain
     * intact, so this is immediate, reversible, and cannot restart level
     * walking.  A later bounded stable pass reconciles the retained cut before
     * removing the ceiling.
     *
     * This applies to pose-only motion as well as zoom.  A prepared pose frame
     * normally lets completed-frame feedback choose the ceiling, but a hard
     * abort supplies no such sample; leaving pose motion excluded retries the
     * same over-budget cut until button-up.  Lowering only this global ceiling
     * keeps every occurrence cut and resident suffix available for immediate
     * quiet restoration. */
    if (state && cadDrawAttempted && !preparationOnlyRetry) {
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	if (interruptedPresentedMaximum > 0) {
	    const int correctedCeiling = interruptedPresentedMaximum - 1;
	    this->d->lodViewDemandPolicy.noteQualityMiss(
		correctedCeiling, state->activeRenderCost());
	    this->d->lodInteractiveProgressiveCeiling = correctedCeiling;
	    presentationState->setCadPresentationProgressiveCutCeiling(
		correctedCeiling);
	    if (!interactive) {
		this->d->lodPresentationPolicy.armHandoff(true);
		/* The camera and policy epochs did not change, but the aborted cut
		 * still requires an owner-thread admission pass to reconcile the
		 * retained occurrences before the temporary render ceiling can be
		 * removed.  Without invalidating this proof, the ordinary scene
		 * submission fast path declares "already submitted" and leaves the
		 * handoff latch pending forever. */
		this->d->lodLastSubmittedViewRevision.reset();
		this->d->lodLastSubmittedPolicyRevision.reset();
	    }
	}
    }

    /* A many-part view can hit a harder floor than PoP triangle detail: one
     * coherent minimum prefix per visible occurrence.  The normal small-part
     * aggregation controller learns only from completed frames, which is an
     * impossible prerequisite when every attempted minimum-prefix frame is
     * itself interrupted.  Treat a quiet deadline abort as conservative
     * feedback for the reversible presentation-only point cut as well.  This
     * keeps all immutable meshes and desired cuts resident, aggregates only
     * occurrences below the camera-local screen threshold, and lets the
     * ordinary completed-frame relaxation probe recover the finest sustainable
     * threshold afterward.
     *
     * Do not apply this to active input: pose/scale interaction already owns
     * a separately measured cut and repeated unapplied motion feedback could
     * otherwise over-coarsen the first post-input frame. */
    /* Unlike the renderer-wide PoP ceiling above, a point threshold changes
     * the occurrence population itself.  No preparation-heavy abort is valid
	 * evidence for that bracket: a lower PoP ceiling may make retained
	 * construction finish, while poisoning the
     * point bracket would preserve visibly coarse multi-pixel objects after
     * the one-time work is gone. */
    const size_t interruptedActiveCost = state ?
	state->activeRenderCost() : 0;
    const size_t interruptedMinimumCost = state ?
	state->minimumActiveRenderCost() : interruptedActiveCost;
    /* activeRenderCost() describes the retained occurrence cuts, while the
     * renderer-wide ceiling describes what this interrupted frame actually
     * attempted.  Once that reversible ceiling reaches the minimum PoP cut,
     * triangle population presented by every progressive occurrence is
     * already irreducible even if rich resident suffixes keep active cost
     * above minimumActiveRenderCost().  Requiring retained cost itself to be
     * minimal deadlocks a many-part software view at the minimum cut: no further
     * triangle correction is possible, yet the small-part aggregation
     * controller is never allowed to act. */
    const SbBool interruptedPopulationIrreducible =
	(interruptedActiveCost <= interruptedMinimumCost ||
	 interruptedPresentedMaximum == 0) ? TRUE : FALSE;
    if (state && state->hasCadPresentationAssemblies() &&
	cadDrawAttempted && !cadPreparationChanged && !interactive &&
	interruptedPopulationIrreducible &&
	this->d->quietAllocationTargetFps() > 0.0f) {
	const BObolLodPointProxyCalibrationPolicy::Decision pressure =
	    this->d->lodPointProxyCalibrationPolicy.interrupted(
		this->d->lodPresentationPointProxyPixelThreshold,
		elapsedNanoseconds, this->d->quietAllocationTargetFps());
	if (pressure.changed) {
	    this->d->lodPresentationPointProxyPixelThreshold =
		pressure.threshold;
	    BObolViewLodState *presentationState =
		this->d->viewAttachment->getViewLodState();
	    presentationState->setCadPresentationPointProxyPixelThreshold(
		pressure.threshold);
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	}
    }

    if (preparationOnlyRetry) {
	/* Preserve every admission and calibration barrier until retained
	 * preparation completes and an exact unchanged replay can measure the
	 * selected cut. */
	this->requestRender("render-preparation-replay");
	return;
    }

    /* An interrupted traversal requires another bounded policy/presentation
     * attempt, but it is not a persistent capacity sample.  Keep the barriers
     * armed while the renderer-only corrections above make that retry
     * cheaper; retained-cut admission remains based on completed frames. */
    this->d->lodBudgetPolicy.resetOverloadRecovery();
    this->d->lodBudgetPolicy.setProbeCount(3);
    this->d->lodBudgetPolicy.resetPass();
    this->markProgressiveWorkPending();
    this->requestRender("render-deadline");
}

uint64_t
BObolViewController::getInterruptedPresentationFrameCount(void) const
{
    return this->d->interruptedPresentationFrameCount;
}

uint64_t
BObolViewController::getLastInterruptedPresentationTimeNanoseconds(void) const
{
    return this->d->lastInterruptedPresentationTimeNanoseconds;
}

void
BObolViewController::armStableLodHeadroomProbeIfReady(void)
{
    const BObolViewLodState *viewLodState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    const size_t activePopulationCost = viewLodState ?
	viewLodState->activeRenderCost() : 0;
    const size_t activeTasks = this->d->lodService ?
	this->d->lodService->activeTaskCountForGeneration(
	    this->d->lodActiveGeneration) : 0;
    const size_t queuedResults = this->d->lodService ?
	this->d->lodService->queuedResultCountForGeneration(
	    this->d->lodActiveGeneration) : 0;
    size_t activePayloads = 0;
    size_t satisfiedPayloads = 0;
    size_t memoryLimitedPayloads = 0;
    if (viewLodState)
	viewLodState->convergencePayloadCounts(activePayloads,
	    satisfiedPayloads, memoryLimitedPayloads);

    size_t presentedSubpixelOccurrences = 0;
    size_t presentedStructuralBoxes = 0;
    const bool exactOccurrenceCoverage = viewLodState &&
	viewLodState->lastCadPresentationOccurrenceCoverage(
	    presentedSubpixelOccurrences, presentedStructuralBoxes);
    const size_t terminalOccurrenceFailures = viewLodState ?
	viewLodState->cadOccurrenceTerminalFailureCount() : 0;

    /* Initial structural coverage and terminal presentation coverage are
     * deliberately different proofs.  The former gets a complete model on
     * screen quickly; the latter permits only meshes, camera-valid subpixel
     * points, or an explicit provider failure.  A bounded source scan can
     * finish its box-first pass and later exhaust a refinement window without
     * visiting every remaining box.  If an exact quiet frame observes that
     * state, resume a normal (not structural-only) current-view pass from the
     * retained population.  No source data or PoP prefix is discarded.
     *
     * This is also the resize recovery edge: a resize may change which boxes
     * are physically subpixel while preserving the view/source epochs.  The
     * completed framebuffer is the authoritative classification, so it must
     * be able to re-arm admission without waiting for another mouse event. */
    const bool unresolvedStructuralPresentation =
	exactOccurrenceCoverage &&
	presentedStructuralBoxes > terminalOccurrenceFailures;
    const bool structuralRepairReady = unresolvedStructuralPresentation &&
	this->d->lodAutoSubmit && this->d->lodService &&
	!this->d->lodInteractive && !this->d->lodGestureActive &&
	this->d->lodCoveragePolicy.coverageComplete() &&
	!this->d->lodSubmissionPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodRefinementAwaitingFrame &&
	!this->d->lodRetainedRefinementPending &&
	!this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodStablePointProxyCalibrationPending &&
	!this->d->lodPointProxyTriangleRecoveryPending &&
	this->d->lodResultsPending.load() == 0 && activeTasks == 0 &&
	queuedResults == 0;
    if (structuralRepairReady) {
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionDeltaActive = FALSE;
	this->d->lodSubmissionDeltaSources.clear();
	this->d->lodSubmissionDeltaPlans.clear();
	this->d->lodSubmissionRefreshMissing = TRUE;
	this->d->lodSubmissionReset = 0;
	this->d->lodSubmissionPending = TRUE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodBudgetPolicy.clearBudgetLimit();
	this->d->lodBudgetPolicy.resetPass();
	this->markProgressiveWorkPending();
	this->requestRender("lod-structural-presentation-repair");
	this->notifyFrameRequest("lod-structural-presentation-repair");
	return;
    }

    /* A complete one-pixel image is the first stable result.  When that exact
     * retained population uses less than one quarter of both the measured
     * scene allowance and resident-memory allowance, continue once to the
     * next subpixel tier.  This is intentionally distinct from the
     * budget-debt retry below: every current request is satisfied, so there is
     * no blocked occurrence capable of arming that mechanism.
     *
     * Honor an explicit LoD scale setting as the user's quality decision.
     * The automatic tier is only the default-scale policy; a scale above one
     * already requests subpixel data, while a scale below one deliberately
     * trades fidelity for capacity. */
    struct bv_lod_policy viewPolicy;
    bv_lod_policy_init(&viewPolicy);
    this->d->viewAttachment->getLodPolicy(&viewPolicy);
    const bool defaultQualityScale = std::isfinite(viewPolicy.scale) &&
	std::fabs(viewPolicy.scale - 1.0) <= 1.0e-6;
    size_t presentedRenderCost = 0;
    const bool exactCompletedFrame = viewLodState &&
	viewLodState->lastCadPresentationFrameExact() &&
	viewLodState->lastCadPresentedRenderCost(presentedRenderCost) &&
	presentedRenderCost > 0;
    bool residentMemoryHeadroom = false;
    if (this->d->lodService) {
	const size_t residentLimit =
	    this->d->lodService->getResidentMeshLimit();
	const size_t residentBytes = this->d->lodService->
	    stableResidentMeshBytesForDiagnostics();
	const size_t reservedGrowth = this->d->lodService->
	    reservedResidentMeshGrowthBytesForDiagnostics();
	residentMemoryHeadroom = residentLimit == SIZE_MAX ||
	    (residentBytes <= residentLimit / 4 &&
	     reservedGrowth <= residentLimit / 4 -
		std::min(residentBytes, residentLimit / 4));
    }
    const bool stableTerminalContext = this->d->lodAutoSubmit &&
	this->d->lodService && !this->d->lodInteractive &&
	!this->d->lodGestureActive && defaultQualityScale &&
	!this->d->lodUseForcedCut &&
	this->d->lodCoveragePolicy.coverageComplete() &&
	activePayloads > 0 && memoryLimitedPayloads == 0 &&
	!this->d->lodSubmissionPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodRefinementAwaitingFrame &&
	!this->d->lodRetainedRefinementPending &&
	!this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodStablePointProxyCalibrationPending &&
	!this->d->lodPointProxyTriangleRecoveryPending &&
	this->d->lodResultsPending.load() == 0 && activeTasks == 0 &&
	queuedResults == 0;
    const bool stableSubpixelContext = stableTerminalContext &&
	activePayloads == satisfiedPayloads &&
	this->d->lodInteractiveProgressiveCeiling < 0;
    if (getenv("BOBOL_LOD_TRACE_HEADROOM") &&
	this->d->lodTargetPixelError > 0.5001f &&
	this->d->lodTargetPixelError <= 1.01f)
	bu_log("BObol LoD subpixel arm eligible=%d exact=%d default=%d "
	       "coverage=%d submit=%d rescan=%d refinement_frame=%d "
	       "retained=%d handoff=%d point=%d recovery=%d ceiling=%d "
	       "results=%d active_tasks=%zu queued_results=%zu "
	       "active_cost=%zu budget=%zu active_payloads=%zu "
	       "satisfied=%zu memory_limited=%zu memory_headroom=%d "
	       "last_ms=%.3f\n",
	       stableSubpixelContext ? 1 : 0, exactCompletedFrame ? 1 : 0,
	       defaultQualityScale ? 1 : 0,
	       this->d->lodCoveragePolicy.coverageComplete() ? 1 : 0,
	       this->d->lodSubmissionPending ? 1 : 0,
	       this->d->lodBudgetPolicy.rescanAfterFrame() ? 1 : 0,
	       this->d->lodRefinementAwaitingFrame ? 1 : 0,
	       this->d->lodRetainedRefinementPending ? 1 : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
	       this->d->lodStablePointProxyCalibrationPending ? 1 : 0,
	       this->d->lodPointProxyTriangleRecoveryPending ? 1 : 0,
	       this->d->lodInteractiveProgressiveCeiling,
	       this->d->lodResultsPending.load(), activeTasks, queuedResults,
	       activePopulationCost, this->d->lodBudgetPolicy.currentBudget(),
	       activePayloads, satisfiedPayloads, memoryLimitedPayloads,
	       residentMemoryHeadroom ? 1 : 0,
	       this->d->lastRenderTimeNanoseconds / 1000000.0);
    const float stablePixelError = stableSubpixelContext ?
	BObolLodQualityPolicy::stablePixelError(
	    this->d->lodTargetPixelError, activePopulationCost,
	    this->d->lodBudgetPolicy.currentBudget(),
	    this->d->lastRenderTimeNanoseconds,
	    this->d->quietAllocationTargetFps(),
	    exactCompletedFrame, residentMemoryHeadroom) :
	this->d->lodTargetPixelError;
    if (stablePixelError + 1.0e-6f < this->d->lodTargetPixelError) {
	this->d->lodTargetPixelError = stablePixelError;
	this->advanceLodPolicyRevision();
	this->d->lodLastSubmittedPolicyRevision.reset();
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPending = TRUE;
	this->d->lodSubmissionRescanPending = FALSE;
	this->markProgressiveWorkPending();
	this->requestRender("lod-stable-subpixel");
	this->notifyFrameRequest("lod-stable-subpixel");
	return;
    }
    const bool actionableQualityDebt =
	activePayloads > satisfiedPayloads &&
	activePayloads - satisfiedPayloads > memoryLimitedPayloads;

    /* The preferred quiet cadence is an initial convergence target, not a
     * permanent fidelity ceiling for an event-driven static image.  Once the
     * ordinary allocation is terminal, exact, fully covered, and free of
     * memory pressure, recompute the scene-wide importance allocation using
     * the separate hard static-frame deadline.  The completed framebuffer is
     * retained without redraw; any later input immediately leaves this mode
     * and restores the motion budget from the existing occurrence cuts.
     *
     * This is one phase transition, not an open-ended headroom probe.  The
     * ordinary deadline recovery bounds a discrete PoP jump which turns out
     * to be too expensive, while pixel demand and resident-memory limits
     * remain authoritative terminal conditions. */
    const uint64_t preferredQuietNanoseconds =
	this->d->lodStableTargetFps > 0.0f ?
	    static_cast<uint64_t>(1000000000.0L /
		static_cast<long double>(this->d->lodStableTargetFps)) : 0;
    const bool staticOverscanEligible = stableTerminalContext &&
	exactCompletedFrame && actionableQualityDebt &&
	this->d->lodBudgetPolicy.stableBudgetLimited() &&
	!this->d->lodStaticOverscanActive &&
	!this->d->lodResourcePolicy.anyPressure() &&
	this->d->stablePresentationFrameDeadlineNanoseconds >
	    preferredQuietNanoseconds;
    if (staticOverscanEligible) {
	this->d->lodStaticOverscanActive = TRUE;
	/* A quiet handoff may have retained the responsive presentation ceiling
	 * learned during zoom.  That ceiling is valuable while input is active,
	 * but leaving it installed here prevents an already resident richer cut
	 * from ever participating in the bounded static-frame trial.  Remove only
	 * the renderer-side ceiling: occurrence demand and immutable resident PoP
	 * storage stay unchanged.  If the richer frame misses the hard static
	 * deadline, notePresentationRenderInterrupted() restores a one-cut-lower
	 * ceiling and the previous completed framebuffer remains available. */
	const int priorCeiling =
	    this->d->lodInteractiveProgressiveCeiling;
	this->d->lodInteractiveProgressiveCeiling = -1;
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	if (presentationState) {
	    presentationState->setCadPresentationProgressiveCutCeiling(-1);
	    presentationState->setCadPresentationCameraMotionFrameReuse(FALSE);
	}
	this->d->lodPresentationPolicy.cancelHandoff();
	this->d->lodHeadroomPolicy.cancelRetry();
	this->d->lodBudgetPolicy.clearBudgetLimit();
	this->d->lodBudgetPolicy.resetProbeSeries();
	this->d->lodBudgetPolicy.resetOverloadRecovery();
	this->d->lodBudgetPolicy.resetPass();
	this->d->lodBudgetPolicy.requestRetainedReallocation(false);
	this->d->lodLastSubmittedViewRevision.reset();
	this->d->lodLastSubmittedPolicyRevision.reset();
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPending = TRUE;
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodDiscretePopulationTrialAvailable = TRUE;
	if (getenv("BOBOL_LOD_TRACE_HEADROOM"))
	    bu_log("BObol LoD static overscan armed prior_ceiling=%d "
		   "active_max=%d active_cost=%zu budget=%zu "
		   "target_fps=%.3f last_ms=%.3f\n",
		   priorCeiling,
		   presentationState ? presentationState->
		       maximumActiveProgressiveCut() : -1,
		   activePopulationCost,
		   this->d->lodBudgetPolicy.currentBudget(),
		   this->d->quietAllocationTargetFps(),
		   this->d->lastRenderTimeNanoseconds / 1000000.0);
	this->markProgressiveWorkPending();
	this->requestRender("lod-static-overscan");
	this->notifyFrameRequest("lod-static-overscan");
	return;
    }
    /* A terminal planning pass can finish while the motion-to-stable or
     * small-part presentation barrier still owns the next frame.  Arm one
     * exact, unchanged replay from either side of that barrier.  The
     * headroom policy makes the (view, policy, active population) witness
     * one-shot, so an actually capacity-limited population cannot repaint in
     * a loop. */
    const bool eligible = this->d->lodAutoSubmit && this->d->lodService &&
	!this->d->lodInteractive && !this->d->lodGestureActive &&
	this->d->lodBudgetPolicy.stableBudgetLimited() &&
	!this->d->lodSubmissionPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodRefinementAwaitingFrame &&
	!this->d->lodRetainedRefinementPending &&
	!this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodStablePointProxyCalibrationPending &&
	!this->d->lodPointProxyTriangleRecoveryPending &&
	this->d->lodInteractiveProgressiveCeiling < 0 &&
	this->d->lodBudgetPolicy.currentBudget() != SIZE_MAX &&
	actionableQualityDebt &&
	this->d->lodResultsPending.load() == 0 && activeTasks == 0 &&
	queuedResults == 0;
    if (getenv("BOBOL_LOD_TRACE_HEADROOM"))
	bu_log("BObol LoD headroom arm eligible=%d interactive=%d gesture=%d "
	       "limited=%d submit=%d rescan=%d refinement_frame=%d "
	       "retained=%d handoff=%d point=%d recovery=%d ceiling=%d "
	       "results=%d active_tasks=%zu queued_results=%zu active_cost=%zu "
	       "budget=%zu quality_debt=%d active_payloads=%zu satisfied=%zu "
	       "memory_limited=%zu pending=%d\n",
	       eligible ? 1 : 0, this->d->lodInteractive ? 1 : 0,
	       this->d->lodGestureActive ? 1 : 0,
	       this->d->lodBudgetPolicy.stableBudgetLimited() ? 1 : 0,
	       this->d->lodSubmissionPending ? 1 : 0,
	       this->d->lodBudgetPolicy.rescanAfterFrame() ? 1 : 0,
	       this->d->lodRefinementAwaitingFrame ? 1 : 0,
	       this->d->lodRetainedRefinementPending ? 1 : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
	       this->d->lodStablePointProxyCalibrationPending ? 1 : 0,
	       this->d->lodPointProxyTriangleRecoveryPending ? 1 : 0,
	       this->d->lodInteractiveProgressiveCeiling,
	       this->d->lodResultsPending.load(), activeTasks, queuedResults,
	       activePopulationCost, this->d->lodBudgetPolicy.currentBudget(),
	       actionableQualityDebt ? 1 : 0, activePayloads,
	       satisfiedPayloads, memoryLimitedPayloads,
	       this->d->lodHeadroomPolicy.retryPending() ? 1 : 0);
    if (!eligible)
	return;

    if (!activePopulationCost ||
	!this->d->lodHeadroomPolicy.armRetry(
	    this->d->lodViewRevision, this->d->lodPolicyRevision,
	    activePopulationCost,
	    this->d->lodBudgetPolicy.currentBudget()))
	return;

    if (getenv("BOBOL_LOD_TRACE_HEADROOM"))
	bu_log("BObol LoD headroom armed view=%llu policy=%llu "
	       "active_cost=%zu budget=%zu\n",
	       (unsigned long long)this->d->lodViewRevision.value(),
	       (unsigned long long)this->d->lodPolicyRevision.value(),
	       activePopulationCost, this->d->lodBudgetPolicy.currentBudget());

    this->markProgressiveWorkPending();
    this->requestRender("lod-calibrated-headroom-probe");
    this->notifyFrameRequest("lod-calibrated-headroom-probe");
}

void
BObolViewController::resumeLodAfterOnePixelRecovery(void)
{
    if (!this->d->lodBudgetPolicy.
	    confirmRetainedRecoveryPresentation(true))
	return;

    /* The recovery ceiling is a one-frame guard, not a fidelity policy.  Its
     * coherent one-pixel successor is now visible, so let measured headroom
     * grow from the retained population in bounded screen-priority waves. */
    this->d->lodLastSubmittedViewRevision.reset();
    this->d->lodLastSubmittedPolicyRevision.reset();
    this->d->lodSubmissionSourceIndex = 0;
    this->d->lodSubmissionEntryOffset = 0;
    this->d->clearLodSubmissionPlan();
    this->d->lodSubmissionPending = TRUE;
    this->d->lodSubmissionRescanPending = FALSE;
    this->markProgressiveWorkPending();
    this->notifyFrameRequest("lod-recovery-headroom");
}

void
BObolViewController::completeRenderTiming(uint64_t startedNanoseconds,
	SbBool lodCapacityRelevant)
{
    const uint64_t now = this->beginRenderTiming();
    if (!startedNanoseconds || now <= startedNanoseconds)
	return;

    const uint64_t elapsed = now - startedNanoseconds;
    if (elapsed >= 30000000000ULL)
	return;

    this->d->lastRenderTimeNanoseconds = elapsed;
    this->d->smoothedRenderTimeNanoseconds =
	this->d->smoothedRenderTimeNanoseconds ?
	(this->d->smoothedRenderTimeNanoseconds * 9 + elapsed) / 10 : elapsed;

    /* Presentation-only frames still contribute user-facing frame-time
     * telemetry and, when exact, are real presentations of the current CAD
     * population.  They therefore retire publication/refinement barriers,
     * but must not satisfy/cancel capacity probes or teach the retained
     * geometry allocator from a one-time style patch. */
    if (!lodCapacityRelevant) {
	this->completePresentationBarrier(elapsed);
	return;
    }

    /* A completed presentation ends any deadline-backoff sequence.  Aborted
     * traversals never enter the persistent capacity estimator: they do not
     * publish an exact presented-work denominator. */
    this->d->consecutiveInterruptedPresentationFrames = 0;

    /* Calibrate geometry throughput from frames that actually traversed the
     * current aggregate cut.  This is deliberately scene-level: per-object
     * projected error says what would look good, but only observed aggregate
     * work says what this renderer, viewport, draw mode, and machine can
     * sustain.
     *
     * Motion and quiet views have separate calibrations.  A settled retained
     * frame can be exceptionally cheap (especially when one mesh is instanced
     * thousands of times), while a motion frame also pays for cut selection
     * and presentation updates.  Mixing those samples taught the interaction
     * policy a budget it could not sustain. */
    const BObolViewLodState *calibrationState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    /* One complete-frame sample supplies triangle, GPU timer, replay, and
     * resource-pressure consumers.  Sampling here avoids four independent
     * source-presentation walks in the controller and HUD paths. */
    if (calibrationState)
	calibrationState->refreshCadPresentationFrameStatus();
    BObolViewLodState::CadGpuResourceStatus frameGpuResources;
    const SbBool haveFrameGpuResources = calibrationState &&
	calibrationState->cadGpuResourceStatus(frameGpuResources);
    const SbBool gpuPressure = haveFrameGpuResources &&
	frameGpuResources.memoryPressure ? TRUE : FALSE;
    SbBool cpuPressure = FALSE;
    if (this->d->lodService) {
	const size_t residentLimit =
	    this->d->lodService->getResidentMeshLimit();
	const size_t stableResident = this->d->lodService->
	    stableResidentMeshBytesForDiagnostics();
	const size_t reservedGrowth = this->d->lodService->
	    reservedResidentMeshGrowthBytesForDiagnostics();
	cpuPressure = residentLimit != SIZE_MAX &&
	    (stableResident > residentLimit ||
	     reservedGrowth > residentLimit -
		std::min(stableResident, residentLimit)) ? TRUE : FALSE;
    }
    const BObolLodResourcePolicy::Decision resourceDecision =
	this->d->lodResourcePolicy.observe(
	    cpuPressure != FALSE, gpuPressure != FALSE,
	    this->d->lodAutoSubmit && this->d->lodService);
    if (resourceDecision.queueRecovery) {
	/* Resource recovery obeys the same coverage and interaction safety gates
	 * as ordinary stable compaction.  The coordinator owns the pressure edge;
	 * this callback only schedules its decision. */
	const int64_t nowMicroseconds = bu_gettime();
	this->d->lodCompactionPolicy.requestImmediate(nowMicroseconds);
	this->markProgressiveWorkPending();
	this->notifyFrameRequest("lod-resource-pressure");
    }
    const size_t activeCalibrationFaces = calibrationState ?
	calibrationState->activeFaceCount() : 0;
    const size_t activeCalibrationCost = calibrationState ?
	calibrationState->activeRenderCost() : 0;
    size_t calibrationCost = activeCalibrationCost;
    const auto costForPresentedFaces = [activeCalibrationFaces,
	activeCalibrationCost](size_t presentedFaces) {
	if (!activeCalibrationFaces || !activeCalibrationCost)
	    return presentedFaces;
	const long double scaled =
	    static_cast<long double>(activeCalibrationCost) *
	    static_cast<long double>(presentedFaces) /
	    static_cast<long double>(activeCalibrationFaces);
	return scaled >= static_cast<long double>(SIZE_MAX) ? SIZE_MAX :
	    std::max(presentedFaces, static_cast<size_t>(scaled));
    };
    size_t presentedCadPrimitives = 0;
    const SbBool measuredCadPresentation = calibrationState &&
	calibrationState->lastCadPresentedPrimitiveCount(
	    presentedCadPrimitives);
    if (measuredCadPresentation)
	calibrationCost = costForPresentedFaces(presentedCadPrimitives);
    size_t presentedCadRenderCost = 0;
    const SbBool measuredCadRenderCost = calibrationState &&
	calibrationState->lastCadPresentedRenderCost(
	    presentedCadRenderCost);
    if (measuredCadRenderCost)
	calibrationCost = presentedCadRenderCost;
    size_t gpuCadFaces = 0;
    uint64_t gpuCadNanoseconds = 0;
    uint64_t gpuCadSampleSerial = 0;
    float gpuCadPointProxyPixelThreshold = 1.0f;
    const SbBool measuredGpuCadPresentation = calibrationState &&
	calibrationState->lastCadGpuMeasurement(
	    gpuCadFaces, gpuCadNanoseconds, gpuCadSampleSerial,
	    gpuCadPointProxyPixelThreshold);
    const SbBool newGpuCadMeasurement =
	measuredGpuCadPresentation &&
	gpuCadSampleSerial != this->d->lodLastCadGpuSampleSerial;
    if (newGpuCadMeasurement) {
	this->d->lodLastCadGpuSampleSerial = gpuCadSampleSerial;
	this->d->lodLastCadGpuTimeNanoseconds = gpuCadNanoseconds;
	calibrationCost = costForPresentedFaces(gpuCadFaces);
    }
    const SbBool haveCadPresentationAssemblies = calibrationState &&
	calibrationState->hasCadPresentationAssemblies();
    /* Admission budgets are denominated in BObolViewLodState's indexed-mesh
     * render-cost units.  SoCADAssembly's completed-work record intentionally
     * describes backend work instead: on the ordinary OSMesa route its
     * position count may be the expanded submitted stream.  Treating that
     * backend counter as allocator currency made a 75 ms Hubble frame appear
     * to contain nearly twice its retained cost; each recovery then removed
     * only a few entries and repeated for hundreds of frames.
     *
     * When no presentation-only PoP ceiling is hiding the retained cut, the
     * exact active population is the authoritative cost in allocator units.
     * A point threshold of one pixel is part of the normal presentation
     * contract, not a different retained allocation.  Keep the expanded work
     * record for diagnostics and GPU/backend analysis, but pair elapsed time
     * with activeCalibrationCost for admission calibration. */
    const SbBool exactRetainedPopulation =
	haveCadPresentationAssemblies && activeCalibrationCost > 0 &&
	calibrationState->lastCadPresentationFrameExact() &&
	this->d->lodInteractiveProgressiveCeiling < 0 &&
	this->d->lodPresentationPointProxyPixelThreshold <= 1.01f;
    if (exactRetainedPopulation) {
	calibrationCost = activeCalibrationCost;
	this->d->lodBudgetPolicy.noteRetainedQualityFloorMet();
    }
    /* activeRenderCost is retained demand, not necessarily work submitted by
     * this completed traversal.  An empty/replaced CAD plan can legitimately
     * publish zero exact primitives for one frame while the prior occurrence
     * population is still active.  Dividing that retained population by the
     * tiny empty-frame time creates fictitious throughput and reopens a cut
     * which just missed its deadline. */
    const SbBool exactCalibrationWork =
	!haveCadPresentationAssemblies || measuredCadRenderCost ||
	newGpuCadMeasurement;
    this->d->lodLastRenderWasPreparedCadReplay =
	!measuredCadPresentation ||
	calibrationState->lastCadPresentationUsedPreparedReplay();
    const auto addUploadBytes = [](uint64_t left, uint64_t right) {
	return right > UINT64_MAX - left ? UINT64_MAX : left + right;
    };
    uint64_t geometryUploadBytes = 0;
    if (haveFrameGpuResources) {
	geometryUploadBytes = addUploadBytes(geometryUploadBytes,
	    frameGpuResources.ordinaryPartFullUploadBytes);
	geometryUploadBytes = addUploadBytes(geometryUploadBytes,
	    frameGpuResources.ordinaryPartSuffixUploadBytes);
	geometryUploadBytes = addUploadBytes(geometryUploadBytes,
	    frameGpuResources.triangleAtlasFullUploadBytes);
	geometryUploadBytes = addUploadBytes(geometryUploadBytes,
	    frameGpuResources.triangleAtlasSuffixUploadBytes);
    }
    const SbBool geometryUploadChanged = haveFrameGpuResources &&
	(!this->d->lodLastCadGpuGeometryUploadBytesValid ||
	 geometryUploadBytes != this->d->lodLastCadGpuGeometryUploadBytes) ?
	    TRUE : FALSE;
    if (haveFrameGpuResources) {
	this->d->lodLastCadGpuGeometryUploadBytes = geometryUploadBytes;
	this->d->lodLastCadGpuGeometryUploadBytesValid = TRUE;
    }
    /* Prepared multi-draw replay is one reusable presentation path, but it
     * is not the only one.  OSMesa deliberately uses retained ordinary VBOs;
     * once their cumulative upload counters stop changing, an unchanged
     * traversal is an equally valid sustainable-throughput sample. */
    this->d->lodLastRenderWasReusableCadPresentation =
	this->d->lodLastRenderWasPreparedCadReplay ||
	(haveFrameGpuResources && !geometryUploadChanged) ? TRUE : FALSE;
    uint64_t calibrationElapsed =
	newGpuCadMeasurement ? gpuCadNanoseconds : elapsed;
    const SbBool gpuPointProxySampleMatchesCurrent =
	newGpuCadMeasurement &&
	std::fabs(gpuCadPointProxyPixelThreshold -
	    this->d->lodPresentationPointProxyPixelThreshold) <= 0.01f ?
	    TRUE : FALSE;
    /* Timer queries complete on a later render.  Use them for the point-cut
     * bracket only when the query carries the same threshold as the cut being
     * evaluated; otherwise the just-completed CPU traversal is the current
     * evidence.  Scene throughput still uses the query's paired face count. */
    const uint64_t pointCalibrationElapsed =
	gpuPointProxySampleMatchesCurrent ? gpuCadNanoseconds : elapsed;
    /* Host presentation cadence is deliberately not renderer capacity.  A
     * slow mouse, scripted pause, compositor cap, or event-loop scheduling
     * gap cannot be improved by discarding geometry.  CPU traversal duration
     * and the asynchronous GL timer above are the actionable work samples. */
    /* Result publication may allocate/upload newly decoded geometry during
     * this traversal.  That one-time cost is relevant to the cooldown below,
     * but not to the sustainable FPS of the retained cut.  Quiet,
     * unchanged calibration probes measure the latter.  Interactive samples
     * keep the full cost because upload latency is part of the user's motion
     * experience. */
    const SbBool transientStablePublication =
	this->d->lodPublicationPolicy.framePending() &&
	!this->d->lodInteractive && !this->d->lodGestureActive;
    const SbBool transientStablePreparation =
	!this->d->lodLastRenderWasReusableCadPresentation &&
	!this->d->lodInteractive && !this->d->lodGestureActive;
    const SbBool reusableCadWorkSample =
	this->d->lodLastRenderWasReusableCadPresentation &&
	((haveCadPresentationAssemblies && measuredCadRenderCost) ||
	 (!haveCadPresentationAssemblies && activeCalibrationCost > 0)) ?
	    TRUE : FALSE;
    if (calibrationCost > 0 && calibrationElapsed > 0 &&
	exactCalibrationWork &&
	(newGpuCadMeasurement ||
	 (!transientStablePublication && !transientStablePreparation))) {
	const long double sample =
	    static_cast<long double>(calibrationCost) * 1000000000.0L /
	    static_cast<long double>(calibrationElapsed);
	if (std::isfinite(sample) && sample > 0.0L) {
	    const SbBool interactiveSample =
		this->d->lodInteractive || this->d->lodGestureActive;
	    long double &calibration = interactiveSample ?
		this->d->lodInteractiveCalibratedRenderCostPerSecond :
		this->d->lodStableCalibratedRenderCostPerSecond;
	    if (calibration <= 0.0L) {
		calibration = sample;
	    } else if (interactiveSample && sample < calibration) {
		/* Missing an interaction target is immediately visible to the
		 * user.  Adopt the observed lower ceiling immediately in that
		 * case.  A fast underloaded cut is not a lower throughput proof:
		 * fixed frame overhead divided into 44 or 118 triangles produces
		 * an arbitrarily small faces/second number and used to lock OSMesa
		 * permanently at that cut despite ample frame headroom. */
		const long double targetNanoseconds =
		    this->d->lodInteractiveTargetFps > 0.0f ?
		    1000000000.0L /
			static_cast<long double>(
			    this->d->lodInteractiveTargetFps) : 0.0L;
		const long double priorFrameBudget =
		    this->d->lodInteractiveTargetFps > 0.0f ?
		    calibration /
			static_cast<long double>(
			    this->d->lodInteractiveTargetFps) : 0.0L;
		const SbBool materiallyLoaded =
		    priorFrameBudget <= 0.0L ||
		    static_cast<long double>(calibrationCost) >=
			priorFrameBudget * 0.10L;
		if (targetNanoseconds > 0.0L &&
		    materiallyLoaded &&
		    static_cast<long double>(calibrationElapsed) >
			targetNanoseconds)
		    calibration = sample;
	    } else if (!interactiveSample && newGpuCadMeasurement &&
		sample < calibration) {
		/* A completed timer query measures the driver's actual queued CAD
		 * work rather than CPU command submission or an unrelated window
		 * composition delay.  A lower GPU ceiling is therefore actionable
		 * immediately; damping it recreates several known-bad long frames
		 * while the estimate slowly catches up. */
		calibration = sample;
	    } else if (!interactiveSample && sample < calibration) {
		/* A quiet event-driven view includes occasional screenshots,
		 * buffer readbacks, window composition, and cache activity.  One
		 * such outlier must not make admitted leaves disappear and then
		 * return on the next ordinary frame.  Reduce stable capacity only
		 * after a frame materially misses its target, and then do so with
		 * damping. */
		const long double targetNanoseconds =
		    this->d->quietAllocationTargetFps() > 0.0f ?
		    1000000000.0L /
			static_cast<long double>(
			    this->d->quietAllocationTargetFps()) : 0.0L;
		if (targetNanoseconds > 0.0L &&
		    static_cast<long double>(calibrationElapsed) >
			targetNanoseconds * 1.20L)
		    calibration = calibration * 0.80L + sample * 0.20L;
	    } else {
		const long double oldWeight = interactiveSample ? 0.98L : 0.90L;
		calibration =
		    calibration * oldWeight + sample * (1.0L - oldWeight);
	    }

	    /* A quality probe is a scale-changing presentation, but an
	     * completed frame which already meets the quieter stable deadline is
	     * a strict lower-bound proof for stable admission.  This remains true
	     * when that frame extended a resident VBO: total elapsed time includes
	     * the one-time upload, making the proof more conservative, and the
	     * completed prefix is resident afterward.  Transfer only exact
	     * submitted work: an occurrence may retain a much richer prefix behind
	     * the renderer ceiling, so activeRenderCost is not a safe substitute
	     * on compatibility renderers.  The scalar policy owns the deadline and
	     * throughput arithmetic. */
	    /* A triangle count is exact in face units but cannot safely stand in
	     * for the scheduler's weighted point/normal/occurrence units.  Direct
	     * renderers publish the complete logical work record for this
	     * cross-mode proof.  GPU-only tiers remain on their ordinary
	     * mode-local calibration until they publish the same record. */
	    const bool exactPresentedCost = measuredCadRenderCost;
	    const bool exactScaleQualityFrame =
		interactiveSample &&
		this->d->lodViewDemandPolicy.scaleChangingInteraction(
		    interactiveSample) &&
		exactPresentedCost;
	    int presentedMaximum = calibrationState ?
		calibrationState->maximumActiveProgressiveCut() : -1;
	    if (this->d->lodInteractiveProgressiveCeiling >= 0 &&
		presentedMaximum >= 0)
		presentedMaximum = std::min(presentedMaximum,
		    this->d->lodInteractiveProgressiveCeiling);
	    /* Stable throughput targets the configured quiet cadence, but the
	     * restoration witness has a different contract: it is the richest
	     * exact scale-quality frame which completed under the explicit 100 ms
	     * hard bound.  Conflating the two remembered an earlier 60 Hz motion
	     * cut and discarded a visibly richer, deadline-safe frame as soon as
	     * wheel input became quiet.  The richer cut is already resident and
	     * was actually presented; retaining it does not teach the 20 Hz
	     * allocator an inflated throughput. */
	    if (exactScaleQualityFrame && elapsed > 0 &&
		elapsed <= BObolLodViewDemandPolicy::
		    qualityFrameDurationNanoseconds())
		this->d->lodPresentationPolicy.noteStableQuality(
		    presentedMaximum,
		    this->d->lodPresentationPointProxyPixelThreshold,
		    controller_lod_presentation_population(calibrationState),
		    this->d->lodViewRevision, calibrationCost);
	    const BObolLodQualityPolicy::StableCapacityEvidence
		stableEvidence = BObolLodQualityPolicy::stableCapacityEvidence(
		    calibrationCost, elapsed,
		    this->d->quietAllocationTargetFps(),
		    exactScaleQualityFrame);
		    if (stableEvidence.proven) {
			this->d->lodStableCalibratedRenderCostPerSecond = std::max(
			    this->d->lodStableCalibratedRenderCostPerSecond,
			    stableEvidence.renderCostPerSecond);
		    }
	}
    }

    /* A terminal budget pass explicitly requests an unchanged replay for its
     * late-headroom witness.  A retained VBO may still perform its one-time
     * upload or command preparation on that first replay; permit a strictly
     * bounded successor replay in that case, then consume only the first
     * reusable frame.  Opportunistic reuse of a later selection/HUD repaint
     * made a view report STABLE and then restart LoD work under unrelated user
     * interaction. */
    if (this->d->lodHeadroomPolicy.retryPending()) {
	const bool stableContext =
	    !this->d->lodInteractive && !this->d->lodGestureActive &&
	    this->d->lodBudgetPolicy.stableBudgetLimited() &&
	    !this->d->lodSubmissionPending &&
	    !this->d->lodBudgetPolicy.rescanAfterFrame() &&
	    !this->d->lodRefinementAwaitingFrame &&
	    this->d->quietAllocationTargetFps() > 0.0f &&
	    this->d->lodStableCalibratedRenderCostPerSecond > 0.0L &&
	    this->d->lodBudgetPolicy.currentBudget() != SIZE_MAX;
	const long double demonstratedBudget = stableContext ?
	    this->d->lodStableCalibratedRenderCostPerSecond /
		static_cast<long double>(
		    this->d->quietAllocationTargetFps()) : 0.0L;
	const uint64_t stableTargetNanoseconds = stableContext ?
	    static_cast<uint64_t>(1000000000.0L /
		static_cast<long double>(
		    this->d->quietAllocationTargetFps())) : 0;
	const bool matchingHeadroomWitness =
	    this->d->lodHeadroomPolicy.pendingMatches(
		this->d->lodViewRevision, this->d->lodPolicyRevision,
		activeCalibrationCost);
	const size_t discreteTrialExcess =
	    BObolLodQualityPolicy::discreteTrialOverBudgetAllowance(
		activeCalibrationCost,
		this->d->lodBudgetPolicy.currentBudget());
	const bool stableDiscreteTrial =
	    matchingHeadroomWitness && stableContext &&
	    measuredCadRenderCost && reusableCadWorkSample &&
	    !transientStablePublication && !transientStablePreparation &&
	    calibrationCost == activeCalibrationCost &&
	    stableTargetNanoseconds > 0 &&
	    elapsed <= stableTargetNanoseconds && discreteTrialExcess > 0;
	const bool transientHeadroomReplay =
	    matchingHeadroomWitness && stableContext &&
	    measuredCadRenderCost &&
	    (transientStablePublication || transientStablePreparation);
	const bool deferredHeadroomReplay = transientHeadroomReplay &&
	    this->d->lodHeadroomPolicy.deferTransientReplay(
		this->d->lodViewRevision, this->d->lodPolicyRevision,
		activeCalibrationCost);
	const bool reusableHeadroom = !deferredHeadroomReplay &&
	    this->d->lodHeadroomPolicy.consumeRetry(
		this->d->lodViewRevision, this->d->lodPolicyRevision,
		this->d->lodBudgetPolicy.currentBudget(), activeCalibrationCost,
		demonstratedBudget,
		calibrationElapsed, stableTargetNanoseconds,
		stableContext && reusableCadWorkSample &&
		    !transientStablePublication &&
		    !transientStablePreparation);
	if (getenv("BOBOL_LOD_TRACE_HEADROOM"))
	    bu_log("BObol LoD headroom complete stable=%d matching=%d "
		   "measured=%d reusable=%d transient_publication=%d "
		   "transient_preparation=%d deferred=%d admitted=%d "
		   "discrete=%d active_cost=%zu calibration_cost=%zu "
		   "budget=%zu demonstrated=%.3Lf elapsed_ms=%.3f "
		   "target_ms=%.3f view=%llu policy=%llu\n",
		   stableContext ? 1 : 0,
		   matchingHeadroomWitness ? 1 : 0,
		   measuredCadRenderCost ? 1 : 0,
		   reusableCadWorkSample ? 1 : 0,
		   transientStablePublication ? 1 : 0,
		   transientStablePreparation ? 1 : 0,
		   deferredHeadroomReplay ? 1 : 0,
		   reusableHeadroom ? 1 : 0,
		   stableDiscreteTrial ? 1 : 0,
		   activeCalibrationCost, calibrationCost,
		   this->d->lodBudgetPolicy.currentBudget(),
		   demonstratedBudget,
		   calibrationElapsed / 1000000.0,
		   stableTargetNanoseconds / 1000000.0,
		   (unsigned long long)this->d->lodViewRevision.value(),
		   (unsigned long long)this->d->lodPolicyRevision.value());
	if (deferredHeadroomReplay) {
	    this->markProgressiveWorkPending();
	    this->requestRender("lod-calibrated-headroom-replay");
	    this->notifyFrameRequest("lod-calibrated-headroom-replay");
	} else if (reusableHeadroom || stableDiscreteTrial) {
	    this->d->lodBudgetPolicy.clearBudgetLimit();
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionPending = TRUE;
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodBudgetPolicy.resetPass();
	    /* The replay proved a different sustainable allowance.  Recompute one
	     * complete screen-importance allocation with that new capacity even
	     * when the old active population is already below it.  An ordinary
	     * first-come pass can strand the unused remainder on a discrete PoP
	     * population and makes cold and warm caches converge differently. */
	    this->d->lodBudgetPolicy.requestRetainedReallocation(false);
	    this->d->lodDiscretePopulationTrialAvailable =
		stableDiscreteTrial ? TRUE : FALSE;
	    this->markProgressiveWorkPending();
	    this->requestRender(stableDiscreteTrial ?
		"lod-discrete-headroom" : "lod-calibrated-headroom");
	}
    }

    /*
     * Consume a quiet small-part aggregation probe without invoking another
     * O(N) provider/admission scan.  If the measured reusable presentation
     * still misses the stable target, ratchet the camera-local threshold and
     * request exactly one successor frame.  The immutable per-part geometry
     * and desired PoP cuts are untouched throughout.
    */
    /* Point aggregation is implemented by SoCADAssembly's occurrence batch.
     * Ordinary retained VBO payloads have no point-cut population to
     * calibrate.  Waiting for one on that route creates an unwitnessable
     * progressive barrier even though all renderer/service work is idle.
     *
     * An assembly may also remain retained after every occurrence under its
     * drawn root has been hidden (for example, `erase all/r.stl`).  The latest
     * complete coverage census is the authority in that case: a zero visible
     * population has no point threshold to measure.  Retire the old bracket
     * and restore the one-pixel default while the scene is empty.  Otherwise
     * an interrupted pre-erase software frame can leave the controller asking
     * an exact zero-work frame for a CAD work record which cannot exist. */
    const SbBool noVisibleCadPopulation =
	this->d->lodCoveragePolicy.hasCompleteVisibleCount() &&
	this->d->lodCoveragePolicy.completeVisibleCount() == 0 &&
	!this->d->lodCoveragePolicy.active() &&
	!this->d->lodSubmissionPending ? TRUE : FALSE;
    if (this->d->lodStablePointProxyCalibrationPending &&
	noVisibleCadPopulation) {
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->lodPointProxyCalibrationPolicy.reset();
	BObolViewLodState *presentationState = this->d->viewAttachment ?
	    this->d->viewAttachment->getViewLodState() : NULL;
	if (presentationState)
	    presentationState->setCadPresentationPointProxyPixelThreshold(1.0f);
	(void)this->d->lodBudgetPolicy.
	    confirmRetainedRecoveryPresentation(true);
    }
    if (this->d->lodStablePointProxyCalibrationPending &&
	!haveCadPresentationAssemblies) {
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	(void)this->d->lodBudgetPolicy.
	    confirmRetainedRecoveryPresentation(true);
    }
    if (this->d->lodStablePointProxyCalibrationPending) {
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	const long double stableTargetNanoseconds =
	    this->d->quietAllocationTargetFps() > 0.0f ?
		1000000000.0L /
		    static_cast<long double>(
			this->d->quietAllocationTargetFps()) : 0.0L;
	const SbBool reusableStableSample =
	    !this->d->lodInteractive &&
	    !this->d->lodGestureActive &&
	    reusableCadWorkSample &&
	    !transientStablePublication &&
	    !transientStablePreparation &&
	    stableTargetNanoseconds > 0.0L &&
	    pointCalibrationElapsed > 0;
	const size_t minimumCalibrationCost = calibrationState ?
	    calibrationState->minimumActiveRenderCost() :
	    activeCalibrationCost;
	const SbBool reducibleProgressiveDetail =
	    activeCalibrationCost > minimumCalibrationCost ? TRUE : FALSE;
	const SbBool stableSampleOverloaded =
	    reusableStableSample &&
	    static_cast<long double>(pointCalibrationElapsed) >
		stableTargetNanoseconds * 1.10L ? TRUE : FALSE;
	const SbBool coarsePointCut =
	    this->d->lodPresentationPointProxyPixelThreshold > 1.01f ?
		TRUE : FALSE;
	SbBool scheduledTriangleRecovery = FALSE;
	SbBool restoredOnePixelCut = FALSE;

	/* Resolve the two quality dimensions from measured work.  A multi-pixel
	 * point cut represents only occurrences below that screen-space size; it
	 * is not permission to discard detail from large arrays, wheels, blades,
	 * tails, or other prominent parts.  The former coarsePointCut condition
	 * collapsed every reducible prefix even after a safe aggregated frame,
	 * restored the one-pixel population, and then repeated forever.  Recover
	 * triangle capacity only from an exact overloaded frame.  Otherwise let
	 * the point bracket converge while the importance allocator spends the
	 * retained triangle budget on visible multi-pixel geometry.
	 */
	if (reusableStableSample &&
	    BObolLodPointProxyCalibrationPolicy::shouldRecoverTriangleDetail(
		reducibleProgressiveDetail != FALSE,
		stableSampleOverloaded != FALSE,
		coarsePointCut != FALSE,
		this->d->lodBudgetPolicy.retainedQualityFloorActive() ||
		    this->d->lodBudgetPolicy.retainedQualityFloorRejected())) {
	    /* recoveryBudget is consumed by the retained allocator, so scale its
	     * own active-cost population by the measured deadline ratio.  Backend
	     * submitted-work counts are not dimensionally interchangeable with
	     * this budget (notably for OSMesa's expanded position stream). */
	    long double affordable =
		static_cast<long double>(activeCalibrationCost) *
		stableTargetNanoseconds * 0.80L /
		static_cast<long double>(pointCalibrationElapsed);
	    if (!std::isfinite(affordable) || affordable <= 0.0L)
		affordable = static_cast<long double>(minimumCalibrationCost);
	    size_t recoveryBudget = affordable >=
		    static_cast<long double>(SIZE_MAX) ? SIZE_MAX :
		static_cast<size_t>(affordable);
	    recoveryBudget = std::max(
		minimumCalibrationCost, recoveryBudget);
	    if (recoveryBudget < activeCalibrationCost) {
		this->d->lodBudgetPolicy.requestRetainedRecovery(
		    recoveryBudget);
		this->d->lodBudgetPolicy.resetProbeSeries();
		this->d->lodBudgetPolicy.resetOverloadRecovery();
		this->d->lodBudgetPolicy.resetPass();
		this->d->lodHeadroomPolicy.cancelRetry();
		/* The handoff's cost floor was proven with the coarse point cut
		 * which hid this population.  It is not evidence for the one-pixel
		 * hierarchy and would immediately re-admit the detail being
		 * compacted here.  Retain the renderer ceiling itself until the
		 * coherent recovery frame below, but retire the stale proof. */
		this->d->lodPresentationPolicy.cancelHandoff();
		this->d->lodLastSubmittedViewRevision.reset();
		this->d->lodLastSubmittedPolicyRevision.reset();
		this->d->lodSubmissionSourceIndex = 0;
		this->d->lodSubmissionEntryOffset = 0;
		this->d->clearLodSubmissionPlan();
		this->d->lodSubmissionRescanPending = FALSE;
		this->d->lodSubmissionPending = TRUE;
		this->d->lodRetainedRefinementPending = FALSE;
		this->d->lodRetainedRefinementCutAdvanced = FALSE;
		this->d->lodRetainedRefinementBudgetBlocked = FALSE;
		this->d->lodPassAdmittedWork = FALSE;
		this->d->lodPointProxyTriangleRecoveryPending = TRUE;
		this->markProgressiveWorkPending();
		scheduledTriangleRecovery = TRUE;
	    } else if (coarsePointCut) {
		/* The complete retained cut already fits the demonstrated capacity.
		 * Test the one-pixel population directly; no provider or admission
		 * traversal is necessary. */
		this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
		this->d->lodPointProxyCalibrationPolicy.reset();
		BObolViewLodState *presentationState =
		    this->d->viewAttachment->getViewLodState();
		if (presentationState)
		    presentationState->
			setCadPresentationPointProxyPixelThreshold(1.0f);
		this->d->lodStablePointProxyCalibrationPending = TRUE;
		this->requestRender("lod-stable-point-restore");
		this->notifyFrameRequest("lod-stable-point-restore");
		restoredOnePixelCut = TRUE;
	    }
	}

	const BObolLodPointProxyCalibrationPolicy::Decision calibration =
	    (!scheduledTriangleRecovery && !restoredOnePixelCut) ?
	    this->d->lodPointProxyCalibrationPolicy.completed(
		this->d->lodPresentationPointProxyPixelThreshold,
		pointCalibrationElapsed,
		this->d->quietAllocationTargetFps(),
		reusableStableSample != FALSE) :
	    BObolLodPointProxyCalibrationPolicy::Decision();
	if (!scheduledTriangleRecovery && !restoredOnePixelCut &&
	    calibration.changed) {
	    this->d->lodPresentationPointProxyPixelThreshold =
		calibration.threshold;
	    BObolViewLodState *presentationState =
		this->d->viewAttachment->getViewLodState();
	    if (presentationState)
		presentationState->
		    setCadPresentationPointProxyPixelThreshold(
			calibration.threshold);
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	    this->requestRender("lod-stable-point-calibration");
	    this->notifyFrameRequest(
		"lod-stable-point-calibration");
	} else if (!scheduledTriangleRecovery && !restoredOnePixelCut &&
	    !reusableStableSample &&
	    !this->d->lodInteractive &&
	    !this->d->lodGestureActive &&
	    haveCadPresentationAssemblies &&
	    stableTargetNanoseconds > 0.0L &&
	    this->d->lodPointProxyCalibrationPolicy.
		requiresReusableConfirmation(
		    this->d->lodPresentationPointProxyPixelThreshold)) {
	    /*
	     * Applying a different point threshold changes the retained
	     * occurrence classification.  Its first frame may rebuild/patch the
	     * prepared plan, publish worker results, or follow an interrupted
	     * traversal and therefore lack a reusable work record.  None of those
	     * cases is evidence that calibration converged.  Preserve the
	     * obligation until an unchanged replay measures the cut it actually
	     * drew.
	     *
	     * Result streaming may extend this handoff across several frames.
	     * Those frames already represent real pending work; retaining the
	     * idempotent frame request neither exposes a richer cut prematurely nor
	     * permits the view to report stable at an unmeasured emergency cut.
	     */
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	    this->requestRender("lod-stable-point-replay");
	    this->notifyFrameRequest("lod-stable-point-replay");
	}
	/* A triangle-prefix recovery which temporarily used a coarser aggregate
	 * point cut reaches this branch twice: once to request the one-pixel frame,
	 * and once after that frame completes.  The latter exact presentation is
	 * the missing transition which retires the measured recovery ceiling. */
	if (!scheduledTriangleRecovery && !restoredOnePixelCut &&
	    exactRetainedPopulation &&
	    this->d->lodPresentationPointProxyPixelThreshold <= 1.01f)
	    this->resumeLodAfterOnePixelRecovery();
    }

    /* A retained-prefix recovery is complete only after its selected cut has
     * been presented and no follow-up admission/handoff barrier remains.
     * Restore the one-pixel occurrence contract now and measure that exact
     * population on the successor frame. */
    if (this->d->lodPointProxyTriangleRecoveryPending &&
	!this->d->lodInteractive && !this->d->lodGestureActive &&
	!this->d->lodSubmissionPending &&
	!this->d->lodSubmissionRescanPending &&
	!this->d->lodRetainedRefinementPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodPresentationPolicy.handoffPending()) {
	this->d->lodPointProxyTriangleRecoveryPending = FALSE;
	const SbBool pointCutChanged =
	    this->d->lodPresentationPointProxyPixelThreshold > 1.01f ?
		TRUE : FALSE;
	this->d->lodInteractiveProgressiveCeiling = -1;
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->lodPointProxyCalibrationPolicy.reset();
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	if (presentationState)
	    presentationState->setCadPresentationProgressiveCutCeiling(-1);
	if (presentationState)
	    presentationState->setCadPresentationPointProxyPixelThreshold(1.0f);
	/* Point calibration measures a changed occurrence population.  A retained
	 * triangle-prefix recovery commonly starts and ends at the one-pixel cut;
	 * arming the point latch in that case asks an unchanged ordinary OSMesa
	 * VBO frame for CAD-batch evidence it cannot publish and can strand an
	 * otherwise idle view forever. */
	this->d->lodStablePointProxyCalibrationPending =
	    pointCutChanged && presentationState &&
	    presentationState->hasCadPresentationAssemblies() ? TRUE : FALSE;
	if (this->d->lodStablePointProxyCalibrationPending) {
	    this->requestRender("lod-stable-point-restore");
	    this->notifyFrameRequest("lod-stable-point-restore");
	} else {
	    /* The recovery ceiling prevents the cheaper preparation frame from
	     * immediately re-admitting the population which just missed its
	     * deadline.  It is a one-witness guard, not a permanent fidelity cap.
	     * Once the coherent one-pixel recovery cut has actually rendered, its
	     * measured throughput must be allowed to grow the budget again in
	     * bounded, screen-priority waves.  Keeping the ceiling indefinitely
	     * made a transient 84 ms Hubble frame settle forever at the 35k-face
	     * minimum even though subsequent 16 ms frames proved ample headroom. */
	    this->resumeLodAfterOnePixelRecovery();
	}
    }

    /*
     * A cold or newly exposed many-leaf stream can spend seconds completing
     * its coverage pass.  Waiting until that pass ends before applying the
     * renderer's measured small-part cut leaves every intervening software
     * frame at 100+ ms and makes console/input handling feel stalled.
     *
     * Result publication time is deliberately excluded from sustainable
     * throughput calibration above, but it is still real user-visible frame
     * latency.  Under severe pressure, raise the presentation-only
     * aggregation cut for the next frame while coverage continues.  The
     * final quiet probe later calibrates it precisely and no retained mesh or
     * cache state is discarded.
     */
    const long double ongoingStableTargetNanoseconds =
	this->d->quietAllocationTargetFps() > 0.0f ?
	    1000000000.0L /
		static_cast<long double>(
		    this->d->quietAllocationTargetFps()) : 0.0L;
    const SbBool ongoingStableWork =
	this->d->lodSubmissionPending ||
	this->hasPendingLodResults() ||
	(this->d->lodService &&
	 (this->d->lodService->activeTaskCountForGeneration(
	      this->d->lodActiveGeneration) > 0 ||
	  this->d->lodService->queuedResultCountForGeneration(
	      this->d->lodActiveGeneration) > 0));
    if (!this->d->lodStablePointProxyCalibrationPending &&
	haveCadPresentationAssemblies &&
	!this->d->lodInteractive &&
	!this->d->lodGestureActive &&
	this->d->lodAutoSubmit &&
	ongoingStableWork &&
	measuredCadPresentation &&
	ongoingStableTargetNanoseconds > 0.0L &&
	static_cast<long double>(elapsed) >
	    ongoingStableTargetNanoseconds * 1.50L) {
	const BObolLodPointProxyCalibrationPolicy::Decision pressure =
	    this->d->lodPointProxyCalibrationPolicy.interrupted(
		this->d->lodPresentationPointProxyPixelThreshold,
		elapsed, this->d->quietAllocationTargetFps());
	if (pressure.changed) {
	    this->d->lodPresentationPointProxyPixelThreshold =
		pressure.threshold;
	    BObolViewLodState *presentationState =
		this->d->viewAttachment->getViewLodState();
	    if (presentationState)
		presentationState->
		    setCadPresentationPointProxyPixelThreshold(
			pressure.threshold);
	    /* The changed SoCADAssembly occurrence population needs one exact
	     * successor frame before the bracket can be adjusted again. */
	    this->d->lodStablePointProxyCalibrationPending =
		this->d->lodPresentationPointProxyPixelThreshold > 1.01f ?
		    TRUE : FALSE;
	    /* This pressure step is deliberately conservative because its frame
	     * includes streaming/publication latency.  The unchanged successor
	     * frame will walk back toward the finest sustainable presentation when
	     * it demonstrates headroom. */
	    this->requestRender("lod-stream-point-pressure");
	    this->notifyFrameRequest("lod-stream-point-pressure");
	}
    }

    /* Point aggregation and the stable-presentation handoff can consume the
     * frame which followed the terminal budget pass.  Revisit the explicit
     * late-headroom witness immediately after those barriers have either
     * converged or re-armed themselves.  This preserves an already resident
     * rich zoom cut and lets a quiet view continue from it instead of
     * compacting to the conservative seed budget. */
    this->armStableLodHeadroomProbeIfReady();

    /* A zoom-quality probe is allowed to spend more than an ordinary motion
     * frame, but never more than the 10 Hz hard responsiveness floor.  If a
     * newly exposed PoP cut misses that floor, correct the next frame with
     * the renderer-wide prefix ceiling.  This is an O(1), reversible
     * presentation change: the richer occurrence cut and its resident bytes
     * stay intact, so continued zoom can refine from them and quiet handoff
     * does not have to walk back up from the minimum mesh.
     *
     * Cut ordinals do not encode a population ratio.  Back off exactly one
     * producer-admissible cut after a proven miss and let the following
     * completed frame decide whether another step is necessary. */
    if ((this->d->lodInteractive || this->d->lodGestureActive) &&
	this->d->lodViewDemandPolicy.interactionScaleChanged() &&
	this->d->lodViewDemandPolicy.qualityProbeActive() &&
	elapsed >
	    BObolLodViewDemandPolicy::qualityFrameDurationNanoseconds()) {
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	const int activeMaximum = presentationState ?
	    presentationState->maximumActiveProgressiveCut() : -1;
	if (activeMaximum > 0) {
	    const int presentedMaximum =
		this->d->lodInteractiveProgressiveCeiling >= 0 ?
		    std::min(activeMaximum,
			this->d->lodInteractiveProgressiveCeiling) :
		    activeMaximum;
	    const int correctedCeiling =
		std::max(0, presentedMaximum - 1);
	    if (correctedCeiling < presentedMaximum) {
		this->d->lodViewDemandPolicy.noteQualityMiss(
		    correctedCeiling,
		    measuredCadRenderCost ? presentedCadRenderCost :
			(activeCalibrationCost ? activeCalibrationCost :
			 presentedCadPrimitives));
		this->d->lodInteractiveProgressiveCeiling =
		    correctedCeiling;
		presentationState->setCadPresentationProgressiveCutCeiling(
		    correctedCeiling);
		this->requestRender("lod-scale-quality-pressure");
		this->notifyFrameRequest("lod-scale-quality-pressure");
	    }
	}
    }

    /* A budget-limited pass may have admitted a bounded first wave and then
     * scanned the remaining boxes without admitting them.  Re-plan from the
     * highest screen-value occurrence after that wave has actually rendered
     * and supplied a throughput sample.  Replanning before the frame would
     * repeatedly spend the same stale allowance. */
    if (this->d->lodBudgetPolicy.rescanAfterFrame()) {
	/*
	 * Calibration probes measure an unchanged retained cut.  They do not
	 * need an intervening O(N) occurrence-planning pass: that pass cannot
	 * admit anything until the samples have changed the aggregate budget.
	 * Present the bounded probe series back-to-back, then scan the sparse
	 * unsatisfied frontier once using the resulting calibration.
	 */
	const BObolLodBudgetPolicy::CompletedFrameDecision calibration =
	    this->d->lodBudgetPolicy.completeCalibrationFrame();
	if (calibration.requestCalibrationFrame) {
	    this->requestRender("lod-budget-calibration");
	    this->notifyFrameRequest("lod-budget-calibration");
	} else if (calibration.restartSubmission) {
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionPending = TRUE;
	    this->d->lodSubmissionRescanPending = FALSE;
	    /* A missing bounded coverage census retires its counter pass while
	     * waiting for this measured-frame admission decision.  The calibrated
	     * retry is another complete current-view census, not merely a quality
	     * pass: reactivate its counters so a cold stream which has since
	     * received every minimum prefix can publish coverageComplete(). */
	    if (this->d->lodCoveragePolicy.required() &&
		!this->d->lodCoveragePolicy.coverageComplete() &&
		!this->d->lodCoveragePolicy.active())
		this->d->lodCoveragePolicy.activate(false);
	    if (this->d->lodCoveragePolicy.active())
		this->d->lodCoveragePolicy.clearPassCounters();
	    this->markProgressiveWorkPending();
	}
    }
    /* Once the scale-quality budget is active, every completed frame is a
	 * measurement of the current render ceiling even if a newer wheel event
	 * has already cleared the one-shot probe label.  Preserve that proof
	 * across the next camera epoch; otherwise the ordinary 60 Hz motion
	 * policy can immediately undo a 10 Hz quality cut which was just shown
	 * successfully. */
    int completedPresentedMaximum = calibrationState ?
	calibrationState->maximumActiveProgressiveCut() : -1;
    if (this->d->lodInteractiveProgressiveCeiling >= 0 &&
	completedPresentedMaximum >= 0)
	completedPresentedMaximum = std::min(completedPresentedMaximum,
	    this->d->lodInteractiveProgressiveCeiling);
    this->d->lodViewDemandPolicy.noteFramePresented(
	completedPresentedMaximum,
	measuredCadRenderCost != FALSE && presentedCadRenderCost > 0,
	elapsed);
    if (this->d->lodViewDemandPolicy.rearmAfterQualityFrame(
	    this->d->lodInteractive || this->d->lodGestureActive,
	    calibrationState ?
		calibrationState->maximumActiveProgressiveCut() : -1,
	    completedPresentedMaximum,
	    measuredCadRenderCost != FALSE && presentedCadRenderCost > 0,
	    elapsed))
	this->markProgressiveWorkPending();
    this->completePresentationBarrier(elapsed);

    /* The frame-completion edge above is what retires the progressive-prefix
     * presentation barrier.  A terminal planning pass may already have tried
     * to arm headroom while that barrier was still active; if no other work
     * remains, there will be no later planning pass from which to try again.
     * Recheck after every completed-frame latch has been reconciled so a
     * cheap, newly presented population gets its one bounded richer probe
     * instead of being misreported as the stable final cut.  The headroom
     * policy keys that witness by view, policy, and active population, which
     * keeps an actually capacity-limited population from repainting in a
     * loop. */
    this->armStableLodHeadroomProbeIfReady();
}

void
BObolViewController::completePresentationBarrier(uint64_t elapsedNanoseconds)
{
    this->d->renderCompletionSerial++;
    if (this->d->renderCompletionSerial == 0)
	this->d->renderCompletionSerial++;

    const bool deadlineHandoffReady =
	this->d->lodPresentationPolicy.noteFramePresented();
    if (deadlineHandoffReady) {
	/* The constrained retry frame is the proof which a deadline-created
	 * handoff was waiting for.  The admission pass normally ran before that
	 * frame and correctly left the ceiling armed; schedule the successor pass
	 * now so it can retire the handoff and restore the quiet presentation. */
	this->d->lodLastSubmittedViewRevision.reset();
	this->d->lodLastSubmittedPolicyRevision.reset();
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPending = TRUE;
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	this->markProgressiveWorkPending();
	this->notifyFrameRequest("lod-deadline-handoff-ready");
    }

    /* Rendering time is not user-idle time.  Start the quiet clock only
     * after a frame for the newest camera has actually completed. */
    if (this->d->lodInteractive &&
	this->d->lodSettleAfterRenderSerial != 0 &&
	this->d->renderCompletionSerial >=
	    this->d->lodSettleAfterRenderSerial) {
	this->d->lodLastViewChangeMicroseconds = bu_gettime();
	this->d->lodSettleAfterRenderSerial = 0;
	if (this->d->lodViewDemandPolicy.noteMotionFrameSettled())
	    this->markProgressiveWorkPending();
    }

    /* A partial PoP result must be presented before the next, potentially
     * much larger prefix is requested.  This is presentation sequencing, not
     * capacity evidence: an exact selection/HUD frame which traversed the
     * current population satisfies the barrier just as an LoD-named frame
     * does. */
    if (this->d->lodRefinementAwaitingFrame &&
	this->d->renderCompletionSerial >=
	    this->d->lodRefinementResumeAfterRenderSerial) {
	const SbBool resumeResidentRefinement =
	    this->d->lodRetainedResidencyPending;
	this->d->lodRefinementAwaitingFrame = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	const uint64_t responsiveFrame = 33333334ULL;
	int64_t cooldownMicroseconds = 0;
	if (elapsedNanoseconds > responsiveFrame) {
	    const uint64_t elapsedMicroseconds =
		elapsedNanoseconds / 1000ULL;
	    cooldownMicroseconds = static_cast<int64_t>(
		std::max<uint64_t>(50000ULL,
		    std::min<uint64_t>(2000000ULL, elapsedMicroseconds)));
	}
	const int64_t nowMicroseconds = bu_gettime();
	this->d->lodRefinementNotBeforeMicroseconds =
	    cooldownMicroseconds > 0 &&
	    nowMicroseconds <=
		std::numeric_limits<int64_t>::max() - cooldownMicroseconds ?
	    nowMicroseconds + cooldownMicroseconds : nowMicroseconds;
	this->d->lodLastSubmittedViewRevision.reset();
	this->d->lodLastSubmittedPolicyRevision.reset();
	if (resumeResidentRefinement) {
	    this->d->lodRetainedResidencyPending = FALSE;
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodSubmissionPending = TRUE;
	}
	this->markProgressiveWorkPending();
	this->notifyFrameRequest(resumeResidentRefinement ?
	    "lod-resident-refinement-resume" : "lod-refinement-resume");
    }

    this->d->lodPublicationPolicy.frameCompleted();
    this->scheduleResidentGrowthReallocationIfReady();
    if (this->d->lodViewDemandPolicy.qualityProbePending()) {
	this->markProgressiveWorkPending();
	this->notifyFrameRequest("lod-scale-quality-resume");
    }
    /* synchronizePresentation() ran before this completed frame, so it
     * contains every immutable result accumulated by the owner thread. */
    this->d->reconcilePhase(BObolLodStateMachine::Event::FRAME_COMPLETED);
}

void
BObolViewController::scheduleLodRefinementFrame(const char *reason)
{
    /* This gate is useful only if a host presentation is guaranteed to
     * follow it.  Merely latching requestRender() is insufficient: a result
     * can be drained while progressiveWorkPending is already set, after
     * which the same advance clears that edge-triggered latch.  With no
     * unrelated Qt event, the controller would then wait forever for a frame
     * it never asked the host to schedule.
     *
     * Latch the render before invoking the callback so synchronous hosts also
     * observe a renderable request.  Do not move an existing gate forward;
     * every cut selected before the pending presentation belongs to the same
     * next frame. */
    if (!this->d->lodRefinementAwaitingFrame) {
	this->d->lodRefinementAwaitingFrame = TRUE;
	this->d->lodRefinementResumeAfterRenderSerial =
	    this->d->renderCompletionSerial + 1;
	if (this->d->lodRefinementResumeAfterRenderSerial == 0)
	    this->d->lodRefinementResumeAfterRenderSerial = 1;
    }
    const char *frameReason = reason ? reason : "lod-refinement-frame";
    /* Preserve a more specific render reason already installed by result
     * application side effects such as residency eviction. */
    if (!this->isRenderRequested())
	this->requestRender(frameReason);
    this->notifyFrameRequest(frameReason);
    this->d->reconcilePhase(BObolLodStateMachine::Event::WORK_SCHEDULED);
}

void
BObolViewController::scheduleResidentGrowthReallocationIfReady(void)
{
    if (!this->d->lodResidentGrowthPolicy.pending() ||
	!this->d->lodService)
	return;

    const bool streamIdle =
	this->d->lodService->activeTaskCountForGeneration(
	    this->d->lodActiveGeneration) == 0 &&
	this->d->lodService->queuedResultCountForGeneration(
	    this->d->lodActiveGeneration) == 0;
    /* Do not interrupt a structural/occurrence census or consume this edge
     * behind a presentation barrier.  The next progressive pump retries it
     * after the complete result wave and its coherent frame are visible. */
    const bool presentationReady =
	!this->d->lodSubmissionPending &&
	!this->d->lodRefinementAwaitingFrame &&
	!this->d->lodPublicationPolicy.pending() &&
	!this->d->lodBudgetPolicy.rescanAfterFrame();
    /* Pose-only motion owns a reversible renderer ceiling and deliberately
     * preserves retained occurrence cuts.  A zoom may spend newly available
     * residency during its bounded quality probes; otherwise retain this
     * obligation until the quiet transition. */
    const bool allocationAllowed =
	!this->d->lodInteractive ||
	this->d->lodViewDemandPolicy.interactionScaleChanged();
    if (getenv("BOBOL_LOD_TRACE_PASS")) {
	static std::atomic<unsigned int> residentGrowthTraceCount(0);
	const unsigned int traceIndex = residentGrowthTraceCount.fetch_add(1);
	if (traceIndex < 256)
	    bu_log("BObol LoD resident growth pending stream_idle=%d "
		   "scene_ready=%d automatic=%d allocation_allowed=%d "
		   "submission=%d coverage=%d refinement_frame=%d "
		   "publication=%d budget_frame=%d\n",
		   streamIdle ? 1 : 0,
		   (presentationReady &&
		    this->d->lodCoveragePolicy.effectiveComplete()) ? 1 : 0,
		   this->d->lodAutoSubmit ? 1 : 0,
		   allocationAllowed ? 1 : 0,
		   this->d->lodSubmissionPending ? 1 : 0,
		   this->d->lodCoveragePolicy.effectiveComplete() ? 1 : 0,
		   this->d->lodRefinementAwaitingFrame ? 1 : 0,
		   this->d->lodPublicationPolicy.pending() ? 1 : 0,
		   this->d->lodBudgetPolicy.rescanAfterFrame() ? 1 : 0);
    }

    /* A newly resident suffix can arrive just after the pass which proved
     * minimum-mesh coverage.  If that pass ended with one or more missing
     * prefixes, coverage is deliberately false—but resident growth is now
     * the event which can satisfy it.  Waiting for coverage before scheduling
     * another census leaves a closed cycle: no task, result, frame, or cursor
     * remains to change either latch.  Resume (or start) the bounded coverage
     * pass first; only its completed proof may consume the scene-wide
     * reallocation below. */
    if (this->d->lodAutoSubmit && streamIdle && presentationReady &&
	allocationAllowed &&
	!this->d->lodCoveragePolicy.effectiveComplete()) {
	if (!this->d->lodCoveragePolicy.active()) {
	    this->d->lodCoveragePolicy.activate(false);
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	}
	this->d->lodSubmissionPending = TRUE;
	this->d->lodSubmissionRescanPending = FALSE;
	this->markProgressiveWorkPending();
	this->notifyFrameRequest("lod-resident-growth-coverage");
	this->d->reconcilePhase(BObolLodStateMachine::Event::WORK_SCHEDULED);
	return;
    }
    if (!this->d->lodResidentGrowthPolicy.consumeIfReady(
	    this->d->lodAutoSubmit != FALSE,
	    streamIdle && presentationReady &&
		this->d->lodCoveragePolicy.effectiveComplete(),
	    allocationAllowed))
	return;

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
    this->d->lodBudgetPolicy.requestRetainedReallocation();
    /* Invalidate only the submission witness.  Camera, source, residency,
     * and calibrated capacity remain valid and are inputs to the new global
     * allocation. */
    this->d->lodLastSubmittedViewRevision.reset();
    this->d->lodLastSubmittedPolicyRevision.reset();
    this->d->resetLodConvergenceFraction();
    if (getenv("BOBOL_LOD_TRACE_PASS"))
	bu_log("BObol LoD resident growth scheduled scene-wide reallocation "
	       "view=%llu policy=%llu budget=%zu\n",
	       (unsigned long long)this->d->lodViewRevision.value(),
	       (unsigned long long)this->d->lodPolicyRevision.value(),
	       this->d->lodBudgetPolicy.currentBudget());
    this->markProgressiveWorkPending();
    this->notifyFrameRequest("lod-resident-growth-reallocation");
    this->d->reconcilePhase(BObolLodStateMachine::Event::WORK_SCHEDULED);
}

uint64_t
BObolViewController::getLastRenderTimeNanoseconds(void) const
{
    return this->d->lastRenderTimeNanoseconds;
}

uint64_t
BObolViewController::getSmoothedRenderTimeNanoseconds(void) const
{
    return this->d->smoothedRenderTimeNanoseconds;
}

uint64_t
BObolViewController::getLastBackgroundRenderTimeNanoseconds(void) const
{
    return this->d->lastBackgroundRenderTimeNanoseconds;
}

uint64_t
BObolViewController::getLastSceneRenderTimeNanoseconds(void) const
{
    return this->d->lastSceneRenderTimeNanoseconds;
}

uint64_t
BObolViewController::getLastProgressiveAdvanceTimeNanoseconds(void) const
{
    return this->d->lastProgressiveAdvanceTimeNanoseconds;
}

uint64_t
BObolViewController::getLastLodResultProcessingTimeNanoseconds(void) const
{
    return this->d->lastLodResultProcessingTimeNanoseconds;
}

uint64_t
BObolViewController::getLastProgressiveProviderTimeNanoseconds(void) const
{
    return this->d->lastProgressiveProviderTimeNanoseconds;
}

uint64_t
BObolViewController::getLastLodSubmissionTimeNanoseconds(void) const
{
    return this->d->lastLodSubmissionTimeNanoseconds;
}

uint64_t
BObolViewController::getLastPresentationSyncTimeNanoseconds(void) const
{
    return this->d->lastPresentationSyncTimeNanoseconds;
}

void
BObolViewController::noteFramePresented(void)
{
    const uint64_t now = this->beginRenderTiming();
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    if (this->d->presentedFrameSerial != UINT64_MAX)
	this->d->presentedFrameSerial++;
    uint64_t generalDisplaySample = 0;
    const auto updateDisplayedCadence = [this](uint64_t interval) {
	uint64_t &displayed =
	    this->d->displayedPresentationIntervalNanoseconds;
	if (!displayed) {
	    displayed = interval;
	    return;
	}

	/* The control EMAs are deliberately quick: they let the LoD policy
	 * observe a change in renderer capacity within a few frames.  That same
	 * response makes a faceplate sampled four times per second look
	 * needlessly nervous.  Give displayed cadence a 750 ms elapsed-time
	 * constant.  Computing alpha from elapsed time, rather than using a
	 * fixed sample weight, gives the same visual response at 10, 60, or
	 * 240 FPS. */
	constexpr double timeConstantNanoseconds = 750000000.0;
	const double alpha = -std::expm1(
	    -static_cast<double>(interval) / timeConstantNanoseconds);
	const double next = static_cast<double>(displayed) +
	    alpha * (static_cast<double>(interval) -
		static_cast<double>(displayed));
	displayed = static_cast<uint64_t>(
	    std::max(1.0, std::floor(next + 0.5)));
    };

    if (this->d->lastPresentationTimestampNanoseconds &&
	now > this->d->lastPresentationTimestampNanoseconds) {
	const uint64_t interval =
	    now - this->d->lastPresentationTimestampNanoseconds;
	/* A retained view is event driven.  The first frame after an idle
	 * second is not a one-FPS rendering sample; it is a discontinuity
	 * between two bursts.  Excluding that gap keeps the faceplate useful
	 * while zooming/rotating without reintroducing an idle repaint timer. */
	if (interval > 1000 && interval <= 250000000ULL) {
	    this->d->smoothedPresentationIntervalNanoseconds =
		this->d->smoothedPresentationIntervalNanoseconds ?
		(this->d->smoothedPresentationIntervalNanoseconds * 4u +
		 interval) / 5u : interval;
	    if (!this->d->lodGestureActive)
		generalDisplaySample = interval;
	} else if (interval > 250000000ULL &&
		   !this->d->lodGestureActive) {
	    /* The next burst must establish its own display baseline.  Retaining
	     * the prior value here would make a newly exposed window or a gesture
	     * after a long idle period report stale historical cadence. */
	    this->d->displayedPresentationIntervalNanoseconds = 0;
	}
    }
    this->d->lastPresentationTimestampNanoseconds = now;

    /* Only an active input gesture supplies a continuous-demand FPS sample.
     * Event-driven stable views can legitimately go 100 ms between a timer,
     * screenshot, or overlay repaint; treating those idle gaps as renderer
     * capacity made the scene budget depend on test/poll cadence. */
    if (this->d->lodGestureActive) {
	if (this->d->lastInteractivePresentationTimestampNanoseconds &&
	    now > this->d->lastInteractivePresentationTimestampNanoseconds) {
	    const uint64_t interval =
		now - this->d->
		    lastInteractivePresentationTimestampNanoseconds;
	    if (interval > 1000 && interval <= 250000000ULL) {
		this->d->smoothedInteractivePresentationIntervalNanoseconds =
		    this->d->
			smoothedInteractivePresentationIntervalNanoseconds ?
		    (this->d->
			smoothedInteractivePresentationIntervalNanoseconds * 4u +
		     interval) / 5u : interval;
		/* Navigation brackets continuous demand explicitly.  Prefer that
		 * clean burst sample over the interval from an unrelated
		 * event-driven frame immediately before the gesture. */
		updateDisplayedCadence(interval);
	    }
	}
	this->d->lastInteractivePresentationTimestampNanoseconds = now;
    } else {
	this->d->lastInteractivePresentationTimestampNanoseconds = 0;
	if (generalDisplaySample)
	    updateDisplayedCadence(generalDisplaySample);
    }
}

uint64_t
BObolViewController::getPresentedFrameSerial(void) const
{
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    return this->d->presentedFrameSerial;
}

uint64_t
BObolViewController::getRenderRequestSerial(void) const
{
    return this->renderRequestSerialGet();
}

uint64_t
BObolViewController::getRenderCompletionSerial(void) const
{
    return this->d->renderCompletionSerial;
}

uint64_t
BObolViewController::getLodSettleAfterRenderSerial(void) const
{
    return this->d->lodSettleAfterRenderSerial;
}

uint64_t
BObolViewController::getLodRefinementResumeAfterRenderSerial(void) const
{
    return this->d->lodRefinementAwaitingFrame ?
	this->d->lodRefinementResumeAfterRenderSerial : 0;
}

uint64_t
BObolViewController::getSmoothedPresentationIntervalNanoseconds(void) const
{
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    return this->d->smoothedPresentationIntervalNanoseconds;
}

uint64_t
BObolViewController::
getSmoothedInteractivePresentationIntervalNanoseconds(void) const
{
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    return this->d->smoothedInteractivePresentationIntervalNanoseconds;
}

uint64_t
BObolViewController::getDisplayedPresentationIntervalNanoseconds(void) const
{
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    return this->d->displayedPresentationIntervalNanoseconds;
}

int
BObolViewController::renderToImage(unsigned char **image,
				     int flip,
				     int alpha,
				     const SbColor *background,
				     SoDB::ContextManager *contextManager,
				     BObolProgressiveStatus *progressiveStatus)
{
    if (!image) {
	return BRLCAD_ERROR;
    }
    *image = NULL;

    if (!this->d->activeCamera || !this->getViewport() ||
	    !this->getRenderRoot()) {
	return BRLCAD_ERROR;
    }

    this->synchronizePresentation();
    /* LoD off -> force-realize the whole scene; LoD on -> coarse-first, letting
     * advanceProgressiveWork stream geometry in (matches the hardware and
     * headless render paths).  Default (no attachment / LoD off) is unchanged. */
    if (this->isForceRealizeDisplay())
	(void)this->realizePending();
    BObolProgressiveStatus localProgressiveStatus;
    (void)this->advanceProgressiveWork(NULL, &localProgressiveStatus);
    this->synchronizePresentation();
    if (progressiveStatus) {
	*progressiveStatus = localProgressiveStatus;
    }

    const SbViewportRegion &region = this->getViewportRegion();
    SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0) {
	return BRLCAD_ERROR;
    }

    SoDB::ContextManager *resolvedManager = contextManager ?
	contextManager : this->getRenderContextManager();
    if (!resolvedManager)
	return BRLCAD_ERROR;

    const bool cacheRenderer =
	resolvedManager == this->getRenderContextManager();
    std::unique_ptr<SoOffscreenRenderer> overrideRenderer;
    SoOffscreenRenderer *renderer = NULL;
    bool configureRenderer = false;
    if (!cacheRenderer) {
	overrideRenderer.reset(new SoOffscreenRenderer(resolvedManager, region));
	renderer = overrideRenderer.get();
	configureRenderer = true;
    } else if (!this->d->imageRenderer ||
	this->d->imageRendererManager != resolvedManager) {
	delete this->d->imageRenderer;
	this->d->imageRenderer = new SoOffscreenRenderer(resolvedManager, region);
	this->d->imageRendererManager = resolvedManager;
	renderer = this->d->imageRenderer;
	configureRenderer = true;
    } else {
	this->d->imageRenderer->setViewportRegion(region);
	renderer = this->d->imageRenderer;
    }
    if (configureRenderer) {
	renderer->getGLRenderAction()->setTransparencyType(
	    this->d->transparencyEnabled ? SoGLRenderAction::BLEND :
	    SoGLRenderAction::NONE);
	renderer->getGLRenderAction()->setSmoothing(
	    this->d->antialiasingEnabled);
	renderer->getGLRenderAction()->setNumPasses(1);
    }
    const SbColor imageBottom = background ? *background :
	this->d->backgroundBottom;
    const SbColor imageTop =
	(background && *background != this->d->backgroundBottom) ?
	*background : this->d->backgroundTop;
    const SbBool gradient = imageBottom != imageTop;
    renderer->setComponents(SoOffscreenRenderer::RGB);
    renderer->setBackgroundColor(imageBottom);
    if (gradient)
	renderer->setBackgroundGradient(imageBottom, imageTop);

    const uint64_t requestSerial = this->renderRequestSerialGet();
    const uint64_t started = this->beginRenderTiming();
    const SbBool rendered = renderer->render(this->getRenderRoot());
    this->completeRenderTiming(started);
    if (!rendered) {
	return BRLCAD_ERROR;
    }

    const unsigned char *buffer = renderer->getBuffer();
    if (!buffer) {
	return BRLCAD_ERROR;
    }

    const size_t width = static_cast<size_t>(size[0]);
    const size_t height = static_cast<size_t>(size[1]);
    const size_t srcStride = width * 3;
    const size_t dstBpp = alpha ? 4 : 3;
    const size_t dstStride = width * dstBpp;
    unsigned char *out = static_cast<unsigned char *>(bu_calloc(
	    height * dstStride, sizeof(unsigned char),
	    "bobol viewport image"));

    for (size_t y = 0; y < height; y++) {
	const size_t srcY = flip ? (height - y - 1) : y;
	const unsigned char *src = buffer + srcY * srcStride;
	unsigned char *dst = out + y * dstStride;
	if (alpha) {
	    for (size_t x = 0; x < width; x++) {
		dst[x * 4 + 0] = src[x * 3 + 0];
		dst[x * 4 + 1] = src[x * 3 + 1];
		dst[x * 4 + 2] = src[x * 3 + 2];
		dst[x * 4 + 3] = 255;
	    }
	} else {
	    std::memcpy(dst, src, srcStride);
	}
    }

    if (!localProgressiveStatus.hasMore)
	this->clearRenderRequestIfUnchanged(requestSerial);
    *image = out;
    return BRLCAD_OK;
}

SbBool
BObolViewController::isRenderRequested(void) const
{
    if (!this->d->endpointGraphicalRenderingEnabled.load(
	    std::memory_order_acquire))
	return FALSE;

    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    return this->d->renderRequested;
}

SbString
BObolViewController::getRenderReason(void) const
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    return this->d->renderReason;
}

static void
controller_accumulate_progressive_status(BObolProgressiveStatus &dst,
	const BObolProgressiveStatus &src)
{
    dst.providerCount += src.providerCount;
    dst.providerAdvanced += src.providerAdvanced;
    dst.lodResultsProcessed += src.lodResultsProcessed;
    dst.lodResultsApplied += src.lodResultsApplied;
    dst.submitted += src.submitted;
    dst.alreadyCached += src.alreadyCached;
    dst.expanded += src.expanded;
    dst.existing += src.existing;
    dst.remaining += src.remaining;
    dst.proxyPublished += src.proxyPublished;
    dst.metadataApplied += src.metadataApplied;
    dst.pendingTasks += src.pendingTasks;
    dst.inFlight += src.inFlight;
    dst.queuedResults += src.queuedResults;
    dst.queuedCacheWrites += src.queuedCacheWrites;
    if (src.changed)
	dst.changed = 1;
    if (src.hasMore)
	dst.hasMore = 1;
}

uint64_t
BObolViewController::registerProgressiveProvider(
    BObolProgressiveAdvanceCallback callback,
    void *userData,
    BObolProgressiveUserDataFreeCallback userDataFree)
{
    if (!callback)
	return 0;

    uint64_t token = this->d->progressiveProviderNextToken++;
    if (token == 0)
	token = this->d->progressiveProviderNextToken++;

    BObolProgressiveProviderRecord record;
    record.token = token;
    record.callback = callback;
    record.userData = userData;
    record.userDataFree = userDataFree;
    this->d->progressiveProviders.push_back(record);
    this->markProgressiveWorkPending();
    this->requestRender("progressive-provider");
    /* Registration is a distinct producer edge even when unrelated startup
     * work has already left both aggregate levels asserted.  Edge-coalesced
     * mark/request calls cannot wake a host in that state, so explicitly
     * announce the new provider once.  Its subsequent pumps remain governed
     * by the ordinary level-triggered pending contract. */
    this->notifyFrameRequest("progressive-provider-registered");
    return token;
}

void
BObolViewController::unregisterProgressiveProvider(uint64_t token)
{
    if (!token)
	return;

    for (std::vector<BObolProgressiveProviderRecord>::iterator it =
	 this->d->progressiveProviders.begin();
	 it != this->d->progressiveProviders.end(); ++it) {
	if (it->token != token)
	    continue;
	if (it->userDataFree)
	    (*it->userDataFree)(it->userData);
	this->d->progressiveProviders.erase(it);
	break;
    }
    if (this->d->progressiveProviders.empty())
	this->clearProgressiveWorkPending();
}

void
BObolViewController::clearProgressiveProviders(void)
{
    for (const BObolProgressiveProviderRecord &record :
	 this->d->progressiveProviders)
	if (record.userDataFree)
	    (*record.userDataFree)(record.userData);
    this->d->progressiveProviders.clear();
    this->clearProgressiveWorkPending();
}

void *
BObolViewController::findProgressiveProviderData(
    BObolProgressiveAdvanceCallback callback) const
{
    if (!callback)
	return NULL;
    for (const BObolProgressiveProviderRecord &record :
	 this->d->progressiveProviders)
	if (record.callback == callback)
	    return record.userData;
    return NULL;
}

uint64_t
BObolViewController::findProgressiveProviderToken(
    BObolProgressiveAdvanceCallback callback) const
{
    if (!callback)
	return 0;
    for (const BObolProgressiveProviderRecord &record :
	 this->d->progressiveProviders)
	if (record.callback == callback)
	    return record.token;
    return 0;
}

SbBool
BObolViewController::hasProgressiveProviders(void) const
{
    return this->d->progressiveProviders.empty() ? FALSE : TRUE;
}

void
BObolViewController::setDefaultProgressiveOptions(
    const BObolProgressiveOptions *options)
{
    if (options)
	this->d->defaultProgressiveOptions = *options;
    else
	this->d->defaultProgressiveOptions = BObolProgressiveOptions();
    this->markProgressiveWorkPending();
    this->requestRender("progressive-options");
}

const BObolProgressiveOptions &
BObolViewController::getDefaultProgressiveOptions(void) const
{
    return this->d->defaultProgressiveOptions;
}

int
BObolViewController::advanceProgressiveWork(
    const BObolProgressiveOptions *options,
    BObolProgressiveStatus *status)
{
    const uint64_t advanceStarted = this->beginRenderTiming();
    if (!options)
	options = &this->d->defaultProgressiveOptions;

    BObolProgressiveStatus localStatus;
    if (options->forceTerminalLodRefinement) {
	/* A deterministic/offline pump asks for the same view-dependent,
	 * pixel-exact terminal cut as an interactive host, but without waiting
	 * through responsiveness debounce tiers. */
	this->d->lodRefinementNotBeforeMicroseconds = 0;
	this->d->lodBudgetPolicy.resetPass();
    }
    this->d->lastLodResultProcessingTimeNanoseconds = 0;
    this->d->lastProgressiveProviderTimeNanoseconds = 0;
    this->d->lastLodSubmissionTimeNanoseconds = 0;

    /* Camera motion keeps a responsive aggregate cut.  A completed scale
     * frame may nevertheless spend one bounded quality step immediately;
     * this keeps continuous zoom progressive rather than magnifying one
     * coarse image until input stops.  Ordinary stable convergence still
     * begins after the 150 ms quiet debounce. */
    if (this->d->lodAutoSubmit && this->d->lodInteractive) {
	const int64_t now = bu_gettime();
	const int64_t elapsed = now - this->d->lodLastViewChangeMicroseconds;
	if (this->d->lodViewDemandPolicy.qualityProbePending()) {
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    const int activeMaximum = viewState ?
		viewState->maximumActiveProgressiveCut() : -1;
	    BObolLodViewDemandPolicy::QualityProbeInputs probeInputs;
	    probeInputs.activeMaximum = activeMaximum;
	    probeInputs.presentationCeiling =
		this->d->lodInteractiveProgressiveCeiling;
	    size_t presentedWork = 0;
	    if (viewState && viewState->lastCadPresentationFrameExact() &&
		viewState->lastCadPresentedRenderCost(presentedWork))
		probeInputs.presentedWork = presentedWork;
	    probeInputs.refinementFramePending =
		this->d->lodRefinementAwaitingFrame != FALSE;
	    probeInputs.publicationFramePending =
		this->d->lodPublicationPolicy.framePending();
	    probeInputs.motionFramePending =
		this->d->lodSettleAfterRenderSerial != 0;
	    const BObolLodViewDemandPolicy::QualityProbeDecision probe =
		this->d->lodViewDemandPolicy.beginQualityProbe(probeInputs);
	    if (probe.begin) {
		/* Admit no more than the next PoP population to presentation.
		 * Residency may already be much richer, but this probe must remain
		 * interruptible and must not jump from a cheap motion frame into the
		 * complete pixel-demanded prefix. */
		this->d->lodInteractiveProgressiveCeiling =
		    probe.progressiveCeiling;
		viewState->setCadPresentationProgressiveCutCeiling(
		    this->d->lodInteractiveProgressiveCeiling);
		this->d->lodBudgetPolicy.resetPass();
		this->advanceLodPolicyRevision(TRUE);
		this->d->lodDiscretePopulationTrialAvailable = TRUE;
		this->markProgressiveWorkPending();
		this->requestRender("lod-scale-interaction-refine");
	    }
	}
	if (!this->d->lodGestureActive &&
	    this->d->lodSettleAfterRenderSerial == 0 &&
	    this->d->lodLastViewChangeMicroseconds > 0 &&
	    elapsed >= 150000) {
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    const bool returnedToStartingScale =
		this->d->lodViewDemandPolicy.interactionScaleChanged() &&
		this->d->lodInteractionStartScaleSignatureValid &&
		this->d->lodInteractionStartScaleSignature.same(
		    this->d->lodViewScaleSignature);
	    const bool scaleChanged =
		this->d->lodViewDemandPolicy.interactionScaleChanged() &&
		!returnedToStartingScale;
	    const bool haveRetainedMeshPayloads =
		this->getActiveLodMeshPayloadCount() > 0;
	    BObolLodPresentationPolicy::RestoreInputs restoreInputs;
	    restoreInputs.orthographic = this->d->activeCamera &&
		this->d->activeCamera->isOfType(
		    SoOrthographicCamera::getClassTypeId());
	    restoreInputs.scaleChanged = scaleChanged;
	    restoreInputs.retainedMeshPayloads = haveRetainedMeshPayloads;
	    restoreInputs.viewEpoch = this->d->lodViewRevision;
	    restoreInputs.population =
		controller_lod_presentation_population(viewState);
	    restoreInputs.currentProgressiveCeiling =
		this->d->lodInteractiveProgressiveCeiling;
	    restoreInputs.currentPointProxyPixelThreshold =
		this->d->lodPresentationPointProxyPixelThreshold;
	    /* A pose-only quiet pass still rechecks current visibility and
	     * resource pressure, but a fully covered orthographic population has
	     * no depth-dependent LoD demand.  Preserve its occurrence cuts unless
	     * measured FPS/resource admission explicitly requires coarsening. */
	    this->d->lodRetainPoseOccurrenceCuts =
		restoreInputs.orthographic && !scaleChanged &&
		this->d->lodInteractionStartedFromReadyView &&
		haveRetainedMeshPayloads ? TRUE : FALSE;
	    /* Orthographic depth is irrelevant, but rotation still changes
	     * projected area and silhouette.  Perspective pose and scale changes
	     * can move importance as well.  Finish one exact demand census before
	     * authorizing a single redistribution of the retained scene budget. */
	    this->d->lodRetainedImportanceCensusPending =
		this->d->lodInteractionStartedFromReadyView &&
		haveRetainedMeshPayloads ? TRUE : FALSE;
	    const BObolLodPresentationPolicy::RestoreDecision restore =
		this->d->lodPresentationPolicy.beginQuiet(restoreInputs);
	    if (restore.apply) {
		this->d->lodInteractiveProgressiveCeiling =
		    restore.progressiveCeiling;
		this->d->lodPresentationPointProxyPixelThreshold =
		    restore.pointProxyPixelThreshold;
		if (viewState) {
		    viewState->setCadPresentationProgressiveCutCeiling(
			this->d->lodInteractiveProgressiveCeiling);
		    viewState->setCadPresentationPointProxyPixelThreshold(
			this->d->lodPresentationPointProxyPixelThreshold);
		}
	    }
	    this->d->lodInteractive = FALSE;
	    const SbBool completedScaleInteraction =
		this->d->lodViewDemandPolicy.finishQuietInteraction(
		    returnedToStartingScale) ?
		    TRUE : FALSE;
	    this->d->lodInteractionStartScaleSignatureValid = FALSE;
	    this->d->lodInteractionStartedFromReadyView = FALSE;
	    /* Motion calibration targets 60 FPS and may intentionally lower the
	     * aggregate allowance behind a renderer-only ceiling.  Stable drawing
	     * targets 20 FPS; seed its handoff from the last proven stable scene
	     * budget, not from the transient motion floor.  Genuine quiet-frame
	     * overload evidence may still reduce this restored value. */
	    if (this->d->lodStableBudgetBeforeInteractionValid)
		this->d->lodBudgetPolicy.raiseCurrentBudget(
		    this->d->lodStableBudgetBeforeInteraction);
	    this->d->lodStableBudgetBeforeInteraction = 0;
	    this->d->lodStableBudgetBeforeInteractionValid = FALSE;
	    this->d->lodDiscretePopulationTrialAvailable = FALSE;
	    this->d->lodReleaseCutFloorActive = FALSE;
	    /* A motion ceiling also applies to native progressive wire stored in
	     * the standing CAD assembly.  It has no view-state mesh payload and
	     * therefore no occurrence pass capable of completing the mesh
	     * handoff below.  Restore its stable range directly; the policy
	     * revision already requests the one frame needed to present it. */
	    this->d->lodBudgetPolicy.resetProbeSeries();
	    this->d->lodBudgetPolicy.resetOverloadRecovery();
	    this->d->lodTargetPixelError = 1.0f;
	    this->d->lodStablePointProxyCalibrationPending = FALSE;
	    this->d->lodPointProxyTriangleRecoveryPending = FALSE;
	    if (viewState) {
		viewState->setCadPresentationCameraMotionFrameReuse(FALSE);
		/*
		 * activeFaceCount() is the retained occurrence population, not
		 * necessarily the cut just presented under the motion ceiling.
		 * Promoting the stable budget to that value authorized the hidden
		 * rich population before the first quiet capacity check.  Keep the
		 * measured presentation and reconcile the retained cuts first.
		 */
	    }
	    /* Moving from the responsive pixel target to the stable target is a
	     * continuation of zoom demand, not a new cold-coverage epoch. */
	    this->advanceLodPolicyRevision(completedScaleInteraction);
	} else {
	    /* Keep the host frame pump alive through the idle debounce. */
	    localStatus.hasMore = 1;
	}
    }

    const uint64_t residentConsumerId =
	static_cast<uint64_t>(reinterpret_cast<uintptr_t>(this));
    const size_t queuedCompactionResults = this->d->lodService ?
	this->d->lodService->
	    queuedResidentMeshCompactionResultCountForDiagnostics(
		residentConsumerId) : 0;
    if (queuedCompactionResults) {
	if (this->d->lodInteractive) {
	    /* Keep the last complete rich generation drawable while the user is
	     * moving.  The bounded completion queue back-pressures compaction
	     * workers until the quiet-view publication phase resumes. */
	    localStatus.hasMore = 1;
	} else {
	    const size_t applyLimit = options->maxLodResults ?
		std::min<size_t>(options->maxLodResults, 1024) : 1024;
	    std::vector<BObolLodResidentCompaction> compacted;
	    this->d->lodService->drainResidentMeshCompactions(
		residentConsumerId, compacted, applyLimit);
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    size_t published = 0;
	    if (viewState) {
		for (const BObolLodResidentCompaction &result : compacted)
		    published +=
			viewState->applyResidentMeshCompaction(result);
	    }
	    localStatus.lodResultsProcessed += compacted.size();
	    localStatus.lodResultsApplied += published;
	    if (published)
		localStatus.changed = 1;
	    if (this->d->lodService->
		    queuedResidentMeshCompactionResultCountForDiagnostics(
			residentConsumerId) > 0)
		localStatus.hasMore = 1;
	}
    }

    SbBool nonLodPresentationChanged = localStatus.changed ? TRUE : FALSE;
    const size_t queuedLodResults = this->d->lodService ?
	this->d->lodService->queuedResultCountForGeneration(
	    this->d->lodActiveGeneration) : 0;
    const SbBool havePendingLodResults =
	this->hasPendingLodResults() || queuedLodResults > 0;
    SbBool coalesceLodResults = FALSE;
    const SbBool holdRicherResultsDuringInteraction =
	havePendingLodResults && this->d->lodInteractive &&
	!this->d->lodViewDemandPolicy.interactionScaleChanged() &&
	this->d->lodCoveragePolicy.coverageComplete() &&
	this->getActiveLodMeshPayloadCount() > 0;
    if (havePendingLodResults && this->d->lodService) {
	/* The first mesh is the time-to-useful-visual boundary and is never
	 * delayed.  Once a scene already has mesh data, coalesce a continuing
	 * many-asset stream until a full apply wave is ready or its latency
	 * bound expires.  A final/isolated result (no outstanding worker work)
	 * publishes immediately, preserving single-huge-mesh and sparse-scene
	 * behavior.
	 *
	 * Stable bulk population uses a 250 ms ceiling: four visibly progressive
	 * updates per second without paying hundreds of increasingly expensive
	 * whole-scene update/render traversals.  Motion uses 100 ms so newly
	 * visible occurrences catch up promptly while the interactive face
	 * budget protects FPS. */
	const size_t applyWave = options->maxLodResults ?
	    options->maxLodResults : 256;
	const SbBool continuingWorkerStream =
	    this->d->lodService->activeTaskCountForGeneration(
		this->d->lodActiveGeneration) > 0;
	const SbBool sceneHasMesh =
	    this->getActiveLodMeshPayloadCount() > 0;
	int64_t firstReady =
	    this->d->lodResultsFirstReadyMicroseconds.load();
	const int64_t now = bu_gettime();
	if (firstReady <= 0) {
	    int64_t expected = 0;
	    (void)this->d->lodResultsFirstReadyMicroseconds.
		compare_exchange_strong(expected, now);
	    firstReady =
		this->d->lodResultsFirstReadyMicroseconds.load();
	}
	const int64_t latencyLimit =
	    this->d->lodInteractive ? 100000 : 250000;
	const int64_t resultAge =
	    now >= firstReady ? now - firstReady : latencyLimit;
	coalesceLodResults =
	    sceneHasMesh && continuingWorkerStream &&
	    queuedLodResults < applyWave && resultAge < latencyLimit;
    }
    if (holdRicherResultsDuringInteraction) {
	/*
	 * Worker-produced arrays are immutable and the result queue is bounded
	 * and coalescing (currently 2048 slots).  Keeping richer generations
	 * there preserves the renderer's last prepared generation during
	 * pose-only camera input.  Scale-changing input is different: zoom must
	 * publish newly loaded suffixes as bounded result waves so detail improves
	 * under the pointer rather than only after the quiet debounce.  Coverage
	 * is reset on every view revision, so a newly exposed box-only leaf is
	 * never held behind this optimization.
	 */
	localStatus.hasMore = 1;
    } else if (coalesceLodResults) {
	localStatus.hasMore = 1;
    } else if (havePendingLodResults) {
	const uint64_t resultStarted = this->beginRenderTiming();
	(void)this->processPendingLodResults(options->maxLodResults,
	    options->maxLodApplyMicroseconds);
	const uint64_t resultCompleted = this->beginRenderTiming();
	this->d->lastLodResultProcessingTimeNanoseconds =
	    resultCompleted > resultStarted ? resultCompleted - resultStarted : 0;
	localStatus.lodResultsProcessed += this->d->lastLodResultCount;
	localStatus.lodResultsApplied += this->d->lastLodAppliedResultCount;
	if (this->d->lastLodAppliedResultCount > 0)
	    localStatus.changed = 1;
	if (!this->d->lodService ||
	    this->d->lodService->queuedResultCountForGeneration(
		this->d->lodActiveGeneration) == 0) {
	    this->d->lodResultsFirstReadyMicroseconds.store(0);
	    /* A result can race the empty observation without producing another
	     * empty->non-empty callback.  Recheck after clearing and install a
	     * fresh age origin if necessary. */
	    if (this->d->lodService &&
		this->d->lodService->queuedResultCountForGeneration(
		    this->d->lodActiveGeneration) > 0) {
		this->d->lodResultsPending.store(1);
		int64_t expected = 0;
		(void)this->d->lodResultsFirstReadyMicroseconds.
		    compare_exchange_strong(expected, bu_gettime());
	    }
	}
    }

    size_t providerLimit = options->maxProviders;
    size_t providerIndex = 0;
    const uint64_t providerStarted = this->beginRenderTiming();
    for (const BObolProgressiveProviderRecord &record :
	 this->d->progressiveProviders) {
	if (!record.callback)
	    continue;
	if (providerLimit && providerIndex >= providerLimit) {
	    localStatus.hasMore = 1;
	    break;
	}

	BObolProgressiveStatus providerStatus;
	providerStatus.providerCount = 1;
	int providerRet = (*record.callback)(this, record.userData, options,
			 &providerStatus);
	if (providerRet > 0) {
	    providerStatus.providerAdvanced++;
	} else if (providerRet < 0) {
	    providerStatus.hasMore = 0;
	}
	controller_accumulate_progressive_status(localStatus,
		providerStatus);
	if (providerStatus.changed)
	    nonLodPresentationChanged = TRUE;
	providerIndex++;
    }
    const uint64_t providerCompleted = this->beginRenderTiming();
    this->d->lastProgressiveProviderTimeNanoseconds =
	providerCompleted > providerStarted ?
	providerCompleted - providerStarted : 0;

    /* Providers publish the latest compact occurrences and their source-mesh
     * requests during this pump.  Submit view demand afterward so a provider
     * that completes in one pass cannot leave a newly published model at its
     * boxes forever simply because there is no reason to schedule a second
     * pass. */
    struct bv_lod_policy lodPolicy;
    bv_lod_policy_init(&lodPolicy);
    this->d->viewAttachment->getLodPolicy(&lodPolicy);
    this->d->lodCoveragePolicy.setRequired(
	this->d->lodAutoSubmit && this->d->lodService &&
	lodPolicy.policy != BV_LOD_OFF && lodPolicy.mesh_enabled);
    if (this->d->lodService) {
	const uint64_t admissionRevision =
	    this->d->lodService->residentMeshAdmissionRevision();
	if (this->d->lodResidentAdmissionRevision != 0 &&
	    admissionRevision !=
		this->d->lodResidentAdmissionRevision) {
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    if (viewState &&
		viewState->memoryLimitedPayloadCount() > 0) {
		/* Reclamation made stable bytes available.  Resume only the
		 * sparse unsatisfied frontier; the fixed camera and every
		 * existing mesh presentation remain untouched meanwhile. */
		this->d->lodSubmissionSourceIndex = 0;
		this->d->lodSubmissionEntryOffset = 0;
		this->d->clearLodSubmissionPlan();
		this->d->lodSubmissionPending = TRUE;
		this->markProgressiveWorkPending();
	    }
	}
	this->d->lodResidentAdmissionRevision =
	    admissionRevision;
    }
    this->scheduleResidentGrowthReallocationIfReady();
    const int64_t refinementNow = bu_gettime();
    const bool refinementCooling =
	!options->forceTerminalLodRefinement &&
	this->d->lodRefinementNotBeforeMicroseconds > refinementNow &&
	!this->d->lodInteractive && !localStatus.changed;
    if (refinementCooling)
	localStatus.hasMore = 1;
    else if (this->d->lodRefinementNotBeforeMicroseconds > 0 &&
	     refinementNow >=
		this->d->lodRefinementNotBeforeMicroseconds)
	this->d->lodRefinementNotBeforeMicroseconds = 0;

    if (this->d->lodAutoSubmit && lodPolicy.policy != BV_LOD_OFF &&
	lodPolicy.mesh_enabled && !refinementCooling) {
	this->d->forceTerminalLodRefinement =
	    options->forceTerminalLodRefinement;
	const uint64_t submissionStarted = this->beginRenderTiming();
	(void)this->submitLodRequestsIfNeeded();
	const uint64_t submissionCompleted = this->beginRenderTiming();
	this->d->lastLodSubmissionTimeNanoseconds =
	    submissionCompleted > submissionStarted ?
	    submissionCompleted - submissionStarted : 0;
	/* A bounded compact-entry scan can have no worker tasks or results at
	 * the boundary between two chunks.  The submission cursor itself is
	 * pending work and must keep the frame pump alive.  Relying on
	 * submitLodRequests() calling markProgressiveWorkPending() is not
	 * sufficient: the status-based epilogue below otherwise clears that
	 * latch in the same advance.  The resulting scene refines only when an
	 * unrelated paint, checkpoint, or input event happens to arrive. */
	if (this->d->lodSubmissionPending)
	    localStatus.hasMore = 1;
	/* A budget-limited quiet cut may intentionally request one unchanged
	 * frame to improve its throughput estimate.  Keep the host pump alive
	 * until completeRenderTiming() turns that presented probe into a new
	 * admission pass. */
	if (this->d->lodBudgetPolicy.rescanAfterFrame())
	    localStatus.hasMore = 1;
	this->d->forceTerminalLodRefinement = FALSE;
    }

    /* Provider status describes database streaming.  Mesh refinement runs on
     * the controller-owned service and is independent of any one provider, so
     * sample it here.  Previously these fields remained zero and the common
     * host pump could declare a frame stable while a PoP task was still
     * queued or running. */
    if (this->d->lodService) {
	localStatus.pendingTasks +=
	    this->d->lodService->pendingTaskCountForGeneration(
		this->d->lodActiveGeneration);
	localStatus.pendingTasks +=
	    this->d->lodService->delayedTaskCountForGeneration(
		this->d->lodActiveGeneration);
	localStatus.inFlight +=
	    this->d->lodService->executingTaskCountForGeneration(
		this->d->lodActiveGeneration);
	localStatus.queuedResults +=
	    this->d->lodService->queuedResultCountForGeneration(
		this->d->lodActiveGeneration);
	localStatus.queuedCacheWrites +=
	    this->d->lodService->queuedCacheWriteCountForGeneration(
		this->d->lodActiveGeneration);
	localStatus.pendingTasks +=
	    this->d->lodService->
		pendingResidentMeshCompactionCountForDiagnostics();
	localStatus.queuedResults +=
	    this->d->lodService->
		queuedResidentMeshCompactionResultCountForDiagnostics(
		    residentConsumerId);
    }

    const int pending_service_work =
	(localStatus.pendingTasks > 0 || localStatus.inFlight > 0 ||
	 localStatus.queuedResults > 0 || localStatus.queuedCacheWrites > 0) ?
	1 : 0;
    if (pending_service_work)
	localStatus.hasMore = 1;
    const int pending_lod_realization_work =
	this->d->lodService &&
	(this->d->lodService->activeTaskCountForGeneration(
	     this->d->lodActiveGeneration) > 0 ||
	 this->d->lodService->queuedResultCountForGeneration(
	     this->d->lodActiveGeneration) > 0 ||
	 this->d->lodService->queuedCacheWriteCountForGeneration(
	     this->d->lodActiveGeneration) > 0) ?
	1 : 0;

    /* Refinement and reclamation are separate phases.  A quiet view first
     * reaches its fast 1 px display target and may then enter the bounded
     * subpixel tiers when exact timing and memory witnesses prove sufficient
     * headroom.  Only after that terminal quality decision and a longer quiet
     * interval do we replace this consumer's complete demand snapshot and
     * trim shared CPU prefixes to the aggregate maximum. */
    const int64_t compactionNow = bu_gettime();
    BObolLodCompactionPolicy::Inputs compactionInputs;
    compactionInputs.automatic = this->d->lodAutoSubmit != FALSE;
    compactionInputs.interactive = this->d->lodInteractive != FALSE;
    compactionInputs.coverageRequired =
	lodPolicy.policy != BV_LOD_OFF && lodPolicy.mesh_enabled;
    compactionInputs.coverageComplete =
	this->d->lodCoveragePolicy.coverageComplete();
    compactionInputs.coverageProgressPending = localStatus.hasMore != 0;
    compactionInputs.settlingPending =
	compactionInputs.coverageRequired &&
	(this->d->lodCoveragePolicy.active() ||
	 this->d->lodSubmissionRescanPending ||
	 this->d->lodRetainedRefinementPending ||
	 this->d->lodRetainedResidencyPending ||
	 this->d->lodBudgetPolicy.rescanAfterFrame() ||
	 this->d->lodPresentationPolicy.handoffPending() ||
	 this->d->lodRefinementAwaitingFrame ||
	 this->d->lodPublicationPolicy.pending() ||
	 this->d->lodBudgetPolicy.calibrationFramesRemaining() != 0 ||
	 this->d->lodStablePointProxyCalibrationPending ||
	 this->d->lodPointProxyTriangleRecoveryPending ||
	 this->d->lodResidentGrowthPolicy.pending() ||
	 this->d->lodHeadroomPolicy.retryPending());
    compactionInputs.nowMicroseconds = compactionNow;
    compactionInputs.realizationPending =
	pending_lod_realization_work != 0;
    compactionInputs.submissionPending =
	this->d->lodSubmissionPending != FALSE;
    compactionInputs.serviceAvailable = this->d->lodService != NULL;
    const BObolLodCompactionPolicy::Decision compactionDecision =
	this->d->lodCompactionPolicy.decide(compactionInputs);
    if (compactionDecision.retiredRequest &&
	this->d->lodResourcePolicy.recoveryPending()) {
	/* The coverage producer disappeared (source erase, cancellation, LoD
	 * disable, or empty view).  Retiring compaction is the complete response
	 * to this pressure revision; leaving the resource edge unacknowledged
	 * otherwise advertises background work with no task, cursor, timer, or
	 * frame capable of discharging it. */
	this->d->lodResourcePolicy.markRecoveryHandled();
    }
    if (compactionDecision.keepPumpAlive)
	localStatus.hasMore = 1;
    if (compactionDecision.admission ==
	BObolLodCompactionPolicy::Admission::PLAN) {
	    std::vector<BObolLodResidentDemand> demands;
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    const uint64_t demandRevision = viewState ?
		viewState->residentMeshDemandRevision() :
		this->d->lodViewRevision.value();
	    const SbBool continuePlanning =
		this->d->lodCompactionPolicy.continues(demandRevision) ?
		    TRUE : FALSE;
	    if (viewState && !continuePlanning)
		viewState->residentMeshDemands(demands);
	    SbBool planningComplete = FALSE;
	    const size_t queued =
		this->d->lodService->scheduleResidentMeshCompaction(
		    residentConsumerId,
		    demandRevision, demands,
		    &planningComplete);
	    this->d->lodCompactionPolicy.finishPlanning(
		planningComplete != FALSE, demandRevision, bu_gettime());
	    if (planningComplete &&
		this->d->lodResourcePolicy.recoveryPending()) {
		/* One complete demand snapshot has now been admitted to the
		 * worker service for this pressure edge.  Persistent pressure can
		 * legitimately describe a visible working set larger than the
		 * configured limit; do not reschedule the identical compaction on
		 * every frame.  A falling/rising completed-frame edge starts a new
		 * recovery revision. */
		this->d->lodResourcePolicy.markRecoveryHandled();
	    }
	    if (queued || !planningComplete)
		localStatus.hasMore = 1;
    }

    BObolLodPublicationPolicy::Inputs publicationInputs;
    publicationInputs.nowMicroseconds = bu_gettime();
    publicationInputs.observedRenderNanoseconds = std::max(
	this->d->lastSceneRenderTimeNanoseconds,
	this->d->lastRenderTimeNanoseconds);
    publicationInputs.interactive = this->d->lodInteractive != FALSE;
    const BObolLodPublicationPolicy::Decision publicationDecision =
	this->d->lodPublicationPolicy.decide(publicationInputs);
    if (publicationDecision.keepPumpAlive)
	localStatus.hasMore = 1;
    if (publicationDecision.requestFrame) {
	this->requestRender("lod-result-batch");
	if (this->d->lodRefinementAwaitingFrame)
	    this->scheduleLodRefinementFrame("lod-result-batch");
    }

    /*
     * A refinement/calibration barrier is, by definition, waiting for a
     * completed presentation.  Result publication and Qt paint coalescing
     * can consume the request which originally installed the barrier before
     * the barrier's target render serial is reached.  Never allow that state
     * to advertise an idle pump with no frame request: reasserting the edge
     * is idempotent, and requestRender() wakes the host only on false->true.
     */
    if (this->hasPendingLodRefinementFrame()) {
	localStatus.hasMore = 1;
	/* A partial result installs its presentation barrier as soon as its
	 * immutable binding is applied, but sparse publication deliberately does
	 * not request that frame until the batch deadline.  Treating the barrier
	 * itself as a due frame collapses every batch to one result and restores
	 * the expensive whole-scene repaint stream this policy exists to avoid.
	 * The publication timer remains the liveness witness until decide()
	 * atomically replaces it with framePending(). */
	const SbBool publicationDeadlinePending =
	    this->d->lodRefinementAwaitingFrame &&
	    this->d->lodPublicationPolicy.awaitingDeadline();
	if (!publicationDeadlinePending && !this->isRenderRequested())
	    this->requestRender("lod-refinement-pending");
    }

    if (localStatus.hasMore)
	this->markProgressiveWorkPending();
    else
	this->clearProgressiveWorkPending();

    /* Pending background work needs a host timer, not a duplicate render of
     * unchanged pixels.  This distinction is especially important for
     * OSMesa, where merely repainting a multi-million-triangle stable cut can
     * consume seconds.  The Qt host pumps pending work independently and
     * actual result/cut changes install their own render requests. */
    if (localStatus.changed &&
	(nonLodPresentationChanged ||
	 this->d->lodPublicationPolicy.framePending()))
	this->requestRender("progressive-update");

    if (status)
	*status = localStatus;

    size_t phasePending = localStatus.pendingTasks;
    phasePending = localStatus.inFlight > SIZE_MAX - phasePending ?
	SIZE_MAX : phasePending + localStatus.inFlight;
    phasePending = localStatus.queuedResults > SIZE_MAX - phasePending ?
	SIZE_MAX : phasePending + localStatus.queuedResults;
    phasePending = localStatus.queuedCacheWrites > SIZE_MAX - phasePending ?
	SIZE_MAX : phasePending + localStatus.queuedCacheWrites;
    this->d->reconcilePhase(BObolLodStateMachine::Event::WORK_PUMPED,
	phasePending);

    const uint64_t advanceCompleted = this->beginRenderTiming();
    this->d->lastProgressiveAdvanceTimeNanoseconds =
	advanceCompleted > advanceStarted ? advanceCompleted - advanceStarted : 0;
    return (localStatus.changed || localStatus.hasMore) ? 1 : 0;
}

void
BObolViewController::markProgressiveWorkPending(void)
{
    if (this->d->progressiveWorkPending.exchange(1) == 0)
	this->notifyFrameRequest("progressive-work");
}

void
BObolViewController::clearProgressiveWorkPending(void)
{
    this->d->progressiveWorkPending.store(0);
}

SbBool
BObolViewController::hasProgressiveWorkPending(void) const
{
    return this->d->progressiveWorkPending.load() != 0 ? TRUE : FALSE;
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
    this->d->lodDiscretePopulationTrialAvailable = FALSE;
    this->d->lodStaticOverscanActive = FALSE;
    this->d->lodCoveragePolicy.reset();
    if (service)
	this->d->lodCoveragePolicy.activate(true);
    this->d->lodRetainedRefinementPending = FALSE;
    this->d->lodRetainedResidencyPending = FALSE;
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
    this->d->lodRefinementAwaitingFrame = FALSE;
    this->d->lodRefinementResumeAfterRenderSerial = 0;
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodPublicationPolicy.reset();

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
    this->d->lodRetainedRefinementPending = FALSE;
    this->d->lodRetainedResidencyPending = FALSE;
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
    this->d->lodRefinementAwaitingFrame = FALSE;
    this->d->lodRefinementResumeAfterRenderSerial = 0;
    this->d->lodPublicationPolicy.reset();
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodInteractiveProgressiveCeiling = -1;
    this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
    this->d->lodPointProxyCalibrationPolicy.reset();
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodPointProxyTriangleRecoveryPending = FALSE;
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
    this->cancelActiveLodGeneration();
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
    this->d->lodAutoSubmit = enabled ? TRUE : FALSE;
    this->d->lodStaticOverscanActive = FALSE;
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
	this->d->lodInteractiveProgressiveCeiling = -1;
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->lodPointProxyCalibrationPolicy.reset();
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	this->d->lodPointProxyTriangleRecoveryPending = FALSE;
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
    return this->d->lodRefinementAwaitingFrame ||
	this->d->lodBudgetPolicy.rescanAfterFrame() ||
	(this->d->lodPresentationPolicy.handoffPending() &&
	 !this->d->lodSubmissionPending) ||
	this->d->lodStablePointProxyCalibrationPending ||
	this->d->lodHeadroomPolicy.retryPending() ||
	this->d->lodResidentGrowthPolicy.pending();
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
	    this->d->lodRefinementAwaitingFrame)
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
	this->d->lodCoveragePolicy.deactivate();
	this->d->lodResidentGrowthPolicy.reset();
	this->d->lodViewDemandPolicy.clearDemandRefresh();
	this->d->lodCoveragePolicy.clearPassCounters();
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
    bool useSourceDelta = false;
    bool extendedPendingDelta = false;
    bool pendingDeltaNeedsFullRescan = false;
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
	    if (!changedSource->getDisplayMeshLodChangedEntries(
		    previousInventory, changedEntries) ||
		changedEntries.empty()) {
		pendingDeltaNeedsFullRescan = true;
		continue;
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
		if (changedSource->getDisplayMeshLodChangedEntries(
			previousInventory, changedEntries) &&
		    !changedEntries.empty()) {
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
    if (sourceSetChanged || inventoryChanged) {
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
	this->d->clearLodConvergenceCandidates();
	this->d->resetLodConvergenceFraction();
    }
    if (sourceSetChanged || useSourceDelta ||
	!this->d->lodSubmissionPending) {
	this->d->lodSubmissionSourceIndex =
	    useSourceDelta ? sourceDeltaFirst : 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	/*
	 * A selective inventory delta is a latency optimization, not a proof
	 * that the complete current view is covered.  In particular, a cold
	 * compact source may publish several large, overlapping deltas while
	 * its structural registry is still growing.  Declaring coverage from
	 * only the last delta made completion timing-dependent: untouched
	 * ranges remained structural boxes until the next camera change.
	 *
	 * Consume the exact delta first so newly arrived leaves become useful
	 * promptly, then perform one bounded all-source pass.  The completion
	 * path clears delta mode before that pass, and later inventory changes
	 * extend the in-flight delta without losing this obligation.
	 */
	this->d->lodSubmissionRescanPending =
	    useSourceDelta ? TRUE : FALSE;
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
    double lodAspect = controller_aspect_from_region(this->d->viewportRegion);
    if (lodAspect <= SMALL_FASTF)
	lodAspect = 1.0;
    const SbViewVolume lodViewVolume =
	this->d->activeCamera->getViewVolume(static_cast<float>(lodAspect));
    const SbMatrix lodViewProjection = lodViewVolume.getMatrix();
    if (!this->d->lodBudgetPolicy.passInitialized()) {
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
	controller_lod_effective_population_cost(sources, lodState,
	    this->d->lodPresentationPointProxyPixelThreshold,
	    &lodViewProjection, activeCost, minimumActiveCost);
	const SbBool scaleQualityProbe =
	    this->d->lodInteractive &&
	    this->d->lodViewDemandPolicy.qualityBudgetActive() ? TRUE : FALSE;
	/* A zoom quality frame deliberately has a lower cadence target than the
	 * ordinary motion frame: coherent PoP populations are discrete, and the
	 * next useful cut can be just beyond a strict 20/60 Hz allowance.  Use the
	 * same 10 Hz hard floor which decides whether that cut remains visible.
	 * Exact hierarchy costs still prevent an unexpectedly huge cut jump. */
	const float targetFps =
	    scaleQualityProbe ?
		BObolLodViewDemandPolicy::
		    qualityTargetFramesPerSecond() :
		(this->d->lodInteractive ?
		    this->d->lodInteractiveTargetFps :
		    this->d->quietAllocationTargetFps());
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
	/* A sparse unsatisfied-refinement plan is not a retained-recovery plan.
	 * Reusing it lets its first subset consume the complete upgrade allowance;
	 * the later all-occurrence pass then sees zero and normalizes everything
	 * else to minimum.  Give the mode transition an explicit plan epoch and
	 * restart at the first source. */
	if (budget.retainedAdmission &&
	    !this->d->lodSubmissionRetainedAdmissionMode) {
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRetainedAdmissionMode = TRUE;
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	} else if (!budget.retainedAdmission &&
	    this->d->lodSubmissionRetainedAdmissionMode) {
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRetainedAdmissionMode = FALSE;
	    this->d->lodRetainedAdmissionMaximumNormalizedError =
		std::numeric_limits<double>::infinity();
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

    float scenePixelError = this->d->lodTargetPixelError;
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

    if (this->d->lodBudgetPolicy.retainedAdmission() &&
	this->d->lodSubmissionSourceIndex == 0 &&
	this->d->lodSubmissionEntryOffset == 0) {
	size_t maximumProtectedBudget =
	    this->d->lodBudgetPolicy.currentBudget();
	if (!this->d->lodInteractive &&
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
	size_t protectedFloorBudget = 0;
	uint64_t protectedFloorSignature = 0;
	this->d->lodRetainedAdmissionMaximumNormalizedError =
	    controller_lod_retained_error_ceiling(
		sources, this->d->viewAttachment->getViewLodState(),
		this->d->lodBudgetPolicy.currentBudget(),
		!this->d->lodInteractive &&
		    !this->d->lodBudgetPolicy.retainedQualityFloorRejected(),
		maximumProtectedBudget, &protectedFloorBudget,
		&protectedFloorSignature,
		this->d->lodViewRevision.value(),
		this->d->lodPolicyRevision.value(),
		this->d->lodPresentationPointProxyPixelThreshold,
		&lodViewProjection);
	this->d->lodBudgetPolicy.setRetainedQualityFloorBudget(
	    protectedFloorBudget, protectedFloorSignature, sceneActiveCost,
	    sceneMinimumActiveCost);
	if (getenv("BOBOL_LOD_TRACE_BUDGET"))
	    bu_log("BObol LoD retained importance ceiling=%.6f "
		   "budget=%zu protected_limit=%zu\n",
		   this->d->lodRetainedAdmissionMaximumNormalizedError,
		   this->d->lodBudgetPolicy.currentBudget(),
		   maximumProtectedBudget);
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
	 * may still contain legacy/direct SoBRLMeshShape children whose request
	 * metadata is owned by the shape, so let the submit action inspect those
	 * sources.  This keeps compatibility at the representation boundary
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
	action.setViewVolume(&lodViewVolume, scenePixelError);
	action.setPointProxyPixelThreshold(
	    this->d->lodPresentationPointProxyPixelThreshold);
	/* Source/inventory discovery already supplied a leaf proxy.  Complete the
	 * useful-coverage census without launching one PoP build per visible leaf;
	 * a fresh budgeted quality pass follows immediately.  A zoom census keeps
	 * existing cuts and is therefore not structural-only. */
	action.setStructuralCoverageOnly(
	    this->d->lodCoveragePolicy.active() &&
	    !this->d->lodViewDemandPolicy.scaleDemandRefreshActive());
	action.setGeneration(generation);
	action.setRevisions(this->d->lodViewRevision.value(),
	    this->d->lodPolicyRevision.value());
	action.setRefreshMissing(refreshMissing);
	action.setReset(reset);
	action.setViewLodState(viewLodState);
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
	const SbBool applyRetainedAdmission =
	    this->d->lodBudgetPolicy.retainedAdmission() &&
	    (!this->d->lodInteractive || scaleInteraction);
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
	    !this->d->lodRefinementAwaitingFrame ? TRUE : FALSE);
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
	    !this->d->lodRefinementAwaitingFrame ? TRUE : FALSE);
	action.setAllowRepresentationRefinement(
	    scaleDemandChanged &&
	    !this->d->lodRefinementAwaitingFrame ? TRUE : FALSE);
	/* A calibrated face budget governs richness above the minimum drawable
	 * mesh floor.  Returning a visible resident PoP payload to its box is
	 * both a visual regression and a false economy: the renderer can batch
	 * tiny minimum prefixes into the aggregate point channel without
	 * discarding the useful retained asset. */
	action.setPreserveMeshCoverage(TRUE);
	action.setRefinementCostBudget(
	    this->d->lodBudgetPolicy.refinementRemaining());
	if (this->d->lodDiscretePopulationTrialAvailable)
	    action.setOneOverBudgetRefinementLimit(
		BObolLodQualityPolicy::discreteTrialOverBudgetAllowance(
		    sceneActiveCost,
		    this->d->lodBudgetPolicy.currentBudget()));
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
	action.setSubmissionTimeLimit(
	    this->d->forceTerminalLodRefinement ? 0 : 3000);
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
	const bool usingUnsatisfiedFrontier =
	    boundedLargeCompact &&
	!this->d->lodSubmissionDeltaActive &&
	!this->d->lodCoveragePolicy.active() &&
	this->d->lodCoveragePolicy.coverageComplete() &&
	    !this->d->lodBudgetPolicy.retainedAdmission() &&
	    viewLodState && sourceMeshRequests > 0 &&
	    viewLodState->cadMeshPayloadCountForSource(source) >=
		sourceMeshRequests;
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
		this->d->lodSubmissionPlanEntries);
	    this->d->lodSubmissionPlanSource = source;
	    this->d->lodSubmissionPlanValid = TRUE;
	    this->d->lodSubmissionPlanRetainedAdmission = TRUE;
	}
	if (usingUnsatisfiedFrontier &&
	    (!this->d->lodSubmissionPlanValid ||
	     this->d->lodSubmissionPlanSource != source)) {
	std::vector<SbString> unsatisfiedKeys;
	viewLodState->unsatisfiedCadOccurrenceKeys(
	    source,
	    this->d->lodService ?
		this->d->lodService->residentMeshAdmissionRevision() : 0,
	    unsatisfiedKeys);
	    controller_prioritize_lod_frontier(source, viewLodState,
		unsatisfiedKeys, this->d->lodSubmissionPlanEntries);
	    this->d->lodSubmissionPlanSource = source;
	    this->d->lodSubmissionPlanValid = TRUE;
	}
	/* The retained frontier is already ordered by visibility, error, emphasis,
	 * and value-per-cost.  In a quiet view, jump each occurrence to its richest
	 * budget-fitting resident prefix in that order.  Requiring one nominal PoP
	 * transition per whole-scene pass recreated level walking at 150k scale.
	 *
	 * Active zoom is different even before its first quality-frame probe has
	 * begun.  Residency is allowed to run ahead of the presentation budget, but
	 * asking one worker for the complete newly demanded suffix can turn a single
	 * large mesh into a multi-gigabyte working-set request.  Admission correctly
	 * rejects that request, leaving zoom refinement dependent on whether the
	 * quality-probe edge happens to run before the next camera event.  Bound all
	 * scale-interaction residency by the next population-bearing band.  This
	 * preserves direct quiet convergence while making active zoom deterministic
	 * and resource bounded. */
	if (scaleInteraction ||
	    this->d->lodViewDemandPolicy.qualityBudgetActive())
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
	action.apply(source);
	if (action.getOneOverBudgetRefinementUsed())
	    this->d->lodDiscretePopulationTrialAvailable = FALSE;
	const bool fullViewCandidatePass =
	    !usingUnsatisfiedFrontier &&
	    !this->d->lodSubmissionDeltaActive;
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
	if (getenv("BOBOL_LOD_TRACE_SUBMISSION")) {
	    static std::atomic<unsigned int> submissionTraceCount(0);
	    const unsigned int traceIndex =
		submissionTraceCount.fetch_add(1);
	    if (traceIndex < 128)
		bu_log("BObol LoD submission at=%lld path=%s entries=%d "
		       "offset=%zu visited=%u visible=%zu covered=%zu "
		       "tasks=%u cuts=%u skipped=%u deferred=%d next=%zu "
		       "coverage_pass=%d diagnostics=%s\n",
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
		       action.getDiagnostics().getString());
	}
	if (fullViewCandidatePass &&
	    !action.hasDeferredCompactEntries())
	    this->d->setLodConvergenceCandidateCount(source,
		this->d->lodSubmissionVisibleCount);
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
    const bool completedPassRetainedAllocation =
	completedPass &&
	this->d->lodSubmissionRetainedAdmissionMode;
    const bool completedPassChangedCut =
	completedPass &&
	this->d->lodRetainedRefinementCutAdvanced;
    const bool completedPassBudgetBlocked =
	completedPass &&
	this->d->lodRetainedRefinementBudgetBlocked;
    if (completedPass && getenv("BOBOL_LOD_TRACE_PASS")) {
	static std::atomic<unsigned int> completedPassTraceCount(0);
	const unsigned int traceIndex =
	    completedPassTraceCount.fetch_add(1);
	if (traceIndex < 512)
	    bu_log("BObol LoD completed pass coverage=%d "
		   "coverage_seen=%d coverage_complete=%d delta=%d "
		   "handoff=%d refinement_pending=%d cut_advanced=%d "
		   "budget_rescan=%d active_faces=%zu active_cost=%zu "
		   "cost_budget=%zu "
		   "visited=%u cuts=%u\n",
		   this->d->lodCoveragePolicy.active() ? 1 : 0,
		   this->d->lodCoveragePolicy.sawBoundedSource() ? 1 : 0,
		   this->d->lodCoveragePolicy.coverageComplete() ? 1 : 0,
		   this->d->lodSubmissionDeltaActive ? 1 : 0,
		   this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
		   this->d->lodRetainedRefinementPending ? 1 : 0,
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
	if (getenv("BOBOL_LOD_TRACE_PASS"))
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
	if (getenv("BOBOL_LOD_TRACE_PASS"))
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
	    this->notifyFrameRequest("lod-coverage-admission");
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
	    this->notifyFrameRequest("lod-coverage-minimum-calibration");
	} else {
	    /* Every projected leaf has a useful structural or mesh presentation. Begin a
	     * fresh bounded pass which may now spend the remaining scene budget
	     * on screen-value-ordered PoP refinement. */
	    this->d->lodSubmissionPending = TRUE;
	    this->markProgressiveWorkPending();
	}
	} else if (coverageCompletion.completed &&
	    !coverageCompletion.missing &&
	    this->d->lodRetainedRefinementPending) {
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
	if (this->d->lodSubmissionDeltaActive) {
	    this->d->lodSubmissionDeltaActive = FALSE;
	    this->d->lodSubmissionDeltaSources.clear();
	    this->d->lodSubmissionDeltaPlans.clear();
	}
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	if (this->d->lodCoveragePolicy.active())
	    this->d->lodCoveragePolicy.clearPassCounters();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = sources.empty() ? FALSE : TRUE;
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
	if (getenv("BOBOL_LOD_TRACE_PASS")) {
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
	    /* A large compact scan can outlive the render request installed by
	     * its first admitted result.  At pass completion there may then be
	     * no worker, result, or fresh progressive-work edge left to wake the
	     * host, even though completeRenderTiming() must consume this frame
	     * before admission may continue.  Install both halves of the frame
	     * contract explicitly. */
	    this->requestRender(frameReason);
	    this->notifyFrameRequest(frameReason);
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
	if (!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	    this->d->lodBudgetPolicy.stableBudgetLimited() &&
	    !this->d->lodInteractive &&
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
		this->notifyFrameRequest("lod-stable-point-calibration");
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
	!this->d->lodRefinementAwaitingFrame) {
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
    handoffInputs.retainedRefinementPending =
	this->d->lodRetainedRefinementPending != FALSE;
    handoffInputs.retainedRefinementBudgetBlocked =
	completedPassBudgetBlocked;
    const BObolLodPresentationPolicy::CompletedPassDecision handoff =
	this->d->lodPresentationPolicy.completePass(handoffInputs);
    if (completedPass && getenv("BOBOL_LOD_TRACE_PASS"))
	bu_log("BObol LoD handoff pass pending=%d rescan=%d changed=%d "
	       "retained=%d blocked=%d presentation_pending=%d "
	       "finish=%d request_rescan=%d retire=%d preserve=%d\n",
	       handoffInputs.submissionPending ? 1 : 0,
	       handoffInputs.rescanAfterFrame ? 1 : 0,
	       handoffInputs.changedCut ? 1 : 0,
	       handoffInputs.retainedRefinementPending ? 1 : 0,
	       handoffInputs.retainedRefinementBudgetBlocked ? 1 : 0,
	       this->d->lodPresentationPolicy.handoffPresentationPending() ? 1 : 0,
	       handoff.finishHandoff ? 1 : 0,
	       handoff.requestRetainedRescan ? 1 : 0,
	       handoff.retireRetainedObservation ? 1 : 0,
	       handoff.preservePresentationLimits ? 1 : 0);
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
	this->notifyFrameRequest("lod-stable-handoff");
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
	this->d->lodPointProxyTriangleRecoveryPending = FALSE;
	const SbBool pointCutChanged =
	    this->d->lodPresentationPointProxyPixelThreshold > 1.01f ?
		TRUE : FALSE;
	this->d->lodInteractiveProgressiveCeiling = -1;
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->lodPointProxyCalibrationPolicy.reset();
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	if (presentationState) {
	    presentationState->setCadPresentationProgressiveCutCeiling(-1);
	    presentationState->setCadPresentationPointProxyPixelThreshold(1.0f);
	}
	this->d->lodStablePointProxyCalibrationPending =
	    pointCutChanged && presentationState &&
	    presentationState->hasCadPresentationAssemblies() ? TRUE : FALSE;
	if (this->d->lodStablePointProxyCalibrationPending) {
	    this->requestRender("lod-stable-point-restore");
	    this->notifyFrameRequest("lod-stable-point-restore");
	} else {
	    this->resumeLodAfterOnePixelRecovery();
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
	!this->d->lodRefinementAwaitingFrame) {
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
     * no compatibility cost on the new CAD path and no duplicated source
     * geometry or presentation state.
     */
    unsigned int directMatched = 0;
    unsigned int directApplied = 0;
    unsigned int directRejected = 0;
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
	 * stale-epoch arbitration must not use the source-agnostic compatibility
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
	    if (drained[i].request.sourceRoutingId != 0)
		return route ? viewState->findCadForResult(route, drained[i]) :
		    static_cast<const BObolViewLodState::CadPayload *>(NULL);
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
	    residentCadBefore->sourceRevision ==
		drained[i].request.sourceRevision.value() &&
	    residentCadBefore->sourceContentHash ==
		drained[i].request.sourceContentHash &&
	    residentCadBefore->drawMode == drained[i].request.drawMode;
	const bool extendsCumulativeCadAsset =
	    sameCumulativeCadAsset &&
	    drained[i].residentCut > residentCadBefore->residentCut;
	/* A worker completion publishes residency, not a coarsening decision.
	 * The current aggregate allocator may already have selected a richer cut
	 * from the same cumulative asset while this suffix was in flight.  Keep
	 * that cut; camera/policy retargeting is the sole authority which lowers
	 * presentation.  A genuine suffix extension necessarily contains the
	 * older drawable prefix. */
	if (extendsCumulativeCadAsset &&
	    drained[i].geometry.activeCut < residentCadBefore->activeCut &&
	    drained[i].progressiveMesh->canDrawCut(
		residentCadBefore->activeCut)) {
	    drained[i].geometry.activeCut = residentCadBefore->activeCut;
	    drained[i].counts.faceCount =
		drained[i].progressiveMesh->faceCount(
		    drained[i].geometry.activeCut);
	    drained[i].counts.pointCount =
		drained[i].progressiveMesh->pointCount(
		    drained[i].geometry.activeCut);
	    drained[i].counts.originalPointCount = drained[i].counts.pointCount;
	    drained[i].counts.normalCount = drained[i].hasNormals ?
		drained[i].counts.faceCount * 3 : 0;
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
	    if (!drained[i].progressiveMesh->canDrawCut(activeCut))
		activeCut = residentCadBefore->activeCut;

	    drained[i].request.viewRevision =
		this->d->lodViewRevision.value();
	    drained[i].request.policyRevision =
		this->d->lodPolicyRevision.value();
	    drained[i].request.requestedCut =
		residentCadBefore->requestedCut;
	    drained[i].resolvedCut = residentCadBefore->requestedCut;
	    drained[i].cacheKey = bobol_lod_cache_key(drained[i].request);
	    drained[i].geometry.activeCut = activeCut;
	    drained[i].counts.faceCount =
		drained[i].progressiveMesh->faceCount(activeCut);
	    drained[i].counts.pointCount =
		drained[i].progressiveMesh->pointCount(activeCut);
	    drained[i].counts.originalPointCount =
		drained[i].counts.pointCount;
	    drained[i].counts.normalCount = drained[i].hasNormals ?
		drained[i].counts.faceCount * 3 : 0;
	    drained[i].terminal =
		drained[i].resolvedCut < 0 ||
		activeCut >= drained[i].resolvedCut ?
		    TRUE : FALSE;
	    drained[i].stale = FALSE;
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
	    drained[i].resolvedCut > availableActiveCut ? TRUE : FALSE;
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
	if (drained[i].request.sourceRoutingId != 0 &&
	    drained[i].request.occurrenceKey.getLength() > 0) {
	    if (route && route->hasCompactInstanceKey(
		    drained[i].request.occurrenceKey.getString())) {
		directMatched++;
		if (viewState &&
		    viewState->consumeSourceResult(
			route, drained[i])) {
		    directApplied++;
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
		"legacy LoD result application requires a scene root");
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
	    if (getenv("BOBOL_LOD_TRACE_PASS"))
		bu_log("BObol LoD resident growth noted applied=%zu "
		       "interactive=%d scale_changed=%d\n",
		       applied, this->d->lodInteractive ? 1 : 0,
		       this->d->lodViewDemandPolicy.interactionScaleChanged() ?
			   1 : 0);
	    /* The missing suffix may complete after this scale epoch has already
	     * spent its one quality probe.  Residency alone is not a visible
	     * improvement, so rearm exactly one bounded probe for the newly
	     * available population.  Its publication frame is requested below;
	     * advanceProgressiveWork() will not consume this latch until that
	     * frame has completed. */
	    (void)this->d->lodViewDemandPolicy.rearmAfterResidentGrowth(
		this->d->lodInteractive != FALSE);
	    this->markProgressiveWorkPending();
	}
	const int64_t now = bu_gettime();
	this->d->lodPublicationPolicy.noteApplied(applied, now);
	const bool firstUsefulMesh =
	    this->getActiveLodMeshPayloadCount() <= applied;
	const bool streamIdle =
	    service->queuedResultCountForGeneration(
		this->d->lodActiveGeneration) == 0 &&
	    service->activeTaskCountForGeneration(
		this->d->lodActiveGeneration) == 0;
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
	    !this->d->lodRefinementAwaitingFrame) {
	    /*
	     * Hold the next prefix immediately, but allow independent coverage
	     * results to accumulate until the adaptive publication deadline.
	     * scheduleLodRefinementFrame() below supplies the host wakeup once
	     * that deadline is reached.
	     */
	    this->d->lodRefinementAwaitingFrame = TRUE;
	    this->d->lodRefinementResumeAfterRenderSerial =
		this->d->renderCompletionSerial + 1;
	    if (this->d->lodRefinementResumeAfterRenderSerial == 0)
		this->d->lodRefinementResumeAfterRenderSerial = 1;
	}
	if (publicationDecision.requestFrame) {
	    this->requestRender("lod-result");
	    /* Result publication may arm the refinement barrier while the generic
	     * progressive latch is already set.  requestRender() records intent but
	     * is not itself an edge-triggered host wakeup; a cold Hubble zoom could
	     * therefore wait forever with render=1, refinement_frame=1, and no
	     * worker capable of producing another event. */
	    this->notifyFrameRequest("lod-result");
	}
	this->d->lodCompactionPolicy.requestAfter(bu_gettime(), 750000);
	(void)this->enforceMeshResidencyBudget();
	if (partialRefinementCandidate && publishNow)
	    this->scheduleLodRefinementFrame("lod-result");
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
    this->d->lodPolicyRevision.set(newRevision);
    this->d->lodStaticOverscanActive = FALSE;
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

BObolSceneController *
BObolViewController::getSceneController(void)
{
    return &this->d->sceneController;
}

const BObolSceneController *
BObolViewController::getSceneController(void) const
{
    return &this->d->sceneController;
}

BObolFeatureStore &
BObolViewController::features(void)
{
    return *this->d->featureStore;
}

const BObolFeatureStore &
BObolViewController::features(void) const
{
    return *this->d->featureStore;
}

BObolPolygonStore &
BObolViewController::polygons(void)
{
    return *this->d->polygonStore;
}

const BObolPolygonStore &
BObolViewController::polygons(void) const
{
    return *this->d->polygonStore;
}

BObolSelectionStore &
BObolViewController::selection(void)
{
    return *this->d->selectionStore;
}

const BObolSelectionStore &
BObolViewController::selection(void) const
{
    return *this->d->selectionStore;
}

SoViewport *
BObolViewController::getViewport(void)
{
    return this->d->viewport;
}

const SoViewport *
BObolViewController::getViewport(void) const
{
    return this->d->viewport;
}

void
BObolViewController::setRenderContextManager(SoDB::ContextManager *manager)
{
    SoGLRenderAction *action = this->d->renderManager ?
	this->d->renderManager->getGLRenderAction() : NULL;
    if (!action || action->getContextManager() == manager)
	return;
    /* The cached image renderer owns provider-specific context state.  Drop
     * it while the old provider is still alive rather than retaining a stale
     * provider through a host switch. */
    delete this->d->imageRenderer;
    this->d->imageRenderer = NULL;
    this->d->imageRendererManager = NULL;
    action->setContextManager(manager);
}

SoDB::ContextManager *
BObolViewController::getRenderContextManager(void) const
{
    SoGLRenderAction *action = this->d->renderManager ?
	this->d->renderManager->getGLRenderAction() : NULL;
    return action ? action->getContextManager() : NULL;
}

SoRenderManager *
BObolViewController::getRenderManager(void)
{
    return this->d->renderManager;
}

const SoRenderManager *
BObolViewController::getRenderManager(void) const
{
    return this->d->renderManager;
}

static int
find_edit_preview_child(SoGroup *group, const char *previewId)
{
    if (!group || !previewId)
	return -1;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(SoBRLEditPreview::getClassTypeId()))
	    continue;
	SoBRLEditPreview *preview = static_cast<SoBRLEditPreview *>(node);
	if (bu_strcmp(preview->previewId.getValue().getString(), previewId) == 0)
	    return i;
    }
    return -1;
}

static int
find_hud_label_overlay_child(SoGroup *group, const char *labelId)
{
    if (!group || !labelId)
	return -1;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(SoBRLHUDLabelOverlay::getClassTypeId()))
	    continue;
	SoBRLHUDLabelOverlay *label = static_cast<SoBRLHUDLabelOverlay *>(node);
	if (bu_strcmp(label->labelId.getValue().getString(), labelId) == 0)
	    return i;
    }
    return -1;
}

static int
find_line_layer_overlay_child(SoGroup *group, const char *overlayId)
{
    if (!group || !overlayId)
	return -1;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(SoBRLLineLayerOverlay::getClassTypeId()))
	    continue;
	SoBRLLineLayerOverlay *overlay = static_cast<SoBRLLineLayerOverlay *>(node);
	if (bu_strcmp(overlay->overlayId.getValue().getString(), overlayId) == 0)
	    return i;
    }
    return -1;
}

int
BObolViewController::replaceEditPreview(const char *previewId,
	const char *identity,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	uint32_t sourceRevision,
	uint32_t inputsRevision)
{
    return this->replaceEditPreviewWithIntent(previewId, identity, NULL, NULL,
	    points, commands, count, sourceRevision, inputsRevision);
}

int
BObolViewController::replaceEditPreviewWithIntent(const char *previewId,
	const char *identity,
	const char *editIntentId,
	const char *editIntentRole,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	uint32_t sourceRevision,
	uint32_t inputsRevision)
{
    if (!previewId || !previewId[0] || !points || !commands || count <= 0)
	return -1;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_edit_preview_child(group, previewId);
    SoBRLEditPreview *preview = NULL;
    if (childIndex >= 0)
	preview = static_cast<SoBRLEditPreview *>(group->getChild(childIndex));
    else
	preview = new SoBRLEditPreview;

    if (sourceRevision == 0)
	sourceRevision = preview->sourceRevision.getValue() + 1;
    if (inputsRevision == 0)
	inputsRevision = preview->inputsRevision.getValue() + 1;

    preview->previewId = previewId;
    preview->setEditIntent(editIntentId ? editIntentId : "",
			   editIntentRole ? editIntentRole : "preview");
    preview->sourceRevision = sourceRevision;
    preview->inputsRevision = inputsRevision;
    SoBRLVListShape *shape = preview->setLineSet(
				 (identity && identity[0]) ? identity : previewId,
				 points, commands, count);
    if (!shape)
	return -1;

    if (childIndex < 0)
	group->addChild(preview);

    this->requestRender("edit-preview");
    return 1;
}

int
BObolViewController::removeEditPreview(const char *previewId)
{
    if (!previewId || !previewId[0])
	return 0;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_edit_preview_child(group, previewId);
    if (childIndex < 0)
	return 0;

    group->removeChild(childIndex);
    this->requestRender("edit-preview");
    return 1;
}

int
BObolViewController::replaceLineLayerOverlay(const char *overlayId,
	const struct bg_line_layer_builder *builder,
	uint32_t sourceId,
	SbBool selectable)
{
    if (!overlayId || !overlayId[0])
	return -1;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_line_layer_overlay_child(group, overlayId);
    if (!builder) {
	if (childIndex >= 0) {
	    group->removeChild(childIndex);
	    this->requestRender("line-layer-overlay");
	}
	return 0;
    }

    SoBRLLineLayerOverlay *overlay = new SoBRLLineLayerOverlay;
    overlay->overlayId = overlayId;
    overlay->sourceId = sourceId;
    overlay->selectable = selectable;
    const int realized = overlay->rebuildGeometry(builder);

    if (childIndex >= 0)
	group->replaceChild(childIndex, overlay);
    else
	group->addChild(overlay);

    this->requestRender("line-layer-overlay");
    return realized;
}

int
BObolViewController::replaceHUDLabelOverlay(const char *labelId,
	const char *text,
	const SbVec2f &position,
	const SbColor &color,
	float fontSize,
	uint32_t sourceId)
{
    if (!labelId || !labelId[0])
	return -1;
    if (!text || !text[0])
	return this->removeHUDLabelOverlay(labelId);

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_hud_label_overlay_child(group, labelId);
    SoBRLHUDLabelOverlay *label = childIndex >= 0 ?
				  static_cast<SoBRLHUDLabelOverlay *>(group->getChild(childIndex)) :
				  new SoBRLHUDLabelOverlay;

    label->labelId = labelId;
    label->sourceId = sourceId;
    label->text = text;
    label->position = position;
    label->color = color;
    label->fontSize = fontSize;
    label->visible = TRUE;
    label->rebuildGeometry();

    if (childIndex < 0)
	group->addChild(label);

    this->requestRender("hud-label-overlay");
    return 1;
}

int
BObolViewController::removeHUDLabelOverlay(const char *labelId)
{
    if (!labelId || !labelId[0])
	return 0;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_hud_label_overlay_child(group, labelId);
    if (childIndex < 0)
	return 0;

    group->removeChild(childIndex);
    this->requestRender("hud-label-overlay");
    return 1;
}

SoGroup *
BObolViewController::findGroup(const char *groupPath) const
{
    return this->d->sceneController.findGroup(groupPath);
}

SoGroup *
BObolViewController::ensureGroup(const char *groupPath)
{
    const uint64_t revision = this->d->sceneController.getStructuralRevision();
    SoGroup *group = this->d->sceneController.ensureGroup(groupPath);
    if (group && this->d->sceneController.getStructuralRevision() != revision)
	this->requestRender("scene-group");
    return group;
}

int
BObolViewController::setGroupDrawIntent(const char *groupPath,
	const char *intentPath,
	int drawMode,
	int fallbackDrawMode,
	SbBool overlayIntent,
	uint32_t revalidationRevision)
{
    const int changed = this->d->sceneController.setGroupDrawIntent(groupPath,
			intentPath, drawMode, fallbackDrawMode, overlayIntent,
			revalidationRevision);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::setGroupDisplayState(const char *groupPath,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    const int changed = this->d->sceneController.setGroupDisplayState(
			    groupPath, visible, selected, highlighted, lineStyle,
			    lineWidth, transparency, colorOverride, color,
			    materialColorValid, materialColor, materialRevision);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::renameGroup(const char *groupPath,
				   const char *newLeafName)
{
    const int changed =
	this->d->sceneController.renameGroup(groupPath, newLeafName);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::appendChildToGroup(const char *groupPath,
	SoNode *child)
{
    const int changed =
	this->d->sceneController.appendChildToGroup(groupPath, child);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::removeChildFromGroup(const char *groupPath,
	SoNode *child)
{
    const int changed =
	this->d->sceneController.removeChildFromGroup(groupPath, child);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::eraseGroupSubpath(const char *parentGroupPath,
	const char *subpath)
{
    const int changed =
	this->d->sceneController.eraseGroupSubpath(parentGroupPath, subpath);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::removeGroup(const char *groupPath)
{
    const int removed = this->d->sceneController.removeGroup(groupPath);
    if (removed > 0)
	this->requestRender("scene-group");
    return removed;
}

int
BObolViewController::clearGroup(const char *groupPath)
{
    const int removed = this->d->sceneController.clearGroup(groupPath);
    if (removed > 0)
	this->requestRender("scene-group");
    return removed;
}

int
BObolViewController::getGroupChildCount(const char *groupPath) const
{
    return this->d->sceneController.getGroupChildCount(groupPath);
}

int
BObolViewController::getGroupDescendantGroupCount(
    const char *groupPath) const
{
    return this->d->sceneController.getGroupDescendantGroupCount(groupPath);
}

int
BObolViewController::getGroupDatabaseSourceCount(
    const char *groupPath) const
{
    return this->d->sceneController.getGroupDatabaseSourceCount(groupPath);
}

SoNode *
BObolViewController::findShape(const char *shapePath) const
{
    return this->d->sceneController.findShape(shapePath);
}

SoGroup *
BObolViewController::findShapeParent(const char *shapePath) const
{
    return this->d->sceneController.findShapeParent(shapePath);
}

int
BObolViewController::moveShapeToGroup(const char *shapePath,
					const char *groupPath)
{
    const int changed =
	this->d->sceneController.moveShapeToGroup(shapePath, groupPath);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BObolViewController::removeShape(const char *shapePath)
{
    const int removed = this->d->sceneController.removeShape(shapePath);
    if (removed > 0)
	this->requestRender("scene-shape");
    return removed;
}

int
BObolViewController::setShapeDrawState(const char *shapePath,
	int drawMode,
	SbBool databaseIntent,
	SbBool overlayIntent,
	SbBool hudIntent)
{
    const int changed = this->d->sceneController.setShapeDrawState(shapePath,
			drawMode, databaseIntent, overlayIntent, hudIntent);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BObolViewController::setShapeDisplayState(const char *shapePath,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    const int changed = this->d->sceneController.setShapeDisplayState(
			    shapePath, visible, selected, highlighted, lineStyle, lineWidth,
			    transparency, colorOverride, color, materialColorValid,
			    materialColor, materialRevision);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BObolViewController::setShapePlacementState(const char *shapePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize)
{
    const int changed = this->d->sceneController.setShapePlacementState(
			    shapePath, drawMatrixValid, drawMatrix, drawCenterValid,
			    drawCenter, drawSizeValid, drawSize);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BObolViewController::setShapeSourceState(const char *shapePath,
	const char *ownerSourcePath,
	uint32_t ownerSourceRevision,
	uint32_t ownerInputsRevision,
	uint32_t ownerViewRevision,
	uint32_t ownerRealizedRevision,
	uint32_t ownerRealizedSourceRevision,
	uint32_t ownerRealizedInputsRevision,
	uint32_t ownerRealizedViewRevision,
	int ownerRealizationStatus,
	const char *ownerRealizationDiagnostic,
	const char *ownerRealizationIdentity,
	SbBool ownerSourceStale,
	uint32_t ownerStaleReason)
{
    const int changed = this->d->sceneController.setShapeSourceState(
			    shapePath, ownerSourcePath, ownerSourceRevision,
			    ownerInputsRevision, ownerViewRevision, ownerRealizedRevision,
			    ownerRealizedSourceRevision, ownerRealizedInputsRevision,
			    ownerRealizedViewRevision, ownerRealizationStatus,
			    ownerRealizationDiagnostic, ownerRealizationIdentity,
			    ownerSourceStale, ownerStaleReason);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BObolViewController::replaceDatabaseSource(const char *sourcePath,
	struct db_i *dbip,
	int drawMode,
	uint32_t sourceRevision)
{
    int changed = this->d->sceneController.replaceDatabaseSource(sourcePath,
		  dbip, drawMode, sourceRevision);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::replaceDatabaseSourceInstance(
    const char *sourceInstanceKey,
    const char *sourcePath,
    struct db_i *dbip,
    int drawMode,
    uint32_t sourceRevision)
{
    int changed = this->d->sceneController.replaceDatabaseSourceInstance(
		      sourceInstanceKey, sourcePath, dbip, drawMode, sourceRevision);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceState(const char *sourcePath,
	SbBool sourceRevisionValid,
	uint32_t sourceRevision,
	uint32_t inputsRevision,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    SoBRLDatabaseSource *source =
	this->d->sceneController.findDatabaseSource(sourcePath);
    const uint32_t previousSourceRevision = source ?
	source->sourceRevision.getValue() : 0;
    const uint32_t previousInputsRevision = source ?
	source->inputsRevision.getValue() : 0;
    const SbBool previousStale = source ? source->stale.getValue() : FALSE;
    const uint32_t previousStaleReason = source ?
	source->staleReason.getValue() : 0;
    const int changed = this->d->sceneController.setDatabaseSourceState(
			    sourcePath, sourceRevisionValid, sourceRevision, inputsRevision,
			    visible, selected, highlighted, lineStyle, lineWidth, transparency,
			    colorOverride, color, materialColorValid, materialColor,
			    materialRevision);
    if (changed > 0) {
	source = this->d->sceneController.findDatabaseSource(sourcePath);
	if (!source ||
	    source->sourceRevision.getValue() != previousSourceRevision ||
	    source->inputsRevision.getValue() != previousInputsRevision ||
	    source->stale.getValue() != previousStale ||
	    source->staleReason.getValue() != previousStaleReason)
	    this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceInstanceState(
    const char *sourceInstanceKey,
    SbBool sourceRevisionValid,
    uint32_t sourceRevision,
    uint32_t inputsRevision,
    SbBool visible,
    SbBool selected,
    SbBool highlighted,
    int lineStyle,
    int lineWidth,
    float transparency,
    SbBool colorOverride,
    const SbColor &color,
    SbBool materialColorValid,
    const SbColor &materialColor,
    uint32_t materialRevision)
{
    SoBRLDatabaseSource *source =
	this->d->sceneController.findDatabaseSourceInstance(sourceInstanceKey);
    const uint32_t previousSourceRevision = source ?
	source->sourceRevision.getValue() : 0;
    const uint32_t previousInputsRevision = source ?
	source->inputsRevision.getValue() : 0;
    const SbBool previousStale = source ? source->stale.getValue() : FALSE;
    const uint32_t previousStaleReason = source ?
	source->staleReason.getValue() : 0;
    const int changed = this->d->sceneController.setDatabaseSourceInstanceState(
			    sourceInstanceKey, sourceRevisionValid, sourceRevision,
			    inputsRevision, visible, selected, highlighted, lineStyle, lineWidth,
			    transparency, colorOverride, color, materialColorValid,
			    materialColor, materialRevision);
    if (changed > 0) {
	source = this->d->sceneController.findDatabaseSourceInstance(
		 sourceInstanceKey);
	if (!source ||
	    source->sourceRevision.getValue() != previousSourceRevision ||
	    source->inputsRevision.getValue() != previousInputsRevision ||
	    source->stale.getValue() != previousStale ||
	    source->staleReason.getValue() != previousStaleReason)
	    this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceDisplayPatch(const char *sourcePath,
	const BObolDatabaseSourceDisplayPatch &patch)
{
    const int changed = this->d->sceneController.setDatabaseSourceDisplayPatch(
			    sourcePath, patch);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceInstanceDisplayPatch(
    const char *sourceInstanceKey,
    const BObolDatabaseSourceDisplayPatch &patch)
{
    const int changed =
	this->d->sceneController.setDatabaseSourceInstanceDisplayPatch(
	    sourceInstanceKey, patch);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceDisplayName(const char *sourcePath,
	const char *displayName)
{
    const int changed = this->d->sceneController.setDatabaseSourceDisplayName(
			    sourcePath, displayName);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceInstanceDisplayName(
    const char *sourceInstanceKey,
    const char *displayName)
{
    const int changed =
	this->d->sceneController.setDatabaseSourceInstanceDisplayName(
	    sourceInstanceKey, displayName);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceBoundsState(const char *sourcePath,
	SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax)
{
    const int changed = this->d->sceneController.setDatabaseSourceBoundsState(
			    sourcePath, boundsValid, boundsMin, boundsMax);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceInstanceBoundsState(
    const char *sourceInstanceKey,
    SbBool boundsValid,
    const SbVec3f &boundsMin,
    const SbVec3f &boundsMax)
{
    const int changed =
	this->d->sceneController.setDatabaseSourceInstanceBoundsState(
	    sourceInstanceKey, boundsValid, boundsMin, boundsMax);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceMaterialPolicy(
    const char *sourcePath,
    int materialPolicy)
{
    const int changed = this->d->sceneController.setDatabaseSourceMaterialPolicy(
			    sourcePath, materialPolicy);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceInstanceMaterialPolicy(
    const char *sourceInstanceKey,
    int materialPolicy)
{
    const int changed =
	this->d->sceneController.setDatabaseSourceInstanceMaterialPolicy(
	    sourceInstanceKey, materialPolicy);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::markDatabaseSourceStale(const char *sourcePath,
	uint32_t staleReason)
{
    const int changed = this->d->sceneController.markDatabaseSourceStale(
			    sourcePath, staleReason);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::markDatabaseSourceInstanceStale(
    const char *sourceInstanceKey,
    uint32_t staleReason)
{
    const int changed = this->d->sceneController.markDatabaseSourceInstanceStale(
			    sourceInstanceKey, staleReason);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::removeDatabaseSource(const char *sourcePath)
{
    int removed = this->d->sceneController.removeDatabaseSource(sourcePath);
    if (removed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return removed;
}

int
BObolViewController::removeDatabaseSourceInstance(
    const char *sourceInstanceKey)
{
    int removed = this->d->sceneController.removeDatabaseSourceInstance(
		      sourceInstanceKey);
    if (removed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return removed;
}

int
BObolViewController::moveDatabaseSourceToGroup(const char *sourcePath,
	const char *groupPath)
{
    int moved = this->d->sceneController.moveDatabaseSourceToGroup(sourcePath,
		groupPath);
    if (moved > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return moved;
}

int
BObolViewController::moveDatabaseSourceInstanceToGroup(
    const char *sourceInstanceKey,
    const char *groupPath)
{
    int moved = this->d->sceneController.moveDatabaseSourceInstanceToGroup(
		    sourceInstanceKey, groupPath);
    if (moved > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return moved;
}

int
BObolViewController::clearDatabaseSources(void)
{
    int removed = this->d->sceneController.clearDatabaseSources();
    if (removed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return removed;
}

SoBRLDatabaseSource *
BObolViewController::getDatabaseSource(int index) const
{
    return this->d->sceneController.getDatabaseSource(index);
}

int
BObolViewController::getDatabaseSourceCount(void) const
{
    return this->d->sceneController.getDatabaseSourceCount();
}

std::vector<SoBRLDatabaseSource *>
BObolViewController::getRenderDatabaseSources(void) const
{
    return controller_render_database_sources(this);
}

SoBRLDatabaseSource *
BObolViewController::findDatabaseSourceInstance(
    const char *sourceInstanceKey) const
{
    return this->d->sceneController.findDatabaseSourceInstance(
	       sourceInstanceKey);
}

SbBool
BObolViewController::getDatabaseSourceSummary(int index,
	BObolDatabaseSourceSummary &summary) const
{
    return this->d->sceneController.getDatabaseSourceSummary(index, summary);
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
    this->d->lodStaticOverscanActive = FALSE;
    this->d->lodHeadroomPolicy.cancelRetry();
    this->d->clearLodConvergenceCandidates();
    this->d->resetLodConvergenceFraction();
    this->d->lodBudgetPolicy.clearBudgetLimit();
    this->d->lodBudgetPolicy.clearRetainedQualityFloorBudget();
    this->d->lodBudgetPolicy.clearRetainedRecoveryCeiling();
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodPointProxyTriangleRecoveryPending = FALSE;
    this->d->lodPointProxyCalibrationPolicy.reset();
    this->d->lodPresentationPolicy.viewInvalidated();
    this->d->lodRetainPoseOccurrenceCuts = FALSE;
    this->d->lodRetainedImportanceCensusPending = FALSE;
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
	    this->d->lodInteractiveProgressiveCeiling,
	    this->d->lodPresentationPointProxyPixelThreshold,
	    controller_lod_presentation_population(snapshotState),
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
    this->d->lodPresentationPointProxyPixelThreshold =
	BObolLodQualityPolicy::pointProxyThreshold(
	    this->d->lodPresentationPointProxyPixelThreshold,
	    interactionTimingSample,
	    this->d->lodInteractiveTargetFps);
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
    if (responsiveCeiling >= 0 &&
	(this->d->lodInteractiveProgressiveCeiling < 0 ||
	 responsiveCeiling <
	    this->d->lodInteractiveProgressiveCeiling))
	this->d->lodInteractiveProgressiveCeiling =
	    responsiveCeiling;
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
	controller_lod_presentation_population(viewState);
    restoreInputs.currentProgressiveCeiling =
	this->d->lodInteractiveProgressiveCeiling;
    restoreInputs.currentPointProxyPixelThreshold =
	this->d->lodPresentationPointProxyPixelThreshold;
    const BObolLodPresentationPolicy::RestoreDecision restore =
	this->d->lodPresentationPolicy.restorePrior(restoreInputs);
    if (restore.apply) {
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
    this->d->lodInteractiveTargetFps = interactiveFps;
    this->d->lodStableTargetFps = stableFps;
    this->d->lodStaticOverscanActive = FALSE;
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
	this->d->lodRefinementAwaitingFrame;
    status.activeGeneration = this->d->lodActiveGeneration;
    status.submissionSourceIndex = this->d->lodSubmissionSourceIndex;
    status.submissionEntryOffset = this->d->lodSubmissionEntryOffset;
    status.budgetCalibrationPending =
	this->d->lodBudgetPolicy.rescanAfterFrame() ||
	this->d->lodHeadroomPolicy.retryPending();
    status.stablePresentationHandoffPending =
	this->d->lodPresentationPolicy.handoffPending() ? TRUE : FALSE;
    status.pointProxyCalibrationPending =
	(this->d->lodStablePointProxyCalibrationPending ||
	 this->d->lodPointProxyTriangleRecoveryPending) ? TRUE : FALSE;
    status.residentGrowthReallocationPending =
	this->d->lodResidentGrowthPolicy.pending() ? TRUE : FALSE;
    status.publicationFramePending = publicationPending;

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
    convergenceInputs.stableBudgetLimited =
	this->d->lodBudgetPolicy.stableBudgetLimited();
    convergenceInputs.presentationLimited =
	this->d->lodInteractiveProgressiveCeiling >= 0 ||
	this->d->lodPresentationPointProxyPixelThreshold > 1.01f;
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
	    BOBOL_LOD_CONVERGENCE_ERROR,
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
		    this->d->lodInteractiveProgressiveCeiling,
		    this->d->lodPresentationPointProxyPixelThreshold,
		    controller_lod_presentation_population(snapshotState),
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
		!reusePriorScalePresentation)
		this->d->lodPresentationPointProxyPixelThreshold =
		    BObolLodQualityPolicy::pointProxyThreshold(
			this->d->lodPresentationPointProxyPixelThreshold,
			interactiveRenderSample,
			this->d->lodInteractiveTargetFps);
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
