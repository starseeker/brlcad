/*          V I E W _ C O N T R O L L E R _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_VIEW_CONTROLLER_PRIVATE_H
#define LIBBOBOL_VIEW_CONTROLLER_PRIVATE_H

#include "BObol/BDatabaseSource.h"
#include "BObol/BSceneController.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "identity_counter_private.h"
#include "lod_telemetry_private.h"
#include "view_lod_coordinator_state_private.h"
#include "bv.h"

#include <Inventor/SoDB.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>

#include <atomic>
#include <array>
#include <cstdint>
#include <exception>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

class SoCamera;
class SoClipPlane;
class SoDirectionalLight;
class SoGroup;
class SoGLRenderAction;
class SoNode;
class SoOffscreenRenderer;
class SoBRLViewLodGroup;

/* Keep the diagnostic transition stack's value type outside the exported
 * controller pimpl.  libstdc++ may emit an out-of-line vector growth helper;
 * nesting this type in BObolViewController::Impl consequently exposed a
 * private implementation name from libBObol. */
struct BObolLodControlTransitionFrame {
    uint64_t token = 0;
    BObolLodControlTransitionEvent event =
	BOBOL_LOD_CONTROL_TRANSITION_UNNAMED;
    SbBool ownerEvent = FALSE;
    SbBool suppressed = FALSE;
    BObolLodControlTraceState state;
};

std::vector<SoBRLDatabaseSource *> controller_render_database_sources(
    const BObolViewController *controller);
std::vector<SoBRLDatabaseSource *> controller_render_database_source_roots(
    const BObolViewController *controller);
bool controller_lod_structural_population_settled(
    const std::vector<SoBRLDatabaseSource *> &sources,
    bool externalProducersSettled);
bool controller_lod_source_inputs_unsubmitted(
    const std::vector<SoBRLDatabaseSource *> &sources,
    const std::vector<BObolLodSourceSnapshot> &submitted);
bool controller_lod_mesh_first_scene_safe(
    const std::vector<SoBRLDatabaseSource *> &sources);
bool controller_lod_adaptive_point_aggregation_allowed(
    const BObolViewController *controller,
    bool staticQualityCapacityRejected);
const char *controller_database_id(const struct db_i *dbip);
bool controller_lod_trace_enabled(const char *name, uint64_t viewRevision);
BObolLodPresentationPolicy::Population controller_lod_presentation_population(
    const BObolViewLodState *state, uint64_t sceneDomainRevision);
std::shared_ptr<BObolLodService> controller_acquire_managed_lod_service(
    size_t workerCount);
SbVec3f bobol_headlight_default_offset(void);
double controller_aspect_from_region(const SbViewportRegion &region);
SbViewportRegion controller_viewport_region_with_size(
    const SbViewportRegion &region, unsigned int width, unsigned int height);
void controller_configure_render_environment(SoViewport *viewport);
int controller_camera_root_index(SoViewport *viewport, SoCamera *camera);
SoClipPlane *controller_clip_plane(SoViewport *viewport, SbBool minimum);
SoClipPlane *controller_cutting_plane(SoViewport *viewport);
bool controller_camera_relative_clip_planes(const double center[3],
    const double viewZ[3], double horizontalSize, double minimum,
    double maximum, SbPlane &minimumPlane, SbPlane &maximumPlane);
static constexpr size_t controller_camera_light_count = 3;
std::array<SoDirectionalLight *, controller_camera_light_count>
controller_camera_lights(SoViewport *viewport);
std::array<SbVec3f, controller_camera_light_count>
controller_camera_light_directions(const SbRotation &orientation,
    const SbVec3f &headlightOffset);
bool controller_cutting_plane_affordance_update_needed(
    SoGroup *presentationRoot, SbBool enabled, bool cameraGeometryChanged);

class BObolPreparedCuttingPlaneAffordance {
public:
    BObolPreparedCuttingPlaneAffordance(SoViewport *viewport,
	SoGroup *presentationRoot, SoCamera *camera, const SbPlane &plane,
	SbBool enabled, double horizontalSize, double aspect);
    ~BObolPreparedCuttingPlaneAffordance();
    BObolPreparedCuttingPlaneAffordance(
	const BObolPreparedCuttingPlaneAffordance &) = delete;
    BObolPreparedCuttingPlaneAffordance &operator=(
	const BObolPreparedCuttingPlaneAffordance &) = delete;
    void commit();
    void notify(std::exception_ptr &failure);
private:
    class Impl;
    std::unique_ptr<Impl> impl;
};
/* The visible section aid is deliberately presentation-only.  It lives in the
 * post-CAD overlay, so the plane that defines the cut cannot clip or occlude
 * its own user cue. */
void controller_update_cutting_plane_affordance(SoViewport *viewport,
	SoGroup *presentationRoot,
	const SbPlane &plane,
	SbBool enabled,
	double horizontalSize,
	double aspect);

/* One level-triggered host render request.  Semantic presentation, LoD
 * planning presentation, and capacity sampling are increasing strengths of
 * the same request, not independent latches.  A weaker repaint cannot
 * downgrade a stronger transaction, and consuming or clearing it retires the
 * diagnostic reason at the same transition.  The controller's
 * renderRequestMutex owns this value. */
class BObolRenderRequestState {
public:
    enum class Kind : uint8_t {
	NONE = 0,
	PRESENTATION,
	LOD_PRESENTATION,
	CAPACITY_SAMPLE
    };

    struct Decision {
	bool changed = false;
	bool wakeEndpoint = false;
    };

    Decision request(SbString &&reason, bool capacityRelevant,
	bool planningRelevant) noexcept
    {
	Decision decision;
	const Kind requested = capacityRelevant ?
	    Kind::CAPACITY_SAMPLE : planningRelevant ?
	    Kind::LOD_PRESENTATION : Kind::PRESENTATION;
	if (static_cast<uint8_t>(requested) < static_cast<uint8_t>(this->kindValue))
	    return decision;
	decision.changed = this->kindValue != requested;
	decision.wakeEndpoint = this->kindValue == Kind::NONE;
	/* Same-strength requests coalesce but retain the newest diagnostic. */
	this->reasonValue = std::move(reason);
	this->kindValue = requested;
	return decision;
    }

    bool retireCapacity(void)
    {
	if (this->kindValue != Kind::CAPACITY_SAMPLE)
	    return false;
	this->kindValue = Kind::PRESENTATION;
	return true;
    }

    bool clear(void)
    {
	const bool changed = this->kindValue != Kind::NONE ||
	    this->reasonValue.getLength() > 0;
	this->kindValue = Kind::NONE;
	this->reasonValue = "";
	return changed;
    }

    bool consume(SbString *reason, SbBool *capacityRelevant,
	SbBool *planningRelevant)
    {
	const bool wasPending = this->pending();
	if (reason)
	    *reason = this->reasonValue;
	if (capacityRelevant)
	    *capacityRelevant = this->capacityRelevant() ? TRUE : FALSE;
	if (planningRelevant)
	    *planningRelevant = this->planningRelevant() ? TRUE : FALSE;
	this->kindValue = Kind::NONE;
	this->reasonValue = "";
	return wasPending;
    }

    bool pending(void) const
    {
	return this->kindValue != Kind::NONE;
    }

    bool capacityRelevant(void) const
    {
	return this->kindValue == Kind::CAPACITY_SAMPLE;
    }

    bool planningRelevant(void) const
    {
	return this->kindValue == Kind::LOD_PRESENTATION ||
	    this->kindValue == Kind::CAPACITY_SAMPLE;
    }

    const SbString &reason(void) const
    {
	return this->reasonValue;
    }

private:
    Kind kindValue = Kind::NONE;
    SbString reasonValue;
};

/* The standing request and its callback each retain a prepared reason. This
 * also makes aliasing the preceding request's reason safe across replacement. */
struct BObolPreparedRenderRequest {
    BObolPreparedRenderRequest() = default;
    BObolPreparedRenderRequest(const char *message, bool capacity, bool planning) :
	reason(message ? message : ""), notificationReason(reason),
	capacityRelevant(capacity), planningRelevant(planning)
    {}

    SbString reason;
    SbString notificationReason;
    bool capacityRelevant = false;
    bool planningRelevant = false;
    BObolRenderRequestState::Decision decision;
    uint64_t serial = 0;
};

struct BObolViewController::PreparedRendererInvalidation {
    BObolPreparedRenderRequest renderRequest;
    bool automatic = false;
    bool requestCommitted = false;
};

struct BObolViewController::PreparedRenderContextPublication {
    PreparedRendererInvalidation rendererInvalidation;
    SoGLRenderAction *action = nullptr;
    SoDB::ContextManager *previousManager = nullptr;
    SoDB::ContextManager *nextManager = nullptr;
    bool changed = false;
    bool bound = false;
    bool complete = false;
};

struct BObolViewController::PreparedViewportPublication {
    BObolPreparedRenderRequest renderRequest;
    bool changed = false;
};

/* Diagnostic RAII boundary around one production effect transaction.  An
 * owner scope derives its event from the currently selected refinement owner;
 * explicit scopes are reserved for external input publication.  The outermost
 * scope is the commit boundary.  Nested scopes are implementation helpers and
 * must not expose a partially updated control ledger as a model state. */
class BObolLodControlTransitionScope {
public:
    explicit BObolLodControlTransitionScope(BObolViewController *controller) :
	controllerValue(controller),
	tokenValue(begin(controller, BOBOL_LOD_CONTROL_TRANSITION_UNNAMED, TRUE))
    {
    }

    BObolLodControlTransitionScope(BObolViewController *controller,
	BObolLodControlTransitionEvent event) :
	controllerValue(controller),
	tokenValue(begin(controller, event, FALSE))
    {
    }

    ~BObolLodControlTransitionScope()
    {
	if (!this->controllerValue || this->tokenValue == 0) return;
	try { this->controllerValue->endLodControlTransition(this->tokenValue); }
	catch (...) { this->controllerValue->discardLodControlTransition(); }
    }

    BObolLodControlTransitionScope(
	const BObolLodControlTransitionScope &) = delete;
    BObolLodControlTransitionScope &operator=(
	const BObolLodControlTransitionScope &) = delete;

private:
    static uint64_t begin(BObolViewController *controller,
	BObolLodControlTransitionEvent event, SbBool ownerEvent)
    {
	if (!controller) return 0;
	try { return controller->beginLodControlTransition(event, ownerEvent); }
	catch (...) {
	    controller->discardLodControlTransition();
	    throw;
	}
    }
    BObolViewController *controllerValue;
    uint64_t tokenValue;
};

/* Retire the host's in-flight frame witness immediately before the enclosing
 * completed/interrupted-frame transition is recorded. */
class BObolRenderClaimCompletionScope {
public:
    explicit BObolRenderClaimCompletionScope(
	BObolViewController *controller) : controllerValue(controller)
    {
    }

    ~BObolRenderClaimCompletionScope()
    {
	if (this->controllerValue)
	    this->controllerValue->retireClaimedRender();
    }

    BObolRenderClaimCompletionScope(
	const BObolRenderClaimCompletionScope &) = delete;
    BObolRenderClaimCompletionScope &operator=(
	const BObolRenderClaimCompletionScope &) = delete;

private:
    BObolViewController *controllerValue;
};

/* Uninstalled state shared by the responsibility-specific controller units.
 * It retains direct storage and one allocation: extracting scene forwarding
 * must not add virtual dispatch or accessors to LoD hot paths. */
struct BObolViewController::Impl : BObolLodCoordinator {
    explicit Impl(BObolViewController *owner) :
	controller(owner),
	headlightOffsetEye(bobol_headlight_default_offset())
    {
	auto viewportOwner = std::make_unique<SoViewport>();
	auto attachmentOwner = std::make_unique<BObolViewAttachment>();
	auto managerOwner = std::make_unique<SoRenderManager>();
	auto featureOwner = std::make_unique<BObolFeatureStore>(owner);
	auto polygonOwner = std::make_unique<BObolPolygonStore>(owner);
	auto selectionOwner = std::make_unique<BObolSelectionStore>();
	/* SoViewport has its own library default. Establish the controller's
	 * three-copy region invariant before any public constructor returns. */
	viewportOwner->setViewportRegion(this->viewportRegion);
	managerOwner->setViewportRegion(this->viewportRegion);
	this->viewport = viewportOwner.release();
	this->viewAttachment = attachmentOwner.release();
	this->renderManager = managerOwner.release();
	this->featureStore = featureOwner.release();
	this->polygonStore = polygonOwner.release();
	this->selectionStore = selectionOwner.release();
    }

    void requireExactPresentationFrame(void)
    {
	const uint64_t mutationNanoseconds = this->controller ?
	    this->controller->beginRenderTiming() : 0;
	this->requireExactPresentationFrameAfter(mutationNanoseconds);
	if (this->controller)
	    this->controller->markProgressiveWorkPending();
    }

    /* Renderer-only presentation limits are mirrored in BObolViewLodState.
     * Keep publication in one owner operation: a controller value without
     * the matching renderer value is not a valid intermediate state. */
    void requireCadPresentationCommitIfChanged(
	BObolViewLodState *state,
	const Obol::CadViewState &previous)
    {
	if (state && state->hasCadPresentationAssemblies() &&
	    previous != state->cadPresentationViewState())
	    this->requireExactPresentationFrame();
    }

    void publishCadProgressiveCeiling(int ceiling,
	float nextFraction = 0.0f)
    {
	ceiling = ceiling < 0 ? -1 : std::min<int>(
	    BOBOL_MESH_LOD_CUT_COUNT_MAX - 1, ceiling);
	nextFraction = ceiling < 0 || !std::isfinite(nextFraction) ? 0.0f :
	    std::max(0.0f, std::min(1.0f, nextFraction));
	this->lodInteractiveProgressiveCeiling = ceiling;
	BObolViewLodState *state = this->viewAttachment ?
	    this->viewAttachment->getViewLodState() : NULL;
	if (state) {
	    const Obol::CadViewState previous =
		state->cadPresentationViewState();
	    state->setCadPresentationProgressiveCutCeiling(
		ceiling, nextFraction);
	    this->requireCadPresentationCommitIfChanged(state, previous);
	}
    }

    bool rendererCeilingQualityPending(void) const
    {
	if (!this->lodAutoSubmit || this->forceTerminalLodRefinement ||
	    this->lodInteractionSession.active() ||
	    this->lodInteractiveProgressiveCeiling < 0)
	    return false;
	const BObolViewLodState *state = this->viewAttachment ?
	    this->viewAttachment->getViewLodState() : NULL;
	if (!state || state->maximumActiveProgressiveCut() <=
		this->lodInteractiveProgressiveCeiling)
	    return false;
	const BObolLodAdmissionRevisionStamp stamp =
	    this->admissionRevisionStamp();
	const float nextFraction =
	    state->cadPresentationProgressiveCutNextFraction();
	return !this->lodStaticQualityTrial.retainsRendererCeilingFor(
		stamp, this->lodInteractiveProgressiveCeiling) &&
	    !this->lodStaticQualityTrial.acceptedRendererCeilingFor(
		stamp, this->lodInteractiveProgressiveCeiling, nextFraction);
    }

    static float normalizedCadPointProxyThreshold(float threshold)
    {
	if (!std::isfinite(threshold) ||
	    threshold < POINT_PROXY_PIXEL_THRESHOLD_MINIMUM)
	    return POINT_PROXY_PIXEL_THRESHOLD_MINIMUM;
	return threshold > POINT_PROXY_PIXEL_THRESHOLD_MAXIMUM ?
	    POINT_PROXY_PIXEL_THRESHOLD_MAXIMUM : threshold;
    }

    void publishCadPointProxyThreshold(float threshold)
    {
	threshold = normalizedCadPointProxyThreshold(threshold);
	const bool changed = std::fabs(
	    this->lodPresentationPointProxyPixelThreshold - threshold) >
	    1.0e-6f;
	this->lodPresentationPointProxyPixelThreshold = threshold;
	BObolViewLodState *state = this->viewAttachment ?
	    this->viewAttachment->getViewLodState() : NULL;
	if (state) {
	    const Obol::CadViewState previous =
		state->cadPresentationViewState();
	    state->setCadPresentationPointProxyPixelThreshold(threshold);
	    this->requireCadPresentationCommitIfChanged(state, previous);
	}
	/* The point classifier is an allocator input.  The formal capacity
	 * machine treats it as part of its immutable revision tuple, so publish
	 * that edge here rather than relying on every calibration caller to
	 * remember a matching revision update. */
	if (changed)
	    this->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::CAPACITY);
    }

    void publishCadDiscoveryPointProxyThreshold(float threshold)
    {
	threshold = normalizedCadPointProxyThreshold(threshold);
	const bool changed = std::fabs(
	    this->lodDiscoveryPointProxyPixelThreshold - threshold) > 1.0e-6f;
	this->lodDiscoveryPointProxyPixelThreshold = threshold;
	BObolViewLodState *state = this->viewAttachment ?
	    this->viewAttachment->getViewLodState() : NULL;
	if (state) {
	    const Obol::CadViewState previous =
		state->cadPresentationViewState();
	    state->setCadPresentationDiscoveryPointProxyPixelThreshold(
		threshold);
	    this->requireCadPresentationCommitIfChanged(state, previous);
	}
	if (changed)
	    this->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::CAPACITY);
    }

    void publishCadPointProxyThresholds(float presentationThreshold,
	float discoveryThreshold)
    {
	presentationThreshold = normalizedCadPointProxyThreshold(
	    presentationThreshold);
	discoveryThreshold = normalizedCadPointProxyThreshold(
	    discoveryThreshold);
	const bool changed = std::fabs(
	    this->lodPresentationPointProxyPixelThreshold -
		presentationThreshold) > 1.0e-6f ||
	    std::fabs(this->lodDiscoveryPointProxyPixelThreshold -
		discoveryThreshold) > 1.0e-6f;
	this->lodPresentationPointProxyPixelThreshold = presentationThreshold;
	this->lodDiscoveryPointProxyPixelThreshold = discoveryThreshold;
	BObolViewLodState *state = this->viewAttachment ?
	    this->viewAttachment->getViewLodState() : NULL;
	if (state) {
	    const Obol::CadViewState previous =
		state->cadPresentationViewState();
	    state->setCadPresentationPointProxyPixelThreshold(
		presentationThreshold);
	    state->setCadPresentationDiscoveryPointProxyPixelThreshold(
		discoveryThreshold);
	    this->requireCadPresentationCommitIfChanged(state, previous);
	}
	if (changed)
	    this->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::CAPACITY);
    }

    void publishCadCameraMotionFrameReuse(SbBool enabled)
    {
	BObolViewLodState *state = this->viewAttachment ?
	    this->viewAttachment->getViewLodState() : NULL;
	if (!state)
	    return;
	const Obol::CadViewState previous =
	    state->cadPresentationViewState();
	state->setCadPresentationCameraMotionFrameReuse(enabled);
	this->requireCadPresentationCommitIfChanged(state, previous);
    }

    void publishCadPresentationLimits(int ceiling, float nextFraction,
	float pointProxyThreshold)
    {
	this->publishCadProgressiveCeiling(ceiling, nextFraction);
	this->publishCadPointProxyThreshold(pointProxyThreshold);
    }

    void resetCadPresentationLimits(void)
    {
	this->publishCadPresentationLimits(-1, 0.0f,
	    POINT_PROXY_PIXEL_THRESHOLD_MINIMUM);
    }

    void retireDisplayEndpointWorkNoexcept(void) noexcept
    {
	/* Destruction cannot depend on diagnostic journaling or observer
	 * delivery. Attempt each presentation reset independently, then retire
	 * every control owner and host-work level unconditionally. */
	try { this->resetCadPresentationLimits(); }
	catch (...) {}
	try { this->controller->resetDiscoveryPointProxyFloor(FALSE); }
	catch (...) {}
	try { this->publishCadCameraMotionFrameReuse(FALSE); }
	catch (...) {}
	this->retireAutomaticLodControlNoexcept();

	std::lock_guard<std::mutex> lock(this->renderRequestMutex);
	const bool renderRetired = this->renderRequest.clear();
	const bool frameRetired = renderRetired || this->renderClaimed;
	const bool changed = frameRetired || this->progressiveWorkPending;
	this->renderClaimed = FALSE;
	this->capacitySampleClaimed = FALSE;
	this->progressiveWorkPending = FALSE;
	if (changed) {
	    bobol_identity_advance(this->hostWorkRevision);
	    if (frameRetired)
		bobol_identity_advance(this->renderRequestSerial);
	}
    }

    static constexpr float POINT_PROXY_PIXEL_THRESHOLD_MINIMUM = 1.0f;
    static constexpr float POINT_PROXY_PIXEL_THRESHOLD_MAXIMUM = 64.0f;

    BObolSceneController sceneController;
    SoViewport *viewport = NULL;
    BObolViewController *controller = NULL;
    SoBRLViewLodGroup *renderLodRoot = NULL;
    SoNode *renderBatchRoot = NULL;
    SoNode *renderPresentationRoot = NULL;
    SoGroup *framebufferUnderlayRoot = NULL;
    SoGroup *framebufferInterlayRoot = NULL;
    SoGroup *framebufferOverlayRoot = NULL;
    BObolViewAttachment *viewAttachment = NULL;
    SoRenderManager *renderManager = NULL;
    SoOffscreenRenderer *imageRenderer = NULL;
    SoDB::ContextManager *imageRendererManager = NULL;
    SoCamera *activeCamera = NULL;
    SbViewportRegion viewportRegion = SbViewportRegion(1, 1);
    SbColor backgroundBottom = SbColor(0.0f, 0.0f, 0.0f);
    SbColor backgroundTop = SbColor(0.0f, 0.0f, 0.0f);
    BObolViewController::SoftwareWireMode softwareWireMode =
	BObolViewController::SOFTWARE_WIRE_AUTO;
    std::atomic<int> endpointGraphicalRenderingEnabled {1};
    std::atomic<int> endpointGraphicalRenderingSyncPending {0};
    SbBool transparencyEnabled = TRUE;
    SbBool antialiasingEnabled = FALSE;
    BObolViewController::LightingProfile lightingProfile =
	BObolViewController::LIGHTING_STUDIO;
    SbBool headlightEnabled = TRUE;
    SbBool headlightCameraTracked = TRUE;
    SbBool sceneLightsEnabled = FALSE;
    SbVec3f headlightOffsetEye;
    SbRotation lastCameraOrientation = SbRotation::identity();
    std::vector<BObolSceneLightRealization> sceneLights;
    /* In orthographic views the image scale may become arbitrarily small
     * without changing the scene depth.  Retain the largest scale observed
     * for this scene structure so camera clipping cannot turn a zoom into an
     * implicit section cut. */
    double orthographicDepthReferenceSize = 0.0;
    uint64_t orthographicDepthReferenceStructuralRevision = 0;
    double clipMinimum = BV_VIEW_MIN;
    double clipMaximum = BV_VIEW_MAX;
    SbBool cuttingPlaneEnabled = FALSE;
    SbPlane cuttingPlane = SbPlane(SbVec3f(0.0f, 0.0f, 1.0f),
	SbVec3f(0.0f, 0.0f, 0.0f));
    /* Last camera synchronization gives the retained plane aid a bounded
     * world-space extent without deriving a second camera state. */
    double cuttingPlaneAffordanceHorizontalSize = 1.0;
    double cuttingPlaneAffordanceAspect = 1.0;
    mutable std::mutex renderRequestMutex;
    SbBool progressiveWorkPending = FALSE;
    BObolRenderRequestState renderRequest;
    SbBool renderClaimed = FALSE;
    SbBool capacitySampleClaimed = FALSE;
    uint64_t hostWorkRevision = 0;
    uint64_t renderRequestSerial = 0;
    std::mutex frameRequestMutex;
    BObolFrameRequestCallback frameRequestCallback = NULL;
    void *frameRequestUserData = NULL;
    mutable std::mutex presentationSyncMutex;
    BObolPresentationSyncCallback presentationSyncCallback = NULL;
    void *presentationSyncUserData = NULL;
    /* The owner-thread transition capture feeds an optional bounded in-memory
     * journal and the environment-selected streaming telemetry sink. */
    bool lodControlTransitionActive(void) const
    {
	return this->lodControlTransitionJournalEnabled != FALSE ||
	    this->lodTelemetry.get() != NULL;
    }
    SbBool lodControlTransitionJournalEnabled = FALSE;
    size_t lodControlTransitionRecordLimit =
	BOBOL_LOD_CONTROL_TRACE_DEFAULT_RECORD_LIMIT;
    uint64_t lodControlTransitionNextToken = 1;
    uint64_t lodControlTransitionNextSerial = 1;
    uint64_t lodControlTransitionDropped = 0;
    SbBool lodControlTransitionHasEndpoint = FALSE;
    BObolLodControlTraceState lodControlTransitionEndpoint;
    std::vector<BObolLodControlTransitionFrame> lodControlTransitionFrames;
    std::vector<BObolLodControlTransitionRecord> lodControlTransitionRecords;
    std::unique_ptr<BObolLodTelemetry> lodTelemetry;
    /* Worker callbacks cannot inspect Coin/controller state.  They publish
     * only the finite event kind; the owner thread captures both endpoints at
     * its next journal boundary. */
    std::atomic<int> lodControlPendingExternalEvent {
	BOBOL_LOD_CONTROL_TRANSITION_UNNAMED};
    BObolFeatureStore *featureStore = NULL;
    BObolPolygonStore *polygonStore = NULL;
    BObolSelectionStore *selectionStore = NULL;
};

#endif /* LIBBOBOL_VIEW_CONTROLLER_PRIVATE_H */
