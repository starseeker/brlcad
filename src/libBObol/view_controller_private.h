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
#include "view_lod_coordinator_state_private.h"
#include "bv.h"

#include <Inventor/SoDB.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>

#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <vector>

class SoCamera;
class SoClipPlane;
class SoGroup;
class SoNode;
class SoOffscreenRenderer;
class SoBRLViewLodGroup;

std::vector<SoBRLDatabaseSource *> controller_render_database_sources(
    const BObolViewController *controller);
std::vector<SoBRLDatabaseSource *> controller_render_database_source_roots(
    const BObolViewController *controller);
const char *controller_database_id(const struct db_i *dbip);
bool controller_lod_trace_enabled(const char *name, uint64_t viewRevision);
BObolLodPresentationPolicy::Population controller_lod_presentation_population(
    const BObolViewLodState *state, uint64_t sceneDomainRevision);
std::shared_ptr<BObolLodService> controller_acquire_managed_lod_service(
    size_t workerCount);
SbVec3f bobol_headlight_default_offset(void);
double controller_aspect_from_region(const SbViewportRegion &region);
void controller_configure_render_environment(SoViewport *viewport);
SoClipPlane *controller_clip_plane(SoViewport *viewport, SbBool minimum);
SoClipPlane *controller_cutting_plane(SoViewport *viewport);
/* The visible section aid is deliberately presentation-only.  It lives in the
 * post-CAD overlay, so the plane that defines the cut cannot clip or occlude
 * its own user cue. */
void controller_update_cutting_plane_affordance(SoViewport *viewport,
	SoGroup *presentationRoot,
	const SbPlane &plane,
	SbBool enabled,
	double horizontalSize,
	double aspect);

/* One level-triggered host render request.  Capacity sampling is a stronger
 * form of the same request, not an independent latch: a weaker repaint cannot
 * downgrade it, and consuming or clearing it retires the diagnostic reason at
 * the same transition.  The controller's renderRequestMutex owns this value. */
class BObolRenderRequestState {
public:
    enum class Kind : uint8_t {
	NONE = 0,
	PRESENTATION,
	CAPACITY_SAMPLE
    };

    struct Decision {
	bool changed = false;
	bool wakeEndpoint = false;
    };

    Decision request(const char *reason, bool capacityRelevant)
    {
	Decision decision;
	const Kind requested = capacityRelevant ?
	    Kind::CAPACITY_SAMPLE : Kind::PRESENTATION;
	if (this->kindValue == Kind::NONE) {
	    this->kindValue = requested;
	    this->reasonValue = reason ? reason : "";
	    decision.changed = true;
	    decision.wakeEndpoint = true;
	} else if (this->kindValue == Kind::PRESENTATION &&
		requested == Kind::CAPACITY_SAMPLE) {
	    this->kindValue = requested;
	    this->reasonValue = reason ? reason : "";
	    decision.changed = true;
	} else if (this->kindValue == requested) {
	    /* Same-strength requests coalesce but retain the newest diagnostic. */
	    this->reasonValue = reason ? reason : "";
	}
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

    bool consume(SbString *reason, SbBool *capacityRelevant)
    {
	const bool wasPending = this->pending();
	if (reason)
	    *reason = this->reasonValue;
	if (capacityRelevant)
	    *capacityRelevant = this->capacityRelevant() ? TRUE : FALSE;
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

    const SbString &reason(void) const
    {
	return this->reasonValue;
    }

private:
    Kind kindValue = Kind::NONE;
    SbString reasonValue;
};

/* Uninstalled state shared by the responsibility-specific controller units.
 * It retains direct storage and one allocation: extracting scene forwarding
 * must not add virtual dispatch or accessors to LoD hot paths. */
struct BObolViewController::Impl : BObolLodCoordinator {
    explicit Impl(BObolViewController *owner) :
	headlightOffsetEye(bobol_headlight_default_offset()),
	featureStore(new BObolFeatureStore(owner)),
	polygonStore(new BObolPolygonStore(owner)),
	selectionStore(new BObolSelectionStore)
    {
    }

    /* Renderer-only presentation limits are mirrored in BObolViewLodState.
     * Keep publication in one owner operation: a controller value without
     * the matching renderer value is not a valid intermediate state. */
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
	if (state)
	    state->setCadPresentationProgressiveCutCeiling(
		ceiling, nextFraction);
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
	this->lodPresentationPointProxyPixelThreshold = threshold;
	BObolViewLodState *state = this->viewAttachment ?
	    this->viewAttachment->getViewLodState() : NULL;
	if (state)
	    state->setCadPresentationPointProxyPixelThreshold(threshold);
    }

    void publishCadDiscoveryPointProxyThreshold(float threshold)
    {
	threshold = normalizedCadPointProxyThreshold(threshold);
	this->lodDiscoveryPointProxyPixelThreshold = threshold;
	BObolViewLodState *state = this->viewAttachment ?
	    this->viewAttachment->getViewLodState() : NULL;
	if (state)
	    state->setCadPresentationDiscoveryPointProxyPixelThreshold(
		threshold);
    }

    void publishCadPointProxyThresholds(float presentationThreshold,
	float discoveryThreshold)
    {
	presentationThreshold = normalizedCadPointProxyThreshold(
	    presentationThreshold);
	discoveryThreshold = normalizedCadPointProxyThreshold(
	    discoveryThreshold);
	this->lodPresentationPointProxyPixelThreshold = presentationThreshold;
	this->lodDiscoveryPointProxyPixelThreshold = discoveryThreshold;
	BObolViewLodState *state = this->viewAttachment ?
	    this->viewAttachment->getViewLodState() : NULL;
	if (state) {
	    state->setCadPresentationPointProxyPixelThreshold(
		presentationThreshold);
	    state->setCadPresentationDiscoveryPointProxyPixelThreshold(
		discoveryThreshold);
	}
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

    static constexpr float POINT_PROXY_PIXEL_THRESHOLD_MINIMUM = 1.0f;
    static constexpr float POINT_PROXY_PIXEL_THRESHOLD_MAXIMUM = 64.0f;

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
    BObolViewController::SoftwareWireMode softwareWireMode =
	BObolViewController::SOFTWARE_WIRE_AUTO;
    std::atomic<int> endpointGraphicalRenderingEnabled {1};
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
    uint64_t hostWorkRevision = 0;
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

#endif /* LIBBOBOL_VIEW_CONTROLLER_PRIVATE_H */
