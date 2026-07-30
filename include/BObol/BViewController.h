/*             B V I E W C O N T R O L L E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BViewController.h */

#ifndef BOBOL_BVIEWCONTROLLER_H
#define BOBOL_BVIEWCONTROLLER_H

#include "BObol/BDefines.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BPickDetail.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbRotation.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec2f.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SoDB.h>
#include <memory>
#include <stddef.h>
#include <stdint.h>
#include <vector>

#include "vmath.h"

class BObolLodService;
class BObolSceneController;
class BObolViewAttachment;
class BObolViewLodState;
class BObolFeatureStore;
class BObolPolygonStore;
class BObolSelectionStore;
class SoBRLExportAction;
class SoBRLMeasureAction;
class SoBRLSnapAction;
class SoBRLViewLodGroup;
class SoCamera;
class SoGroup;
class SoNode;
class SoOffscreenRenderer;
class SoRenderManager;
class SoViewport;
struct bg_line_layer_builder;
struct bobol_display_endpoint;
struct db_i;
struct bv_view_info;

class BObolViewController;

BOBOL_EXPORT SbMatrix bobol_sbmatrix_from_brl_mat(const mat_t mat);
BOBOL_EXPORT SbRotation bobol_camera_orientation_from_brl_rotation(
    const mat_t rotation);

struct BOBOL_EXPORT BObolProgressiveOptions {
    BObolProgressiveOptions(void);

    size_t maxLodResults;
    uint64_t maxLodApplyMicroseconds;
    size_t maxProviders;
    /** Maximum streamed database occurrences merged by one provider pump.
     * Zero removes the item limit, but maxProviderMicroseconds still applies
     * between merge batches. */
    size_t maxProviderItems;
    /** Cooperative host-thread provider budget.  A single merge batch is
     * atomic, so this is checked between bounded batches.  Zero disables the
     * time limit. */
    uint64_t maxProviderMicroseconds;
    /** Admit the requested pixel-exact PoP cut without the interactive
     * frame-time/quiet-time ceilings.  Offline captures and deterministic
     * tests still present every selected cut before advancing, but need not
     * spend eight seconds emulating a human pause. */
    SbBool forceTerminalLodRefinement;
};

struct BOBOL_EXPORT BObolProgressiveStatus {
    BObolProgressiveStatus(void);
    void clear(void);

    size_t providerCount;
    size_t providerAdvanced;
    size_t lodResultsProcessed;
    size_t lodResultsApplied;
    size_t submitted;
    size_t alreadyCached;
    size_t expanded;
    size_t existing;
    size_t remaining;
    size_t proxyPublished;
    size_t metadataApplied;
    size_t pendingTasks;
    size_t inFlight;
    size_t queuedResults;
    size_t queuedCacheWrites;
    int changed;
    int hasMore;
};

enum BObolLodConvergencePhase {
    BOBOL_LOD_CONVERGENCE_IDLE = 0,
    BOBOL_LOD_CONVERGENCE_DISCOVERING,
    BOBOL_LOD_CONVERGENCE_INTERACTIVE,
    BOBOL_LOD_CONVERGENCE_REFINING,
    BOBOL_LOD_CONVERGENCE_CALIBRATING,
    BOBOL_LOD_CONVERGENCE_BACKGROUND,
    BOBOL_LOD_CONVERGENCE_ERROR
};

/** User-facing progress for one view epoch.
 *
 * The fraction describes progress toward the current view's terminal,
 * frame-rate-aware presentation.  A value of one never promises that every
 * source triangle is resident.  Cache writes are reported separately because
 * they may continue after the visible view is ready. */
struct BOBOL_EXPORT BObolLodConvergenceStatus {
    BObolLodConvergenceStatus(void);
    void clear(void);

    int phase;
    uint64_t viewRevision;
    size_t expectedLeafCount;
    size_t availableLeafCount;
    size_t visibleTargetCount;
    size_t activePayloadCount;
    size_t satisfiedPayloadCount;
    size_t pendingTasks;
    size_t inFlight;
    size_t queuedResults;
    size_t queuedCacheWrites;
    size_t activeFaces;
    size_t faceBudget;
    size_t residentMeshBytes;
    size_t stableResidentMeshBytes;
    size_t reservedResidentMeshGrowthBytes;
    size_t residentMeshLimitBytes;
    size_t memoryLimitedPayloadCount;
    size_t activeWorkingSetBytes;
    size_t peakWorkingSetBytes;
    uint64_t residentCompactionCount;
    float fraction;
    SbBool viewReady;
    SbBool backgroundPending;
    SbBool performanceLimited;
    SbBool memoryLimited;
    /* Individual presentation barriers.  These are intentionally exposed in
     * the aggregate status so hosts and regression reports can distinguish a
     * pending PoP frame from scene-budget calibration or motion handoff. */
    SbBool refinementFramePending;
    SbBool budgetCalibrationPending;
    SbBool stablePresentationHandoffPending;
    SbBool pointProxyCalibrationPending;
    unsigned int failedSourceCount;
};

typedef int (*BObolProgressiveAdvanceCallback)(
    BObolViewController *controller,
    void *userData,
    const BObolProgressiveOptions *options,
    BObolProgressiveStatus *status);
typedef void (*BObolProgressiveUserDataFreeCallback)(void *userData);

/* A positive callback result records provider progress.  Set status->changed
 * only when that progress published a visible scene update; status->hasMore
 * independently requests a subsequent frame for pending work. */

typedef void (*BObolFrameRequestCallback)(void *userData,
    const char *reason);

/* Runs on the controller/host owner thread immediately before a scene is
 * rendered.  Producers may use it to apply thread-safe image-stream updates
 * to retained Coin presentation nodes without letting worker threads touch
 * Coin or an OpenGL context. */
typedef void (*BObolPresentationSyncCallback)(void *userData);

struct BOBOL_EXPORT BObolProgressiveProviderRecord {
    BObolProgressiveProviderRecord(void);

    uint64_t token;
    BObolProgressiveAdvanceCallback callback;
    void *userData;
    BObolProgressiveUserDataFreeCallback userDataFree;
};

/**
 * Narrow application-facing controller for an Obol-backed BRL-CAD view.
 *
 * The controller records view services that GUI and command layers need while
 * leaving scene hierarchy, camera fields, visibility, and renderable geometry
 * in Obol.  It is the migration target for code that currently carries a
 * legacy view/display-manager pair.
 */
class BOBOL_EXPORT BObolViewController
{
public:
    enum SoftwareWireMode {
	SOFTWARE_WIRE_AUTO = 0,
	SOFTWARE_WIRE_QUALITY = 1,
	SOFTWARE_WIRE_FAST = 2
    };

    BObolViewController(void);
    explicit BObolViewController(SoNode *root, SoCamera *camera = NULL);
    ~BObolViewController(void);

    void setSceneRoot(SoNode *root);
    SoNode *getSceneRoot(void) const;
    void setRenderSceneRoot(SoNode *root);
    SoNode *getRenderSceneRoot(void) const;
    SoNode *getRenderRoot(void) const;
    /** Stable retained framebuffer layers.  Underlay and overlay surround
     * the CAD render batch.  Interlay is inserted by the hosted GED render
     * composition between model geometry and view-local screen features. */
    SoGroup *getFramebufferUnderlayRoot(void) const;
    SoGroup *getFramebufferInterlayRoot(void) const;
    SoGroup *getFramebufferOverlayRoot(void) const;
    void setViewAttachment(BObolViewAttachment *attachment);
    BObolViewAttachment *getViewAttachment(void) const;
    BObolViewLodState *getViewLodState(void) const;
    void clearViewLodState(void);

    void setCamera(SoCamera *camera);
    SoCamera *getCamera(void) const;

    void setViewportRegion(const SbViewportRegion &region);
    const SbViewportRegion &getViewportRegion(void) const;
    void setViewportSize(unsigned int width, unsigned int height);
    void setBackgroundColors(const SbColor &bottom, const SbColor &top);
    const SbColor &getBackgroundBottomColor(void) const;
    const SbColor &getBackgroundTopColor(void) const;
    void setDepthTestEnabled(SbBool enabled);
    SbBool isDepthTestEnabled(void) const;
    void setLightingEnabled(SbBool enabled);
    SbBool isLightingEnabled(void) const;
    void setNormalStyle(BObolViewLodState::NormalStyle style,
	float creaseAngleDegrees = 60.0f);
    BObolViewLodState::NormalStyle getNormalStyle(void) const;
    float getNormalCreaseAngle(void) const;
    /** Enable/disable the camera-driven headlight (layered under the master
     * setLightingEnabled()). */
    void setHeadlightEnabled(SbBool enabled);
    SbBool isHeadlightEnabled(void) const;
    /** When TRUE the headlight direction tracks the camera each frame (old
     * main-branch style); when FALSE the direction stays scene-fixed. */
    void setHeadlightCameraTracked(SbBool tracked);
    SbBool isHeadlightCameraTracked(void) const;
    /** Eye-space headlight offset direction (normalized).  Not straight-on, to
     * avoid washed-out shading. */
    void setHeadlightOffset(const SbVec3f &eyeDir);
    SbVec3f getHeadlightOffset(void) const;
    /** Current world-space headlight travel direction (as last aimed). */
    SbVec3f getHeadlightDirection(void) const;
    /** Enable/disable in-scene (database-derived) light sources. */
    void setSceneLightsEnabled(SbBool enabled);
    SbBool isSceneLightsEnabled(void) const;
    /** Supply the in-scene lights (world-space) for this view.  The GED layer
     * derives these from the database's "light"-shader regions and pushes them
     * here (independent of the geometry realize/LoD cache). */
    void setSceneLights(const std::vector<BObolSceneLightRealization> &lights);
    /** Root group holding realized in-scene lights (NULL if none), so other
     * render paths (e.g. the software raytrace preview) can honor them. */
    SoNode *getSceneLightsRoot(void) const;
    /** Rebuild the in-scene light group from every database source's realized
     * light snapshots.  Safe to call after any realization path (the endpoint
     * realize, the render-time realize action, etc.). */
    void rebuildSceneLights(void);
    void setHeadlightColor(const SbColor &color);
    SbColor getHeadlightColor(void) const;
    void setHeadlightIntensity(float intensity);
    float getHeadlightIntensity(void) const;
    void setTransparencyEnabled(SbBool enabled);
    SbBool isTransparencyEnabled(void) const;
    /** Enable Obol's single-pass line and point smoothing.  This deliberately
     * does not enable expensive accumulation-buffer multipass rendering. */
    void setAntialiasingEnabled(SbBool enabled);
    SbBool isAntialiasingEnabled(void) const;
    SbBool setClipBounds(double minimum, double maximum);
    void getClipBounds(double &minimum, double &maximum) const;
    size_t getActiveClipPlanes(SbPlane planes[2]) const;
    void setDepthCueEnabled(SbBool enabled);
    SbBool isDepthCueEnabled(void) const;
    void renderBackground(void) const;
    void setSoftwareWireMode(SoftwareWireMode mode);
    SoftwareWireMode getSoftwareWireMode(void) const;
    SbBool syncCameraFromViewContext(const void *viewCtx,
				     SbBool createCamera = TRUE,
				     SbBool *changedOut = NULL);
    SbBool getViewInfo(struct bv_view_info *info) const;

    SbBool realizePending(void);

    /* Derived display mode: when the view's LoD policy has neither mesh nor CSG
     * LoD enabled, drawing behaves like the classic force-realize (whole tree
     * before first frame) path; when LoD is enabled the render paths stay on the
     * progressive coarse-first pipeline.  Consulted by renderPending() and the
     * headless host so the single `view lod` setting controls both without a
     * separate progressive-display toggle. */
    SbBool isForceRealizeDisplay(void) const;

    unsigned int getLastVisitedSourceCount(void) const;
    unsigned int getLastRealizedSourceCount(void) const;
    unsigned int getLastFailedSourceCount(void) const;
    const SbString &getLastDiagnostics(void) const;

    void requestRender(const char *reason = NULL);
    void setFrameRequestCallback(BObolFrameRequestCallback callback,
	void *userData);
    void clearFrameRequestCallback(void *userData);
    void setPresentationSyncCallback(BObolPresentationSyncCallback callback,
	void *userData);
    void clearPresentationSyncCallback(void *userData);
    void synchronizePresentation(void);
    void clearRenderRequest(void);
    SbBool consumeRenderRequest(SbString *reason = NULL);
    /** Snapshot/retire protocol for hosts which traverse the render root
     * directly instead of using renderPending().  Capture the serial after
     * all state synchronized into a frame, then clear it after presentation
     * only if no newer request was published while rendering. */
    uint64_t renderRequestSerialGet(void) const;
    void clearRenderRequestIfUnchanged(uint64_t serial);
    SbBool renderPending(SbBool clearWindow = TRUE,
			 SbBool clearZBuffer = TRUE,
			 SbString *reason = NULL);
    uint64_t beginRenderTiming(void) const;
    void completeRenderTiming(uint64_t startedNanoseconds);
    uint64_t getLastRenderTimeNanoseconds(void) const;
    uint64_t getSmoothedRenderTimeNanoseconds(void) const;
    /** Host-side phase timings for the most recently completed render. */
    uint64_t getLastBackgroundRenderTimeNanoseconds(void) const;
    uint64_t getLastSceneRenderTimeNanoseconds(void) const;
    /** Host-thread time spent preparing the most recent progressive frame.
     * These diagnostics deliberately exclude the GL traversal reported by
     * getLastRenderTimeNanoseconds(), making event-loop stalls attributable
     * without a sampling profiler. */
    uint64_t getLastProgressiveAdvanceTimeNanoseconds(void) const;
    uint64_t getLastLodResultProcessingTimeNanoseconds(void) const;
    uint64_t getLastProgressiveProviderTimeNanoseconds(void) const;
    uint64_t getLastLodSubmissionTimeNanoseconds(void) const;
    uint64_t getLastPresentationSyncTimeNanoseconds(void) const;
    /** Record one completed host/offscreen presentation.  This cadence is
     * intentionally separate from render work duration. */
    void noteFramePresented(void);
    /** Short-horizon presentation cadence used by adaptive rendering policy
     * and diagnostic telemetry.  This intentionally reacts faster than the
     * user-facing FPS indication. */
    uint64_t getSmoothedPresentationIntervalNanoseconds(void) const;
    uint64_t getSmoothedInteractivePresentationIntervalNanoseconds(void) const;
    /** Human-facing presentation cadence.  This uses an elapsed-time EMA so
     * its response is independent of frame rate and short scheduling bursts
     * do not make the FPS faceplate flicker. */
    uint64_t getDisplayedPresentationIntervalNanoseconds(void) const;
    /** Capture with the controller-bound provider, or with a transient
     * explicit override whose lifetime needs to extend only through this call. */
    int renderToImage(unsigned char **image,
		      int flip = 0,
		      int alpha = 0,
		      const SbColor *background = NULL,
		      SoDB::ContextManager *contextManager = NULL,
		      BObolProgressiveStatus *progressiveStatus = NULL);
    /** True only when presentation state has explicitly requested a frame.
     * Progressive work is reported independently by
     * hasProgressiveWorkPending() and must not make this query true: a host
     * may need to service a provider without repainting an unchanged scene. */
    SbBool isRenderRequested(void) const;
    SbString getRenderReason(void) const;
    uint64_t registerProgressiveProvider(
	BObolProgressiveAdvanceCallback callback,
	void *userData,
	BObolProgressiveUserDataFreeCallback userDataFree = NULL);
    void unregisterProgressiveProvider(uint64_t token);
    void clearProgressiveProviders(void);
    void *findProgressiveProviderData(
	BObolProgressiveAdvanceCallback callback) const;
    uint64_t findProgressiveProviderToken(
	BObolProgressiveAdvanceCallback callback) const;
    SbBool hasProgressiveProviders(void) const;
    void setDefaultProgressiveOptions(
	const BObolProgressiveOptions *options);
    const BObolProgressiveOptions &getDefaultProgressiveOptions(void) const;
    int advanceProgressiveWork(
	const BObolProgressiveOptions *options = NULL,
	BObolProgressiveStatus *status = NULL);
    void markProgressiveWorkPending(void);
    void clearProgressiveWorkPending(void);
    SbBool hasProgressiveWorkPending(void) const;

    void setLodService(BObolLodService *service);
    BObolLodService *getLodService(void) const;
    BObolLodService *ensureManagedLodService(size_t workerCount);
    void stopManagedLodService(void);
    size_t getManagedLodWorkerCount(void) const;
    void setLodAutoSubmit(SbBool enabled);
    SbBool isLodAutoSubmitEnabled(void) const;
    void setLodForcedLevel(int level);
    void clearLodForcedLevel(void);
    SbBool hasLodForcedLevel(void) const;
    int getLodForcedLevel(void) const;
    void setExactFullDetailBudget(uint64_t maxFaceCount,
				  uint64_t maxPointCount);
    uint64_t getMaxExactFullDetailFaceCount(void) const;
    uint64_t getMaxExactFullDetailPointCount(void) const;
    int consumeExportSourceFullDetail(SoBRLExportAction &exportAction,
				      uint64_t generation = 0,
				      int *submittedRequestCount = NULL);
    int consumeMeasureSourceFullDetail(SoBRLMeasureAction &measureAction,
				       uint64_t generation = 0,
				       int *submittedRequestCount = NULL);
    int consumeSnapSourceFullDetail(SoBRLSnapAction &snapAction,
				    uint64_t generation = 0,
				    int *submittedRequestCount = NULL);
    int prepareRtPickCaches(void);
    int getRtPickCacheCount(void) const;
    BObolRtPickCache *getRtPickCache(int index) const;
    uint32_t getRtPickCacheSourceRevision(int index) const;
    int pickSourceMeshExactRay(BObolSourceMeshPickResult &pick,
			       const SbVec3f &rayOrigin,
			       const SbVec3f &rayDirection,
			       uint64_t generation = 0,
			       int *submittedRequestCount = NULL);
    int pickRtExactRay(std::vector<BObolRtPickResult> &results,
		       const SbVec3f &rayOrigin,
		       const SbVec3f &rayDirection,
		       SbBool pickAll = FALSE);
    void clearRtPickCaches(void);
    void setMeshResidencyBudget(size_t maxResidentMeshBytes,
				SbBool evictDisplayPayloads = TRUE);
    void clearMeshResidencyBudget(void);
    SbBool hasMeshResidencyBudget(void) const;
    size_t getMaxResidentMeshBytes(void) const;
    SbBool isMeshResidencyDisplayEvictionEnabled(void) const;
    size_t evictMeshPayloadsToBudget(size_t maxBytes,
				     SbBool evictDisplayPayloads = TRUE);
    size_t getLastMeshBudgetInitialResidentBytes(void) const;
    size_t getLastMeshBudgetFinalResidentBytes(void) const;
    size_t getLastMeshBudgetFreedResidentBytes(void) const;
    size_t getLastMeshBudgetFreedFullDetailBytes(void) const;
    size_t getLastMeshBudgetFreedDisplayBytes(void) const;
    unsigned int getLastMeshBudgetVisitedMeshCount(void) const;
    unsigned int getLastMeshBudgetEvictedFullDetailMeshCount(void) const;
    unsigned int getLastMeshBudgetEvictedDisplayMeshCount(void) const;
    SbBool hasPendingLodResults(void) const;
    SbBool hasPendingLodSubmissions(void) const;
    /** True while a retained PoP cut or unchanged calibration probe must be
     * presented before another refinement step may be submitted. */
    SbBool hasPendingLodRefinementFrame(void) const;
    size_t processPendingLodResults(size_t maxResults = 0,
	uint64_t maxMicroseconds = 0);
    int submitLodRequestsIfNeeded(SbBool refreshMissing = TRUE,
				  int reset = 0);
    int submitLodRequests(BObolLodService *service = NULL,
			  uint64_t generation = 0,
			  SbBool refreshMissing = TRUE,
			  int reset = 0);
    int applyLodResults(BObolLodService *service = NULL,
			size_t maxResults = 0,
			size_t maxEstimatedBytes = 0);
    uint64_t getLodViewRevision(void) const;
    void setLodPolicyRevision(uint64_t revision);
    uint64_t getLodPolicyRevision(void) const;
    void beginLodInteraction(void);
    void endLodInteraction(void);
    SbBool isLodInteractionActive(void) const;
    SbBool isLodGestureActive(void) const;
    /** True when the newest camera revision changes projected scale (zoom,
     * viewport, or projection) rather than only camera pose.  Existing PoP
     * cuts may track scale changes immediately; pose-only interaction keeps
     * its current cuts unless measured frame pressure requires coarsening. */
    SbBool isLodScaleChangingInteraction(void) const;
    float getLodTargetPixelError(void) const;
    int getLodInteractiveProgressiveCeiling(void) const;
    /** Configure aggregate scene frame-rate goals.  Projected per-object
     * error remains the quality demand; these targets calibrate a total
     * displayed face/segment budget from measured frames so thousands of
     * individually reasonable cuts cannot collectively collapse FPS. */
    void setLodFrameRateTargets(float interactiveFps, float stableFps);
    float getLodInteractiveTargetFps(void) const;
    float getLodStableTargetFps(void) const;
    size_t getCurrentLodFaceBudget(void) const;
    size_t getActiveLodFaceCount(void) const;
    /** Return the calibration selected for the controller's current
     * interaction state. */
    double getCalibratedLodFacesPerSecond(void) const;
    double getInteractiveCalibratedLodFacesPerSecond(void) const;
    double getStableCalibratedLodFacesPerSecond(void) const;
    unsigned int getLastLodVisitedMeshCount(void) const;
    unsigned int getLastLodSubmittedTaskCount(void) const;
    unsigned int getLastLodUpdatedCutCount(void) const;
    unsigned int getLastLodSkippedMeshCount(void) const;
    size_t getLastLodResultCount(void) const;
    unsigned int getLastLodMatchedResultCount(void) const;
    unsigned int getLastLodAppliedResultCount(void) const;
    unsigned int getLastLodRejectedResultCount(void) const;
    unsigned int getLastLodUnmatchedResultCount(void) const;
    const SbString &getLastLodDiagnostics(void) const;
    size_t getActiveLodMeshPayloadCount(void) const;
    size_t getActiveLodProxyPayloadCount(int proxyKind) const;
    size_t getActiveLodCadPayloadCount(void) const;
    void getLodConvergenceStatus(
	BObolLodConvergenceStatus &status) const;

    BObolSceneController *getSceneController(void);
    const BObolSceneController *getSceneController(void) const;
    BObolFeatureStore &features(void);
    const BObolFeatureStore &features(void) const;
    BObolPolygonStore &polygons(void);
    const BObolPolygonStore &polygons(void) const;
    BObolSelectionStore &selection(void);
    const BObolSelectionStore &selection(void) const;
    SoViewport *getViewport(void);
    const SoViewport *getViewport(void) const;
    /** Bind the concrete rendering provider used by direct drawing and by
     * renderToImage() when it has no explicit per-call override.  A NULL
     * provider leaves scene/action services available but disables rendering;
     * BRL-CAD never consults SoDB's process-global provider as a fallback. */
    void setRenderContextManager(SoDB::ContextManager *manager);
    SoDB::ContextManager *getRenderContextManager(void) const;
    SoRenderManager *getRenderManager(void);
    const SoRenderManager *getRenderManager(void) const;

    int replaceLineLayerOverlay(const char *overlayId,
				const struct bg_line_layer_builder *builder,
				uint32_t sourceId = 0,
				SbBool selectable = TRUE);
    int replaceHUDLabelOverlay(const char *labelId,
			       const char *text,
			       const SbVec2f &position,
			       const SbColor &color,
			       float fontSize = 12.0f,
			       uint32_t sourceId = 0);
    int removeHUDLabelOverlay(const char *labelId);
    int replaceEditPreview(const char *previewId,
			   const char *identity,
			   const SbVec3f *points,
			   const int32_t *commands,
			   int count,
			   uint32_t sourceRevision = 0,
			   uint32_t inputsRevision = 0);
    int replaceEditPreviewWithIntent(const char *previewId,
				     const char *identity,
				     const char *editIntentId,
				     const char *editIntentRole,
				     const SbVec3f *points,
				     const int32_t *commands,
				     int count,
				     uint32_t sourceRevision = 0,
				     uint32_t inputsRevision = 0);
    int removeEditPreview(const char *previewId);

    SoGroup *findGroup(const char *groupPath) const;
    SoGroup *ensureGroup(const char *groupPath);
    int setGroupDrawIntent(const char *groupPath,
			   const char *intentPath,
			   int drawMode,
			   int fallbackDrawMode,
			   SbBool overlayIntent,
			   uint32_t revalidationRevision);
    int setGroupDisplayState(const char *groupPath,
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
			     uint32_t materialRevision);
    int renameGroup(const char *groupPath, const char *newLeafName);
    int appendChildToGroup(const char *groupPath, SoNode *child);
    int removeChildFromGroup(const char *groupPath, SoNode *child);
    int eraseGroupSubpath(const char *parentGroupPath,
			  const char *subpath);
    int removeGroup(const char *groupPath);
    int clearGroup(const char *groupPath);
    int getGroupChildCount(const char *groupPath) const;
    int getGroupDescendantGroupCount(const char *groupPath) const;
    int getGroupDatabaseSourceCount(const char *groupPath) const;

    SoNode *findShape(const char *shapePath) const;
    SoGroup *findShapeParent(const char *shapePath) const;
    int moveShapeToGroup(const char *shapePath, const char *groupPath);
    int removeShape(const char *shapePath);
    int setShapeDrawState(const char *shapePath,
			  int drawMode,
			  SbBool databaseIntent,
			  SbBool overlayIntent,
			  SbBool hudIntent);
    int setShapeDisplayState(const char *shapePath,
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
			     uint32_t materialRevision);
    int setShapeSourceState(const char *shapePath,
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
			    uint32_t ownerStaleReason);
    int setShapePlacementState(const char *shapePath,
			       SbBool drawMatrixValid,
			       const SbMatrix &drawMatrix,
			       SbBool drawCenterValid,
			       const SbVec3f &drawCenter,
			       SbBool drawSizeValid,
			       float drawSize);

    int replaceDatabaseSource(const char *sourcePath,
			      struct db_i *dbip,
			      int drawMode = SoBRLDatabaseSource::WIREFRAME,
			      uint32_t sourceRevision = 0);
    int replaceDatabaseSourceInstance(const char *sourceInstanceKey,
				      const char *sourcePath,
				      struct db_i *dbip,
				      int drawMode = SoBRLDatabaseSource::WIREFRAME,
				      uint32_t sourceRevision = 0);
    int setDatabaseSourceState(const char *sourcePath,
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
			       uint32_t materialRevision);
    int setDatabaseSourceInstanceState(const char *sourceInstanceKey,
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
				       uint32_t materialRevision);
    int setDatabaseSourceDisplayPatch(const char *sourcePath,
				      const BObolDatabaseSourceDisplayPatch &patch);
    int setDatabaseSourceInstanceDisplayPatch(const char *sourceInstanceKey,
	    const BObolDatabaseSourceDisplayPatch &patch);
    int setDatabaseSourceDisplayName(const char *sourcePath,
				     const char *displayName);
    int setDatabaseSourceInstanceDisplayName(const char *sourceInstanceKey,
	    const char *displayName);
    int setDatabaseSourceBoundsState(const char *sourcePath,
				     SbBool boundsValid,
				     const SbVec3f &boundsMin,
				     const SbVec3f &boundsMax);
    int setDatabaseSourceInstanceBoundsState(const char *sourceInstanceKey,
	    SbBool boundsValid,
	    const SbVec3f &boundsMin,
	    const SbVec3f &boundsMax);
    int setDatabaseSourceMaterialPolicy(const char *sourcePath,
					int materialPolicy);
    int setDatabaseSourceInstanceMaterialPolicy(const char *sourceInstanceKey,
	    int materialPolicy);
    int markDatabaseSourceStale(const char *sourcePath,
				uint32_t staleReason);
    int markDatabaseSourceInstanceStale(const char *sourceInstanceKey,
					uint32_t staleReason);
    int moveDatabaseSourceToGroup(const char *sourcePath,
				  const char *groupPath);
    int moveDatabaseSourceInstanceToGroup(const char *sourceInstanceKey,
					  const char *groupPath);
    int removeDatabaseSource(const char *sourcePath);
    int removeDatabaseSourceInstance(const char *sourceInstanceKey);
    int clearDatabaseSources(void);
    SoBRLDatabaseSource *getDatabaseSource(int index) const;
    int getDatabaseSourceCount(void) const;
    /** Return the authoritative sources under the active render scene.  This
     * includes a shared GED render root when it differs from the controller's
     * locally managed scene, without transferring source ownership.  Call only
     * on the Coin/controller owner thread. */
    std::vector<SoBRLDatabaseSource *> getRenderDatabaseSources(void) const;
    SoBRLDatabaseSource *findDatabaseSourceInstance(
	const char *sourceInstanceKey) const;
    SbBool getDatabaseSourceSummary(int index,
				    BObolDatabaseSourceSummary &summary) const;

private:
    friend struct bobol_display_endpoint;

    /* Display endpoints own renderer selection.  These private hooks let an
     * endpoint retain this controller's scene while preventing every host
     * path (including direct toolkit repaints) from treating NONE or
     * DIAGNOSTIC as a graphical renderer. */
    void setEndpointGraphicalRenderingEnabled(SbBool enabled);
    void notifyFrameRequest(const char *reason);
    void setViewportSceneGraphWithLod(SoNode *root);
    void cancelActiveLodGeneration(void);
    void invalidateDatabaseSourceLodState(void);
    void syncRenderManager(void);
    void advanceLodViewRevision(void);
    void advanceLodPolicyRevision(void);
    void syncLodViewSignature(SbBool advanceOnChange = TRUE);
    void scheduleLodRefinementFrame(const char *reason);
    size_t enforceMeshResidencyBudget(void);
    static void lodResultReadyCB(BObolLodService *service, void *userData);
    /* Rewrite the headlight direction from the last camera orientation so it
     * tracks the viewer (no-op unless the headlight is enabled and tracked). */
    void applyTrackedHeadlight(void);

    struct Impl;
    std::unique_ptr<Impl> d;
};

#endif /* BOBOL_BVIEWCONTROLLER_H */
