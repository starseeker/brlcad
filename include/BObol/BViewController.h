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
    /** Admit the raster-stable (at most quarter-pixel error) PoP cut without
     * interactive frame-time/quiet-time ceilings.  The mode remains active
     * across controller-owned render pumps until an explicit options call
     * disables it.  Offline captures and deterministic tests still present
     * every selected cut before advancing, but need not emulate human pauses
     * or depend on host rendering speed. */
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
    BOBOL_LOD_CONVERGENCE_ERROR,
    /* Keep established diagnostic values stable when adding HUD phases. */
    BOBOL_LOD_CONVERGENCE_PREPARING
};

/** Internal coordinator phase at the last owner-thread transition boundary.
 *
 * Unlike BObolLodConvergencePhase, which is a user-facing progress category,
 * this reports the actual retained-display state machine. */
enum BObolLodCoordinatorPhase {
    BOBOL_LOD_COORDINATOR_FALLBACK = 0,
    BOBOL_LOD_COORDINATOR_COVERAGE,
    BOBOL_LOD_COORDINATOR_INTERACTIVE,
    BOBOL_LOD_COORDINATOR_SETTLING,
    BOBOL_LOD_COORDINATOR_STABLE,
    BOBOL_LOD_COORDINATOR_COMPACTING
};

/** Owner-thread event which most recently reconciled the LoD coordinator. */
enum BObolLodCoordinatorEvent {
    BOBOL_LOD_EVENT_INITIALIZE = 0,
    BOBOL_LOD_EVENT_FRAME_COMPLETED,
    BOBOL_LOD_EVENT_WORK_SCHEDULED,
    BOBOL_LOD_EVENT_WORK_PUMPED,
    BOBOL_LOD_EVENT_RESULT_PUBLISHED,
    BOBOL_LOD_EVENT_SERVICE_CHANGED,
    BOBOL_LOD_EVENT_GENERATION_CANCELLED,
    BOBOL_LOD_EVENT_AUTO_SUBMIT_CHANGED,
    BOBOL_LOD_EVENT_VIEW_INVALIDATED,
    BOBOL_LOD_EVENT_POLICY_CHANGED,
    BOBOL_LOD_EVENT_INTERACTION_STARTED,
    BOBOL_LOD_EVENT_INTERACTION_ENDED,
    BOBOL_LOD_EVENT_VIEW_OBSERVED
};

/** User-facing progress for one view epoch.
 *
 * The fraction is a cost-weighted estimate of progress toward the current
 * view's terminal, frame-rate-aware presentation.  Structural discovery,
 * initial useful representation, and detailed mesh resolution are distinct
 * work classes: one completed subpixel proxy must not hide a comparatively
 * expensive outstanding mesh.  A value of one never promises that every
 * source triangle is resident.  Cache writes are reported separately because
 * they may continue after the visible view is ready. */
struct BOBOL_EXPORT BObolLodConvergenceStatus {
    BObolLodConvergenceStatus(void);
    void clear(void);

    int phase;
    int coordinatorPhase;
    int coordinatorEvent;
    uint64_t coordinatorTransitionSerial;
    uint64_t coordinatorProgressSequence;
    uint64_t coordinatorDispatchSerial;
    uint64_t coordinatorStagnantDispatchCount;
    uint64_t coordinatorInvariantViolationCount;
    uint32_t coordinatorInvariantMask;
    uint32_t coordinatorInvariantHistoryMask;
    size_t viewQualityHistoryEntryCount;
    size_t viewQualityHistoryRememberCount;
    size_t viewQualityHistoryRecallCount;
    uint64_t viewRevision;
    uint64_t activeGeneration;
    size_t submissionSourceIndex;
    size_t submissionEntryOffset;
    size_t expectedLeafCount;
    size_t availableLeafCount;
    size_t visibleTargetCount;
    size_t activePayloadCount;
    size_t satisfiedPayloadCount;
    size_t presentedSubpixelOccurrenceCount;
    size_t presentedStructuralBoxCount;
    size_t terminalOccurrenceFailureCount;
    size_t pendingTasks;
    size_t inFlight;
    size_t queuedResults;
    size_t queuedCacheWrites;
    size_t activeFaces;
    size_t activeRenderCost;
    size_t renderCostBudget;
    size_t residentMeshBytes;
    size_t stableResidentMeshBytes;
    size_t reservedResidentMeshGrowthBytes;
    size_t residentMeshLimitBytes;
    size_t memoryLimitedPayloadCount;
    size_t activeWorkingSetBytes;
    size_t peakWorkingSetBytes;
    uint64_t residentCompactionCount;
    size_t gpuTrackedBufferBytes;
    size_t gpuOrdinaryPartBufferBytes;
    size_t gpuProgressiveCutBufferBytes;
    size_t gpuProgressiveActiveCutBytes;
    size_t gpuBatchBufferBytes;
    size_t gpuTriangleAtlasAllocatedBytes;
    size_t gpuTriangleAtlasLiveBytes;
    size_t gpuTriangleAtlasConfiguredCapacityBytes;
    size_t gpuTriangleAtlasPartCount;
    size_t gpuTriangleAtlasPageCount;
    uint64_t gpuOrdinaryPartFullUploadBytes;
    uint64_t gpuOrdinaryPartSuffixUploadBytes;
    uint64_t gpuOrdinaryPartGpuCopyBytes;
    uint64_t gpuOrdinaryPartLineageReuseCount;
    uint64_t gpuOrdinaryPartLineageReplacementCount;
    uint64_t gpuTriangleAtlasFullUploadBytes;
    uint64_t gpuTriangleAtlasSuffixUploadBytes;
    uint64_t gpuTriangleAtlasLineageReuseCount;
    size_t gpuPressureProxyCount;
    uint64_t gpuProgressiveEvictionCount;
    uint64_t gpuTriangleAtlasReclamationCount;
    uint64_t gpuResourceSampleSerial;
    float fraction;
    SbBool viewReady;
    SbBool backgroundPending;
    SbBool performanceLimited;
    SbBool memoryLimited;
    SbBool gpuMemoryPressure;
    /* Individual convergence obligations.  These are intentionally exposed
     * in the aggregate status so hosts and regression reports can distinguish
     * a pending PoP presentation frame from nonvisual scene reallocation. */
    SbBool refinementFramePending;
    SbBool budgetCalibrationPending;
    SbBool stablePresentationHandoffPending;
    SbBool pointProxyCalibrationPending;
    SbBool pointProxyAdmissionFramePending;
    SbBool stablePointProxyCalibrationPending;
    SbBool pointProxyTriangleRecoveryPending;
    SbBool residentGrowthReallocationPending;
    SbBool publicationFramePending;
    /* Foreground cold-start work which is still publishing the immutable
     * population against which terminal view allocation will be proved. */
    SbBool sourcePreparationPending;
    size_t sourcePreparationProviderCount;
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

/** Level-triggered work advertised by a view controller to its presentation
 * host.  Pumping and rendering are independent: background/provider progress
 * may need bounded owner-thread service without traversing an unchanged
 * scene, while a render request may be satisfied from already published
 * retained data. */
enum BObolHostWorkFlag {
    BOBOL_HOST_WORK_NONE = 0,
    BOBOL_HOST_WORK_PUMP = 1u << 0,
    BOBOL_HOST_WORK_RENDER = 1u << 1,
    BOBOL_HOST_WORK_CAPACITY_SAMPLE = 1u << 2
};

/** Immutable observation of the controller/host work boundary.  Revision is
 * advanced by every level transition; renderRevision identifies the render
 * transaction and prevents an older completed frame from clearing a request
 * published while that frame was in flight. */
struct BOBOL_EXPORT BObolHostWorkSnapshot {
    BObolHostWorkSnapshot(void);

    uint64_t revision;
    uint64_t renderRevision;
    uint32_t flags;

    SbBool pumpPending(void) const;
    SbBool renderPending(void) const;
    SbBool capacitySampleRequested(void) const;
};

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

    /** Named, renderer-independent camera-lighting policies.  Values match
     * enum bv_lighting_profile so GED can synchronize without translation
     * tables or client-specific defaults. */
    enum LightingProfile {
	LIGHTING_STUDIO = 0,
	LIGHTING_MGED = 1
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
    void setLightingProfile(LightingProfile profile);
    LightingProfile getLightingProfile(void) const;
    float getLightingAmbientIntensity(void) const;
    void setNormalStyle(BObolViewLodState::NormalStyle style,
	float creaseAngleDegrees = 60.0f);
    BObolViewLodState::NormalStyle getNormalStyle(void) const;
    float getNormalCreaseAngle(void) const;
    /** Enable/disable the camera-driven light rig (layered under the master
     * setLightingEnabled()).  The historical name remains the concise public
     * command vocabulary for this independent lighting source. */
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
    /** Snapshot enabled camera-rig lights in world space for retained
     * non-Obol renderers such as the librt preview. */
    void getCameraLights(
	std::vector<BObolSceneLightRealization> &lights) const;
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

    /** Request a frame whose completed CAD traversal is valid evidence for
     * LoD capacity and deadline control. */
    void requestRender(const char *reason = NULL);
    /** Request a presentation-only frame.  Selection, highlighting, HUD, and
     * other style changes must become visible, but their one-time patch cost
     * is not evidence that the retained geometry cut is unsustainable.  If a
     * capacity-relevant request is already pending, the combined frame remains
     * capacity relevant. */
    void requestPresentationRender(const char *reason = NULL);
    void setFrameRequestCallback(BObolFrameRequestCallback callback,
	void *userData);
    void clearFrameRequestCallback(void *userData);
    void setPresentationSyncCallback(BObolPresentationSyncCallback callback,
	void *userData);
    void clearPresentationSyncCallback(void *userData);
    void synchronizePresentation(void);
    /** Obtain the complete level-triggered host work contract in one atomic
     * observation.  Hosts should keep their bounded service loop scheduled
     * while either PUMP or RENDER remains asserted. */
    BObolHostWorkSnapshot getHostWorkSnapshot(void) const;
    void clearRenderRequest(void);
    SbBool consumeRenderRequest(SbString *reason = NULL,
	SbBool *lodCapacityRelevant = NULL);
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
    void completeRenderTiming(uint64_t startedNanoseconds,
	SbBool lodCapacityRelevant = TRUE);
    /** Bound a graphical traversal before it can monopolize the endpoint
     * thread.  Zero disables the corresponding deadline.  Interrupted
     * frames are not presentation samples; hosts preserve the last completed
     * image and call notePresentationRenderInterrupted(). */
    void setPresentationFrameDeadlines(uint64_t interactiveNanoseconds,
	uint64_t stableNanoseconds);
    uint64_t getInteractivePresentationFrameDeadline(void) const;
    uint64_t getStablePresentationFrameDeadline(void) const;
    uint64_t getCurrentPresentationFrameDeadline(void) const;
    /** TRUE only when the current retained population has nonzero
     * LoD-managed work which deadline recovery can make cheaper.  Hosts that
     * traverse the render root directly must use this together with the
     * capacity-relevance bit returned by consumeRenderRequest(). */
    SbBool isLodPresentationCapacityRelevant(void) const;
    void notePresentationRenderInterrupted(uint64_t elapsedNanoseconds,
	SbBool cadDrawAttempted = TRUE,
	SbBool cadPreparationChanged = FALSE,
	SbBool lodCapacityRelevant = TRUE);
    uint64_t getInterruptedPresentationFrameCount(void) const;
    uint64_t getLastInterruptedPresentationTimeNanoseconds(void) const;
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
    /** Monotonic count of completed host/offscreen presentations.  Hosts and
     * test drivers may use it to distinguish published scene/camera state
     * from a frame which has actually reached the presentation surface. */
    uint64_t getPresentedFrameSerial(void) const;
    /** Monotonic controller-side render-request mutation token.  This is
     * diagnostic/liveness state: it lets a host distinguish an unchanged
     * pending level from a newer request published while a prior frame was
     * being consumed. */
    uint64_t getRenderRequestSerial(void) const;
    /** Monotonic count of complete CAD presentation barriers.  Unlike
     * getPresentedFrameSerial(), this advances only for traversals accepted
     * by the LoD state machine (not retained images from interrupted frames). */
    uint64_t getRenderCompletionSerial(void) const;
    /** Completion serial required by the current camera-settle gate, or zero
     * when no camera frame is outstanding. */
    uint64_t getLodSettleAfterRenderSerial(void) const;
    /** Completion serial required before the next progressive cut may be
     * admitted, or zero when no refinement presentation is outstanding. */
    uint64_t getLodRefinementResumeAfterRenderSerial(void) const;
    /** Short-horizon presentation cadence used by diagnostic telemetry.  LoD
     * capacity uses measured CPU/GPU work instead: an event-driven host gap
     * is not a renderer cost. */
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
    /** Begin and claim a result generation on the configured service.
     * External producers must submit through this generation if this
     * controller is to consume their results; unscoped FIFO draining is
     * intentionally unsupported because it lets one view steal another's
     * work. */
    uint64_t beginLodGeneration(void);
    void setLodAutoSubmit(SbBool enabled);
    SbBool isLodAutoSubmitEnabled(void) const;
    void setLodForcedCut(int cut);
    void clearLodForcedCut(void);
    SbBool hasLodForcedCut(void) const;
    int getLodForcedCut(void) const;
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
			size_t maxEstimatedBytes = 0,
			uint64_t generation = 0);
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
     * render-cost budget from measured frames so shaded faces, wire segments,
     * points, and repeated occurrences compete for one measured resource. */
    void setLodFrameRateTargets(float interactiveFps, float stableFps);
    float getLodInteractiveTargetFps(void) const;
    float getLodStableTargetFps(void) const;
    size_t getCurrentLodRenderCostBudget(void) const;
    size_t getActiveLodFaceCount(void) const;
    size_t getActiveLodRenderCost(void) const;
    /** Return the calibration selected for the controller's current
     * interaction state. */
    double getCalibratedLodRenderCostPerSecond(void) const;
    double getInteractiveCalibratedLodRenderCostPerSecond(void) const;
    double getStableCalibratedLodRenderCostPerSecond(void) const;
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
			     const SbVec3f &boundsMax,
			     SbBool boundsExact = FALSE);
    int setDatabaseSourceInstanceBoundsState(const char *sourceInstanceKey,
	    SbBool boundsValid,
	    const SbVec3f &boundsMin,
	    const SbVec3f &boundsMax,
	    SbBool boundsExact = FALSE);
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
    void invalidateRendererPerformanceHistory(void);
    void notifyFrameRequest(const char *reason);
    void setViewportSceneGraphWithLod(SoNode *root);
    void cancelActiveLodGeneration(void);
    void resetDiscoveryPointProxyFloor(SbBool requestFrame);
    void invalidateDatabaseSourceLodState(void);
    void syncRenderManager(void);
    void advanceLodViewRevision(void);
    void advanceLodPolicyRevision(
	SbBool preserveScaleDemandRefresh = FALSE);
    void syncLodViewSignature(SbBool advanceOnChange = TRUE);
    void scheduleLodRefinementFrame(const char *reason);
    void completePresentationBarrier(uint64_t elapsedNanoseconds,
	size_t provenRenderCost = 0);
    void scheduleResidentGrowthReallocationIfReady(void);
    void armStableLodHeadroomProbeIfReady(void);
    void resumeLodAfterRetainedRecovery(void);
    size_t enforceMeshResidencyBudget(void);
    static void lodResultReadyCB(BObolLodService *service, void *userData);
    /* Rewrite the headlight direction from the last camera orientation so it
     * tracks the viewer (no-op unless the headlight is enabled and tracked). */
    void applyTrackedHeadlight(SbBool force = FALSE);

    struct Impl;
    std::unique_ptr<Impl> d;
};

#endif /* BOBOL_BVIEWCONTROLLER_H */
