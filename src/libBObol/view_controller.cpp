/*                 V I E W _ C O N T R O L L E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/log.h"
#include "bu/str.h"
#include "bu/time.h"

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
#include <string.h>
#include <type_traits>
#include <utility>
#include <vector>

#include <Inventor/SbName.h>
#include <Inventor/SbPlane.h>
#include <Inventor/SoDB.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>
#include <Inventor/SoOffscreenRenderer.h>
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

/* Build the large-scene retained frontier without projecting or sorting every
 * occurrence on the GUI thread.  requestedLevel is the last exact projected
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
	const int activeLevel = payload ? payload->activeLevel : -1;
	const int requestedLevel = payload ? payload->requestedLevel : 15;
	const int errorBucket = std::max(0, std::min(15,
	    requestedLevel - std::max(0, activeLevel)));

	/* Secondary marginal value uses retained hierarchy counts but no mesh
	 * scan.  The requested level also supplies a bounded approximation of
	 * projected footprint for ties at the same worst-error tier. */
	double value = static_cast<double>(errorBucket + 1);
	if (payload && payload->progressiveMesh && activeLevel >= 0 &&
	    requestedLevel > activeLevel) {
	    const size_t activeFaces = payload->progressiveMesh->
		hierarchyFaceCount(activeLevel);
	    int nextLevel = requestedLevel;
	    for (int level = activeLevel + 1; level <= requestedLevel;
		 ++level) {
		if (payload->progressiveMesh->hierarchyFaceCount(level) >
			activeFaces) {
		    nextLevel = level;
		    break;
		}
	    }
	    const size_t nextFaces = payload->progressiveMesh->
		hierarchyFaceCount(nextLevel);
	    const size_t addedFaces = nextFaces > activeFaces ?
		nextFaces - activeFaces : 1;
	    const double footprint = std::ldexp(1.0,
		std::max(0, std::min(15, requestedLevel)));
	    const double currentError = std::ldexp(1.0,
		std::max(0, std::min(15, requestedLevel - activeLevel)));
	    const double nextError = std::ldexp(1.0,
		std::max(0, std::min(15, requestedLevel - nextLevel)));
	    value = footprint * std::max(0.0, currentError - nextError) /
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
    this->coordinatorTransitionSerial = 0;
    this->coordinatorProgressSequence = 0;
    this->viewRevision = 0;
    this->expectedLeafCount = 0;
    this->availableLeafCount = 0;
    this->visibleTargetCount = 0;
    this->activePayloadCount = 0;
    this->satisfiedPayloadCount = 0;
    this->pendingTasks = 0;
    this->inFlight = 0;
    this->queuedResults = 0;
    this->queuedCacheWrites = 0;
    this->activeFaces = 0;
    this->faceBudget = 0;
    this->residentMeshBytes = 0;
    this->stableResidentMeshBytes = 0;
    this->reservedResidentMeshGrowthBytes = 0;
    this->residentMeshLimitBytes = 0;
    this->memoryLimitedPayloadCount = 0;
    this->activeWorkingSetBytes = 0;
    this->peakWorkingSetBytes = 0;
    this->residentCompactionCount = 0;
    this->fraction = 0.0f;
    this->viewReady = FALSE;
    this->backgroundPending = FALSE;
    this->performanceLimited = FALSE;
    this->memoryLimited = FALSE;
    this->refinementFramePending = FALSE;
    this->budgetCalibrationPending = FALSE;
    this->stablePresentationHandoffPending = FALSE;
    this->pointProxyCalibrationPending = FALSE;
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

/* Default headlight direction expressed in eye space (camera looks down -Z,
 * +X right, +Y up).  Match BRL-CAD's established shaded-display convention:
 * the light travels from the viewer into the scene.  An oblique direction can
 * make form easier to read, but on unoriented plate-mode BoTs it also
 * exaggerates triangle-to-triangle normal changes enough to look like missing
 * geometry.  Keep the compatibility direction as the default and let the
 * per-view lighting policy select an oblique direction when desired. */
static SbVec3f
bobol_headlight_default_offset(void)
{
    return SbVec3f(0.0f, 0.0f, -1.0f);
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

/* A fixed-function OpenGL light is transformed into eye space when its
 * position/direction is submitted, not when geometry is drawn.  Keep the
 * headlight field in world space for the shader/RT paths, but traverse its node
 * immediately after the camera so OSMesa/fixed GL applies the view transform
 * before freezing the light. */
static void
controller_place_headlight(SoViewport *viewport,
			    SoDirectionalLight *headlight)
{
    if (!viewport || !viewport->getRoot() || !headlight)
	return;

    SoSeparator *root = viewport->getRoot();
    const int oldIndex = root->findChild(headlight);
    if (oldIndex >= 0) {
	headlight->ref();
	root->removeChild(oldIndex);
    }

    int cameraIndex = -1;
    for (int i = 0; i < root->getNumChildren(); i++) {
	SoNode *child = root->getChild(i);
	if (child && child->isOfType(SoCamera::getClassTypeId())) {
	    cameraIndex = i;
	    break;
	}
    }
    root->insertChild(headlight, cameraIndex >= 0 ? cameraIndex + 1 :
	root->getNumChildren());
    if (oldIndex >= 0)
	headlight->unref();
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
	int hasClipMinimum = 0;
	int hasClipMaximum = 0;
	SoDirectionalLight *legacyHeadlight = NULL;
	for (int i = 0; i < renderEnvironment->getNumChildren(); i++) {
	    SoNode *child = renderEnvironment->getChild(i);
	    if (child && child->isOfType(SoDepthBuffer::getClassTypeId())) {
		hasDepthBuffer = 1;
	    }
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
	    controller_place_headlight(viewport, legacyHeadlight);
	    legacyHeadlight->unref();
	}
	const int index = root->findChild(renderEnvironment);
	if (index > 0) {
	    renderEnvironment->ref();
	    root->removeChild(renderEnvironment);
	    root->insertChild(renderEnvironment, 0);
	    renderEnvironment->unref();
	}
	return;
    }

    renderEnvironment = new SoGroup;
    renderEnvironment->setName(SbName(controller_render_environment_name()));

    SoDepthBuffer *depthBuffer = new SoDepthBuffer;
    depthBuffer->test = TRUE;
    depthBuffer->write = TRUE;
    renderEnvironment->addChild(depthBuffer);

    SoEnvironment *environment = new SoEnvironment;
    environment->ambientIntensity = 0.3f;
    environment->ambientColor = SbColor(1.0f, 1.0f, 1.0f);
    renderEnvironment->addChild(environment);

    SoLightModel *lightModel = new SoLightModel;
    lightModel->model = SoLightModel::PHONG;
    renderEnvironment->addChild(lightModel);

    SoDirectionalLight *headlight = new SoDirectionalLight;
    headlight->setName(SbName(controller_headlight_name()));
    headlight->color = SbColor(1.0f, 1.0f, 1.0f);
    headlight->intensity = 1.0f;
    headlight->direction = SbVec3f(0.0f, 0.0f, -1.0f);

    SoClipPlane *clipMinimum = new SoClipPlane;
    clipMinimum->setName(SbName(controller_clip_plane_name(TRUE)));
    clipMinimum->on = FALSE;
    renderEnvironment->addChild(clipMinimum);

    SoClipPlane *clipMaximum = new SoClipPlane;
    clipMaximum->setName(SbName(controller_clip_plane_name(FALSE)));
    clipMaximum->on = FALSE;
    renderEnvironment->addChild(clipMaximum);

    root->insertChild(renderEnvironment, 0);
    controller_place_headlight(viewport, headlight);
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
    SoSeparator *root = viewport->getRoot();
    for (int i = 0; i < root->getNumChildren(); i++) {
	SoNode *child = root->getChild(i);
	if (child && child->isOfType(SoDirectionalLight::getClassTypeId()) &&
	    bu_strcmp(child->getName().getString(),
		controller_headlight_name()) == 0)
	    return static_cast<SoDirectionalLight *>(child);
    }
    return NULL;
}

static const char *
controller_scene_lights_group_name(void)
{
    return "BObolSceneLights";
}

/* Locate (creating if needed) the in-scene lights group, always positioned
 * after the camera and headlight in the viewport root so fixed-function GL
 * transforms their world-space positions/directions into eye space when the
 * light nodes are traversed. */
static SoGroup *
controller_scene_lights_group(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return NULL;
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
    SoDirectionalLight *headlight = controller_headlight(viewport);
    const int headlightIndex = headlight ? root->findChild(headlight) : -1;
    if (headlightIndex >= insertAt)
	insertAt = headlightIndex + 1;
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

struct BObolLodCoordinator : BObolLodPhaseTracker {
    Phase derivedPhase(void) const
    {
	Inputs inputs;
	inputs.interactive = lodInteractive || lodGestureActive;
	inputs.compacting = lodResidentCompactionPlanning;
	inputs.coverageActive = lodCoveragePassActive;
	inputs.coverageComplete = lodVisibleCoverageComplete;
	inputs.generationActive = lodActiveGeneration != 0;
	inputs.settlingWork =
	    lodSubmissionPending || lodSubmissionRescanPending ||
	    lodRetainedRefinementPending || lodBudgetRescanAfterFrame ||
	    lodStablePresentationHandoffActive ||
	    lodRefinementAwaitingFrame ||
	    lodResultPublicationAwaitingFrame ||
	    lodStableCalibrationFramesRemaining ||
	    lodStablePointProxyCalibrationPending;
	return phaseFor(inputs);
    }

    void reconcilePhase(size_t pendingCount = 0)
    {
	Witness witness;
	witness.viewEpoch = lodViewRevision.value();
	witness.policyEpoch = lodPolicyRevision.value();
	witness.renderSerial = renderCompletionSerial;
	witness.activeGeneration = lodActiveGeneration;
	witness.residentDemandRevision =
	    lodResidentCompactionDemandRevision;
	witness.visibleCount = lodCoveragePassVisibleCount;
	witness.completedCount = lodCoveragePassCoveredCount;
	witness.pendingCount = pendingCount;

	switch (derivedPhase()) {
	    case Phase::FALLBACK:
		enterFallback(witness);
		break;
	    case Phase::COVERAGE:
		enterCoverage(witness);
		break;
	    case Phase::INTERACTIVE:
		enterInteractive(witness);
		break;
	    case Phase::SETTLING:
		enterSettling(witness);
		break;
	    case Phase::STABLE:
		enterStable(witness);
		break;
	    case Phase::COMPACTING:
		enterCompacting(witness);
		break;
	    case Phase::COUNT:
		break;
	}
    }

    uint64_t lastRenderTimeNanoseconds = 0;
    uint64_t smoothedRenderTimeNanoseconds = 0;
    uint64_t lastBackgroundRenderTimeNanoseconds = 0;
    uint64_t lastSceneRenderTimeNanoseconds = 0;
    /* Full CAD plan/atlas construction is a real latency problem, but not a
     * measurement of steady triangle throughput.  Only an unchanged prepared
     * replay may drive quiet-view capacity cuts. */
    SbBool lodLastRenderWasPreparedCadReplay = TRUE;
    uint64_t lodLastCadGpuSampleSerial = 0;
    uint64_t lodLastCadGpuTimeNanoseconds = 0;
    uint64_t lastProgressiveAdvanceTimeNanoseconds = 0;
    uint64_t lastLodResultProcessingTimeNanoseconds = 0;
    uint64_t lastProgressiveProviderTimeNanoseconds = 0;
    uint64_t lastLodSubmissionTimeNanoseconds = 0;
    uint64_t lastPresentationSyncTimeNanoseconds = 0;
    uint64_t renderCompletionSerial = 0;
    mutable std::mutex presentationTimingMutex;
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
    std::unique_ptr<BObolLodService> managedLodService;
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
    SoBRLDatabaseSource *lodSubmissionVisibleCountSource = NULL;
    size_t lodSubmissionVisibleCount = 0;
    /* Large compact scenes are consumed in bounded GUI-thread windows.  A
     * coverage pass suppresses upward refinement until one complete current-
     * view scan has offered every projected leaf its minimum drawable mesh.
     * This prevents early entries from exhausting the scene budget while
     * later entries remain boxes. */
    SbBool lodCoveragePassActive = TRUE;
    SbBool lodCoveragePassSawBoundedSource = FALSE;
    size_t lodCoveragePassVisibleCount = 0;
    size_t lodCoveragePassCoveredCount = 0;
    /*
     * The newest completed camera pass proved that every projected mesh
     * occurrence already has a retained mesh presentation.  While input is
     * active, richer worker results may then remain queued without exposing
     * boxes or missing newly visible geometry.
     */
    SbBool lodVisibleCoverageComplete = FALSE;
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
    size_t lodConvergenceCompleteVisibleCount = 0;
    SbBool lodConvergenceCompleteVisibleCountValid = FALSE;
    SbBool lodStableBudgetLimited = FALSE;
    mutable BObolLodViewEpoch lodConvergenceFractionViewRevision;
    mutable BObolLodPolicyEpoch lodConvergenceFractionPolicyRevision;
    mutable float lodConvergenceFractionFloor = 0.0f;
    void clearLodSubmissionPlan(void)
    {
	lodSubmissionPlanSource = NULL;
	lodSubmissionPlanEntries.clear();
	lodSubmissionPlanValid = FALSE;
	lodSubmissionVisibleCountSource = NULL;
	lodSubmissionVisibleCount = 0;
    }
    void clearLodConvergenceCandidates(void)
    {
	lodConvergenceCandidateCounts.clear();
	lodConvergenceCompleteVisibleCount = 0;
	lodConvergenceCompleteVisibleCountValid = FALSE;
    }
    void setLodConvergenceCompleteVisibleCount(size_t count)
    {
	lodConvergenceCompleteVisibleCount = count;
	lodConvergenceCompleteVisibleCountValid = TRUE;
    }
    void resetLodConvergenceFraction(void)
    {
	lodConvergenceFractionViewRevision.reset();
	lodConvergenceFractionPolicyRevision.reset();
	lodConvergenceFractionFloor = 0.0f;
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
	if (lodConvergenceCompleteVisibleCountValid)
	    return lodConvergenceCompleteVisibleCount;
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
    SbBool lodRetainedRefinementCutAdvanced = FALSE;
    SbBool lodRetainedRefinementBudgetBlocked = FALSE;
    SbBool lodBudgetRescanAfterFrame = FALSE;
    unsigned int lodStableCalibrationFramesRemaining = 0;
    SbBool lodPassAdmittedWork = FALSE;
    SbBool forceTerminalLodRefinement = FALSE;
    SbBool lodRefinementFaceBudgetInitialized = FALSE;
    size_t lodRefinementFaceBudgetRemaining = SIZE_MAX;
    SbBool lodRetainedAdmissionActive = FALSE;
    size_t lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
    int64_t lodRefinementBudgetRetryMicroseconds = 0;
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
    BObolLodViewEpoch lodViewRevision {1};
    BObolLodPolicyEpoch lodPolicyRevision {1};
    int64_t lodLastViewChangeMicroseconds = 0;
    SbBool lodInteractive = FALSE;
    SbBool lodGestureActive = FALSE;
    SbBool lodViewScaleChanging = FALSE;
    /* Keep the complete held-button presentation through button-release
     * debounce.  This preserves the measured aggregate presentation, not the
     * richer retained per-occurrence cuts which may be hidden by the
     * renderer-wide motion ceiling. */
    SbBool lodReleaseCutFloorActive = FALSE;
    /*
     * Motion may leave rich immutable occurrence cuts resident while a
     * renderer-only PoP ceiling and small-part aggregation cut provide the
     * responsive presentation.  On the quiet transition, keep those
     * presentation limits until a bounded retained-admission pass has
     * reconciled the occurrence cuts with measured stable capacity.  Clearing
     * them first can turn one release event into a multi-million-triangle
     * flattened rebuild on the GUI thread.
     */
    SbBool lodStablePresentationHandoffActive = FALSE;
    /* Renderer-only stable presentation retained across a camera gesture.
     * Occurrence cuts and resident arrays remain authoritative; this snapshot
     * restores their last measured presentation immediately after pose-only
     * motion instead of walking back through the motion ceiling. */
    SbBool lodStablePresentationSnapshotValid = FALSE;
    SbBool lodStablePresentationSnapshotRestored = FALSE;
    SbBool lodInteractionScaleChanged = FALSE;
    int lodStablePresentationProgressiveCeiling = -1;
    float lodStablePresentationPointProxyPixelThreshold = 1.0f;
    size_t lodStablePresentationActiveFaces = 0;
    uint64_t lodStablePresentationCadRevision = 0;
    uint64_t lodSettleAfterRenderSerial = 0;
    SbBool lodResidentCompactionPending = FALSE;
    SbBool lodResidentCompactionPlanning = FALSE;
    uint64_t lodResidentCompactionDemandRevision = 0;
    uint64_t lodResidentAdmissionRevision = 0;
    /* A resident prefix can grow without a camera change (notably when an
     * erased path is redrawn).  Reclamation therefore needs an independent
     * quiet-period deadline rather than borrowing the view-change timer. */
    int64_t lodResidentCompactionAfterMicroseconds = 0;
    SbBool lodRefinementAwaitingFrame = FALSE;
    uint64_t lodRefinementResumeAfterRenderSerial = 0;
    /* A result-publication frame may include one-time CPU/GPU buffer
     * construction.  It controls refinement pacing, but must not teach the
     * steady-state face-capacity estimator that every later frame costs the
     * same amount. */
    SbBool lodResultPublicationAwaitingFrame = FALSE;
    /*
     * Applied worker results are immutable CPU bindings and need not force
     * one expensive scene render per owner-thread drain quantum.  Retain a
     * bounded-latency publication batch while the worker stream is busy.
     * Any ordinary camera/input frame publishes the complete accumulated
     * batch and clears these counters.
     */
    size_t lodUnpresentedResultCount = 0;
    int64_t lodUnpresentedResultFirstMicroseconds = 0;
    int64_t lodResultPresentationIntervalMicroseconds(void) const
    {
	const uint64_t observedNanoseconds =
	    std::max(lastSceneRenderTimeNanoseconds,
		lastRenderTimeNanoseconds);
	const int64_t minimum = lodInteractive ? 33333 : 50000;
	const int64_t maximum = lodInteractive ? 100000 : 250000;
	if (!observedNanoseconds)
	    return minimum;
	const uint64_t observedMicroseconds =
	    observedNanoseconds / 1000ULL;
	const uint64_t scaled = observedMicroseconds >
		static_cast<uint64_t>(
		    std::numeric_limits<int64_t>::max() / 2) ?
	    static_cast<uint64_t>(maximum) : observedMicroseconds * 2ULL;
	return static_cast<int64_t>(std::max<uint64_t>(
	    static_cast<uint64_t>(minimum),
	    std::min<uint64_t>(static_cast<uint64_t>(maximum), scaled)));
    }
    SbBool lodUnpresentedResultsDue(int64_t now) const
    {
	if (!lodUnpresentedResultCount ||
	    lodUnpresentedResultFirstMicroseconds <= 0)
	    return FALSE;
	const int64_t elapsed =
	    now >= lodUnpresentedResultFirstMicroseconds ?
	    now - lodUnpresentedResultFirstMicroseconds :
	    lodResultPresentationIntervalMicroseconds();
	return elapsed >= lodResultPresentationIntervalMicroseconds() ?
	    TRUE : FALSE;
    }
    int64_t lodRefinementNotBeforeMicroseconds = 0;
    float lodTargetPixelError = 1.0f;
    float lodInteractiveTargetFps = 60.0f;
    float lodStableTargetFps = 20.0f;
    size_t lodSeedFaceBudget = 50000;
    size_t lodCurrentFaceBudget = 50000;
    long double lodInteractiveCalibratedFacesPerSecond = 0.0L;
    long double lodStableCalibratedFacesPerSecond = 0.0L;
    size_t lodStableCalibrationProbeActiveFaces = 0;
    unsigned int lodStableCalibrationProbeCount = 0;
    /* One direct measured backoff per camera epoch bounds convergence when
     * draw-call/instance overhead, rather than triangle count, dominates.
     * The resulting coherent cut is reported performance-limited; the next
     * view epoch starts from its retained timing and may back off again. */
    SbBool lodStableOverloadRecoveryPerformed = FALSE;
    size_t lodStableOverloadRecoveryActiveFaces = 0;
    int lodInteractiveProgressiveCeiling = -1;
    /*
     * Current renderer-side small-occurrence aggregation cut.  Interaction
     * raises it immediately when measured frames miss their target.  A quiet
     * performance-limited view may also raise it when every retained mesh is
     * already at its minimum prefix and that irreducible population still
     * exceeds the calibrated scene budget.
     */
    float lodPresentationPointProxyPixelThreshold = 1.0f;
    /*
     * A changed quiet aggregation threshold needs one measured presentation
     * before convergence is authoritative.  Keep this distinct from PoP
     * admission: no provider scan or retained-array mutation is required.
     */
    SbBool lodStablePointProxyCalibrationPending = FALSE;
    /* True only while a motion-to-stable handoff is cautiously lowering the
     * interaction aggregation threshold toward one pixel.  A slow probe ends
     * relaxation after restoring a safe threshold, avoiding oscillation
     * across a discrete many-occurrence boundary. */
    SbBool lodStablePointProxyRelaxationActive = FALSE;
    uint64_t lodInteractiveCeilingFeedbackRenderSerial = 0;
    SbBool lodUseForcedLevel = FALSE;
    int lodForcedLevel = 0;
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
    /* Camera-driven headlight state.  When headlightCameraTracked is TRUE the
     * headlight direction is rewritten from the camera orientation each sync so
     * the light follows the viewer (old main-branch behavior); when FALSE the
     * direction is left constant (scene-fixed).  headlightEnabled is the finer
     * per-view toggle layered under the master setLightingEnabled(). */
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
    if (this->d->managedLodService)
	this->d->managedLodService->stop();
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

/* Rewrite the headlight's world-space direction from the stored camera
 * orientation so the light tracks the viewer.  No-op unless the headlight is
 * enabled and in camera-tracked mode. */
void
BObolViewController::applyTrackedHeadlight(void)
{
    if (!this->d->headlightEnabled || !this->d->headlightCameraTracked)
	return;
    SoDirectionalLight *light = controller_headlight(this->d->viewport);
    if (!light)
	return;
    SbVec3f worldDir;
    this->d->lastCameraOrientation.multVec(this->d->headlightOffsetEye, worldDir);
    if (worldDir.normalize() > 0.0f && light->direction.getValue() != worldDir)
	light->direction = worldDir;
}

void
BObolViewController::setLightingEnabled(SbBool enabled)
{
    SoLightModel *model = controller_light_model(this->d->viewport);
    SoDirectionalLight *light = controller_headlight(this->d->viewport);
    if (!model || !light)
	return;
    const int requested = enabled ? SoLightModel::PHONG :
	SoLightModel::BASE_COLOR;
    /* Master lighting picks PHONG vs flat shading.  The headlight only actually
     * contributes when both master lighting and the per-view headlight toggle
     * are on, so restore the headlight to its finer-grained state rather than
     * unconditionally forcing it on. */
    const SbBool lightOn = enabled && this->d->headlightEnabled;
    if (model->model.getValue() == requested &&
	light->on.getValue() == lightOn)
	return;
    model->model = requested;
    light->on = lightOn;
    this->requestRender("lighting");
}

SbBool
BObolViewController::isLightingEnabled(void) const
{
    SoLightModel *model = controller_light_model(this->d->viewport);
    return model && model->model.getValue() == SoLightModel::PHONG;
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
    SoDirectionalLight *light = controller_headlight(this->d->viewport);
    if (light) {
	/* Only illuminate when master lighting (PHONG) is also on. */
	const SbBool lightOn = enabled && this->isLightingEnabled();
	if (light->on.getValue() != lightOn)
	    light->on = lightOn;
	if (lightOn)
	    this->applyTrackedHeadlight();
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
    this->applyTrackedHeadlight();
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
    this->d->renderReason = "";
    this->d->renderRequestSerial++;
}

SbBool
BObolViewController::consumeRenderRequest(SbString *reason)
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    const SbBool ret = this->d->renderRequested;
    if (reason)
	*reason = this->d->renderReason;
    this->d->renderRequested = FALSE;
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
    this->d->renderReason = "";
    this->d->renderRequestSerial++;
}

uint64_t
BObolViewController::renderRequestSerialGet(void) const
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    return this->d->renderRequestSerial;
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
    if (!this->consumeRenderRequest(&renderReasonCopy))
	return FALSE;

    if (reason)
	*reason = renderReasonCopy;

    const uint64_t started = this->beginRenderTiming();
    if (clearWindow)
	this->renderBackground();
    else if (clearZBuffer)
	glClear(GL_DEPTH_BUFFER_BIT);
    const uint64_t backgroundCompleted = this->beginRenderTiming();
    this->d->renderManager->render(static_cast<SbBool>(FALSE), static_cast<SbBool>(FALSE));
    const uint64_t sceneCompleted = this->beginRenderTiming();
    this->d->lastBackgroundRenderTimeNanoseconds =
	backgroundCompleted >= started ?
	    backgroundCompleted - started : 0;
    this->d->lastSceneRenderTimeNanoseconds =
	sceneCompleted >= backgroundCompleted ?
	    sceneCompleted - backgroundCompleted : 0;
    this->completeRenderTiming(started);
    return TRUE;
}

uint64_t
BObolViewController::beginRenderTiming(void) const
{
    return static_cast<uint64_t>(
	std::chrono::duration_cast<std::chrono::nanoseconds>(
	    std::chrono::steady_clock::now().time_since_epoch()).count());
}

static float controller_point_proxy_threshold(
    float currentThreshold, uint64_t renderNanoseconds,
    uint64_t presentationNanoseconds, float targetFps);

void
BObolViewController::completeRenderTiming(uint64_t startedNanoseconds)
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
    size_t calibrationFaces = calibrationState ?
	calibrationState->activeFaceCount() : 0;
    size_t presentedCadFaces = 0;
    const SbBool measuredCadPresentation = calibrationState &&
	calibrationState->lastCadPresentedTriangleCount(presentedCadFaces);
    if (measuredCadPresentation)
	calibrationFaces = presentedCadFaces;
    size_t gpuCadFaces = 0;
    uint64_t gpuCadNanoseconds = 0;
    uint64_t gpuCadSampleSerial = 0;
    const SbBool measuredGpuCadPresentation = calibrationState &&
	calibrationState->lastCadGpuMeasurement(
	    gpuCadFaces, gpuCadNanoseconds, gpuCadSampleSerial);
    const SbBool newGpuCadMeasurement =
	measuredGpuCadPresentation &&
	gpuCadSampleSerial != this->d->lodLastCadGpuSampleSerial;
    if (newGpuCadMeasurement) {
	this->d->lodLastCadGpuSampleSerial = gpuCadSampleSerial;
	this->d->lodLastCadGpuTimeNanoseconds = gpuCadNanoseconds;
	calibrationFaces = gpuCadFaces;
    }
    this->d->lodLastRenderWasPreparedCadReplay =
	!measuredCadPresentation ||
	calibrationState->lastCadPresentationUsedPreparedReplay();
    uint64_t calibrationElapsed =
	newGpuCadMeasurement ? gpuCadNanoseconds : elapsed;
    if (this->d->lodGestureActive && !newGpuCadMeasurement) {
	uint64_t presentationInterval = 0;
	{
	    std::lock_guard<std::mutex> lock(
		this->d->presentationTimingMutex);
	    presentationInterval =
		this->d->
		    smoothedInteractivePresentationIntervalNanoseconds;
	}
	/* System GL may queue work which is paid at presentation rather than
	 * inside Coin's traversal.  During an interactive burst the larger of
	 * render and observed presentation cost is the defensible capacity
	 * sample.  Stable event-driven gaps are intentionally excluded. */
	if (presentationInterval > calibrationElapsed &&
	    presentationInterval <= 250000000ULL)
	    calibrationElapsed = presentationInterval;
    }
    /* Result publication may allocate/upload newly decoded geometry during
     * this traversal.  That one-time cost is relevant to the cooldown below,
     * but not to the sustainable FPS of the retained cut.  Quiet,
     * unchanged calibration probes measure the latter.  Interactive samples
     * keep the full cost because upload latency is part of the user's motion
     * experience. */
    const SbBool transientStablePublication =
	this->d->lodResultPublicationAwaitingFrame &&
	!this->d->lodInteractive && !this->d->lodGestureActive;
    const SbBool transientStablePreparation =
	!this->d->lodLastRenderWasPreparedCadReplay &&
	!this->d->lodInteractive && !this->d->lodGestureActive;
    if (calibrationFaces > 0 && calibrationElapsed > 0 &&
	(newGpuCadMeasurement ||
	 (!transientStablePublication && !transientStablePreparation))) {
	const long double sample =
	    static_cast<long double>(calibrationFaces) * 1000000000.0L /
	    static_cast<long double>(calibrationElapsed);
	if (std::isfinite(sample) && sample > 0.0L) {
	    const SbBool interactiveSample =
		this->d->lodInteractive || this->d->lodGestureActive;
	    long double &calibration = interactiveSample ?
		this->d->lodInteractiveCalibratedFacesPerSecond :
		this->d->lodStableCalibratedFacesPerSecond;
	    if (calibration <= 0.0L) {
		calibration = sample;
	    } else if (interactiveSample && sample < calibration) {
		/* Missing an interaction target is immediately visible to the
		 * user.  Adopt the observed lower ceiling now; budget growth is
		 * intentionally much slower and therefore cannot oscillate back
		 * to a one-frame optimistic spike. */
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
		    this->d->lodStableTargetFps > 0.0f ?
		    1000000000.0L /
			static_cast<long double>(
			    this->d->lodStableTargetFps) : 0.0L;
		if (targetNanoseconds > 0.0L &&
		    static_cast<long double>(calibrationElapsed) >
			targetNanoseconds * 1.20L)
		    calibration = calibration * 0.80L + sample * 0.20L;
	    } else {
		const long double oldWeight = interactiveSample ? 0.98L : 0.90L;
		calibration =
		    calibration * oldWeight + sample * (1.0L - oldWeight);
	    }
	}
    }

    /*
     * Consume a quiet small-part aggregation probe without invoking another
     * O(N) provider/admission scan.  If the measured reusable presentation
     * still misses the stable target, ratchet the camera-local threshold and
     * request exactly one successor frame.  The immutable per-part geometry
     * and desired PoP cuts are untouched throughout.
    */
    if (this->d->lodStablePointProxyCalibrationPending) {
	const SbBool relaxationProbe =
	    this->d->lodStablePointProxyRelaxationActive;
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	this->d->lodStablePointProxyRelaxationActive = FALSE;
	const long double stableTargetNanoseconds =
	    this->d->lodStableTargetFps > 0.0f ?
		1000000000.0L /
		    static_cast<long double>(
			this->d->lodStableTargetFps) : 0.0L;
	const SbBool reusableStableSample =
	    !this->d->lodInteractive &&
	    !this->d->lodGestureActive &&
	    measuredCadPresentation &&
	    this->d->lodLastRenderWasPreparedCadReplay &&
	    !transientStablePublication &&
	    !transientStablePreparation &&
	    stableTargetNanoseconds > 0.0L &&
	    calibrationElapsed > 0;
	float nextThreshold =
	    this->d->lodPresentationPointProxyPixelThreshold;
	SbBool continueRelaxation = FALSE;
	if (reusableStableSample &&
	    static_cast<long double>(calibrationElapsed) >
		stableTargetNanoseconds * 1.10L) {
	    nextThreshold = controller_point_proxy_threshold(
		this->d->lodPresentationPointProxyPixelThreshold,
		calibrationElapsed, 0, this->d->lodStableTargetFps);
	    /* A discrete relaxation step was too aggressive.  Restore the
	     * measured safe side, then stop probing downward for this camera
	     * epoch once that correction has been presented. */
	    continueRelaxation = FALSE;
	} else if (reusableStableSample && relaxationProbe &&
	    this->d->lodPresentationPointProxyPixelThreshold > 1.01f &&
	    static_cast<long double>(calibrationElapsed) <
		stableTargetNanoseconds * 0.80L) {
	    /*
	     * Small projected occurrences grow approximately with screen area.
	     * Move no more than 25 percent toward a finer threshold per measured
	     * frame, and use observed headroom for a still smaller step near the
	     * target.  This cannot expose the complete aggregated population in
	     * one frame.
	     */
	    const long double ratio = std::sqrt(
		static_cast<long double>(calibrationElapsed) /
		(stableTargetNanoseconds * 0.80L));
	    const long double factor =
		std::max<long double>(0.75L,
		    std::min<long double>(1.0L, ratio));
	    nextThreshold = static_cast<float>(std::max<long double>(
		1.0L,
		static_cast<long double>(
		    this->d->lodPresentationPointProxyPixelThreshold) *
		    factor));
	    continueRelaxation = nextThreshold > 1.01f ? TRUE : FALSE;
	}
	const SbBool thresholdChanged =
	    fabsf(nextThreshold -
		this->d->lodPresentationPointProxyPixelThreshold) > 0.01f ?
		TRUE : FALSE;
	if (thresholdChanged) {
	    this->d->lodPresentationPointProxyPixelThreshold =
		nextThreshold;
	    BObolViewLodState *presentationState =
		this->d->viewAttachment->getViewLodState();
	    if (presentationState)
		presentationState->
		    setCadPresentationPointProxyPixelThreshold(
			nextThreshold);
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	    this->d->lodStablePointProxyRelaxationActive =
		continueRelaxation;
	    this->requestRender("lod-stable-point-calibration");
	    this->notifyFrameRequest(
		"lod-stable-point-calibration");
	} else if (!reusableStableSample &&
	    !this->d->lodInteractive &&
	    !this->d->lodGestureActive &&
	    measuredCadPresentation &&
	    stableTargetNanoseconds > 0.0L &&
	    (transientStablePublication ||
	     transientStablePreparation)) {
	    /*
	     * Applying a different point threshold changes the retained
	     * occurrence classification.  Its first frame may rebuild/patch the
	     * prepared plan and is intentionally excluded from throughput
	     * calibration above, but that is a one-frame handoff—not evidence
	     * that calibration has converged.  Preserve the probe and request
	     * one unchanged replay which can measure the cut it actually drew.
	     *
	     * Result streaming may extend this handoff across several frames.
	     * Those frames already represent real pending work; retaining this
	     * latch neither busy-spins nor exposes a richer cut prematurely.
	     */
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	    this->d->lodStablePointProxyRelaxationActive =
		relaxationProbe;
	    this->requestRender("lod-stable-point-replay");
	    this->notifyFrameRequest("lod-stable-point-replay");
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
	this->d->lodStableTargetFps > 0.0f ?
	    1000000000.0L /
		static_cast<long double>(
		    this->d->lodStableTargetFps) : 0.0L;
    const SbBool ongoingStableWork =
	this->d->lodSubmissionPending ||
	this->hasPendingLodResults() ||
	(this->d->lodService &&
	 (this->d->lodService->pendingTaskCountForDiagnostics() > 0 ||
	  this->d->lodService->inFlightCount() > 0 ||
	  this->d->lodService->queuedResultCountForDiagnostics() > 0));
    if (!this->d->lodStablePointProxyCalibrationPending &&
	!this->d->lodInteractive &&
	!this->d->lodGestureActive &&
	this->d->lodAutoSubmit &&
	ongoingStableWork &&
	measuredCadPresentation &&
	ongoingStableTargetNanoseconds > 0.0L &&
	static_cast<long double>(elapsed) >
	    ongoingStableTargetNanoseconds * 1.50L) {
	const float nextThreshold =
	    controller_point_proxy_threshold(
		this->d->lodPresentationPointProxyPixelThreshold,
		elapsed, 0, this->d->lodStableTargetFps);
	if (nextThreshold >
	    this->d->lodPresentationPointProxyPixelThreshold + 0.01f) {
	    this->d->lodPresentationPointProxyPixelThreshold =
		nextThreshold;
	    BObolViewLodState *presentationState =
		this->d->viewAttachment->getViewLodState();
	    if (presentationState)
		presentationState->
		    setCadPresentationPointProxyPixelThreshold(
			nextThreshold);
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	    /*
	     * This pressure step is deliberately conservative because its frame
	     * includes streaming/publication latency.  Mark the successor as a
	     * relaxation probe: once an unchanged prepared frame demonstrates
	     * headroom, walk back toward the finest sustainable presentation
	     * instead of permanently retaining the transient high-load cut.
	     */
	    this->d->lodStablePointProxyRelaxationActive = TRUE;
	    this->requestRender("lod-stream-point-pressure");
	    this->notifyFrameRequest("lod-stream-point-pressure");
	}
    }

    /* A budget-limited pass may have admitted a bounded first wave and then
     * scanned the remaining boxes without admitting them.  Re-plan from the
     * highest screen-value occurrence after that wave has actually rendered
     * and supplied a throughput sample.  Replanning before the frame would
     * repeatedly spend the same stale allowance. */
    if (this->d->lodBudgetRescanAfterFrame) {
	/*
	 * Calibration probes measure an unchanged retained cut.  They do not
	 * need an intervening O(N) occurrence-planning pass: that pass cannot
	 * admit anything until the samples have changed the aggregate budget.
	 * Present the bounded probe series back-to-back, then scan the sparse
	 * unsatisfied frontier once using the resulting calibration.
	 */
	if (this->d->lodStableCalibrationFramesRemaining > 0) {
	    this->d->lodStableCalibrationFramesRemaining--;
	    this->d->lodStableCalibrationProbeCount++;
	    this->requestRender("lod-budget-calibration");
	    this->notifyFrameRequest("lod-budget-calibration");
	} else {
	    this->d->lodBudgetRescanAfterFrame = FALSE;
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionPending = TRUE;
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodRefinementFaceBudgetInitialized = FALSE;
	    this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	    this->d->lodRetainedAdmissionActive = FALSE;
	    this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	    this->d->lodRefinementBudgetRetryMicroseconds = 0;
	    if (this->d->lodCoveragePassActive) {
		this->d->lodCoveragePassSawBoundedSource = FALSE;
		this->d->lodCoveragePassVisibleCount = 0;
		this->d->lodCoveragePassCoveredCount = 0;
	    }
	    this->markProgressiveWorkPending();
	}
    }
    this->d->renderCompletionSerial++;
    if (this->d->renderCompletionSerial == 0)
	this->d->renderCompletionSerial++;

    /* Rendering time is not user-idle time.  In particular an OSMesa motion
     * frame can take longer than the quiet debounce by itself.  Start the
     * quiet clock only after the coarse frame for the newest camera has
     * actually completed; otherwise wheel/trackpad input repeatedly jumps
     * back to an expensive stable cut between delivered events. */
    if (this->d->lodInteractive && !this->d->lodGestureActive &&
	this->d->lodSettleAfterRenderSerial != 0 &&
	this->d->renderCompletionSerial >=
	    this->d->lodSettleAfterRenderSerial) {
	this->d->lodLastViewChangeMicroseconds = bu_gettime();
	this->d->lodSettleAfterRenderSerial = 0;
    }

    /* A partial PoP result must be presented before the next, potentially
     * much larger prefix is requested.  This makes population staging an
     * actual user-visible progression and gives the frame-timing feedback a
     * chance to observe every step. */
    if (this->d->lodRefinementAwaitingFrame &&
	this->d->renderCompletionSerial >=
	    this->d->lodRefinementResumeAfterRenderSerial) {
	this->d->lodRefinementAwaitingFrame = FALSE;
	this->d->lodRefinementFaceBudgetInitialized = FALSE;
	this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	this->d->lodRetainedAdmissionActive = FALSE;
	this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	this->d->lodRefinementBudgetRetryMicroseconds = 0;
	/* Refinement is background quality work, not an invitation to monopolize
	 * the GUI thread.  A cheap hardware frame may proceed immediately.  Once
	 * a software or high-load frame exceeds two 60 Hz budgets, leave the
	 * useful cut visible for one measured render interval before exposing the
	 * next prefix (bounded at two seconds).
	 *
	 * The former four-render interval reduced a many-leaf warm stream from
	 * roughly 450 to 60-120 assets/second after publication frames reached
	 * 200 ms.  Interaction already switches to its coarse aggregate budget;
	 * imposing an additional 80% idle duty cycle on unchanged-view
	 * refinement only prolonged boxes.  One-render rest preserves a maximum
	 * 50% background presentation duty cycle and still adapts to OSMesa,
	 * System GL, viewport size, and current scene complexity without
	 * backend-specific policy. */
	const uint64_t responsiveFrame = 33333334ULL;
	int64_t cooldownMicroseconds = 0;
	if (elapsed > responsiveFrame) {
	    const uint64_t elapsedMicroseconds = elapsed / 1000ULL;
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
	this->markProgressiveWorkPending();
    }
    this->d->lodResultPublicationAwaitingFrame = FALSE;
    /* synchronizePresentation() ran before this completed frame, so it
     * contains every immutable result accumulated by the owner thread—not
     * merely the drain quantum which originally requested the frame. */
    this->d->lodUnpresentedResultCount = 0;
    this->d->lodUnpresentedResultFirstMicroseconds = 0;
    this->d->reconcilePhase();
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
    this->d->reconcilePhase();
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
	this->d->lodRefinementFaceBudgetInitialized = FALSE;
	this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	this->d->lodRetainedAdmissionActive = FALSE;
	this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	this->d->lodRefinementBudgetRetryMicroseconds = 0;
    }
    this->d->lastLodResultProcessingTimeNanoseconds = 0;
    this->d->lastProgressiveProviderTimeNanoseconds = 0;
    this->d->lastLodSubmissionTimeNanoseconds = 0;

    /* During camera motion favor responsiveness (4 px projected error).  Once
     * the view has been quiet for 150 ms, issue a new 1 px demand and let the
     * same per-view cache state refine in place. */
    if (this->d->lodAutoSubmit && this->d->lodInteractive) {
	const int64_t now = bu_gettime();
	const int64_t elapsed = now - this->d->lodLastViewChangeMicroseconds;
	if (!this->d->lodGestureActive &&
	    this->d->lodSettleAfterRenderSerial == 0 &&
	    this->d->lodLastViewChangeMicroseconds > 0 &&
	    elapsed >= 150000) {
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    const bool orthographicPresentation = this->d->activeCamera &&
		this->d->activeCamera->isOfType(
		    SoOrthographicCamera::getClassTypeId());
	    const bool unchangedOccurrencePopulation = viewState ?
		viewState->activeFaceCount() ==
		    this->d->lodStablePresentationActiveFaces &&
		viewState->cadRevision() ==
		    this->d->lodStablePresentationCadRevision :
		this->d->lodStablePresentationActiveFaces == 0;
	    if (orthographicPresentation && unchangedOccurrencePopulation &&
		!this->d->lodInteractionScaleChanged &&
		this->d->lodStablePresentationSnapshotValid &&
		!this->d->lodStablePresentationSnapshotRestored) {
		this->d->lodInteractiveProgressiveCeiling =
		    this->d->lodStablePresentationProgressiveCeiling;
		this->d->lodPresentationPointProxyPixelThreshold =
		    this->d->lodStablePresentationPointProxyPixelThreshold;
		this->d->lodStablePresentationSnapshotRestored = TRUE;
	    }
	    this->d->lodInteractive = FALSE;
	    this->d->lodViewScaleChanging = FALSE;
	    this->d->lodReleaseCutFloorActive = FALSE;
	    const bool haveRetainedMeshPayloads =
		this->getActiveLodMeshPayloadCount() > 0;
	    /* A motion ceiling also applies to native progressive wire stored in
	     * the standing CAD assembly.  It has no view-state mesh payload and
	     * therefore no occurrence pass capable of completing the mesh
	     * handoff below.  Restore its stable range directly; the policy
	     * revision already requests the one frame needed to present it. */
	    if (!haveRetainedMeshPayloads) {
		this->d->lodInteractiveProgressiveCeiling = -1;
		this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
		if (viewState) {
		    viewState->setCadPresentationProgressiveLodCeiling(-1);
		    viewState->setCadPresentationPointProxyPixelThreshold(1.0f);
		    viewState->setCadPresentationCameraMotionFrameReuse(FALSE);
		}
	    }
	    this->d->lodStablePresentationHandoffActive =
		haveRetainedMeshPayloads &&
		!this->d->lodStablePresentationSnapshotRestored &&
		(this->d->lodInteractiveProgressiveCeiling >= 0 ||
		 this->d->lodPresentationPointProxyPixelThreshold > 1.01f) ?
		    TRUE : FALSE;
	    this->d->lodStableCalibrationProbeActiveFaces = 0;
	    this->d->lodStableCalibrationProbeCount = 0;
	    this->d->lodStableCalibrationFramesRemaining = 0;
	    this->d->lodStableOverloadRecoveryPerformed = FALSE;
	    this->d->lodStableOverloadRecoveryActiveFaces = 0;
	    this->d->lodTargetPixelError = 1.0f;
	    this->d->lodStablePointProxyCalibrationPending = FALSE;
	    this->d->lodStablePointProxyRelaxationActive = FALSE;
	    if (viewState) {
		viewState->setCadPresentationCameraMotionFrameReuse(FALSE);
		if (this->d->lodStablePresentationSnapshotRestored) {
		    viewState->setCadPresentationProgressiveLodCeiling(
			this->d->lodInteractiveProgressiveCeiling);
		    viewState->setCadPresentationPointProxyPixelThreshold(
			this->d->lodPresentationPointProxyPixelThreshold);
		}
		/*
		 * activeFaceCount() is the retained occurrence population, not
		 * necessarily the cut just presented under the motion ceiling.
		 * Promoting the stable budget to that value authorized the hidden
		 * rich population before the first quiet capacity check.  Keep the
		 * measured presentation and reconcile the retained cuts first.
		 */
	    }
	    this->d->lodStablePresentationSnapshotValid = FALSE;
	    this->d->lodStablePresentationSnapshotRestored = FALSE;
	    this->d->lodInteractionScaleChanged = FALSE;
	    this->d->lodStablePresentationActiveFaces = 0;
	    this->d->lodStablePresentationCadRevision = 0;
	    this->advanceLodPolicyRevision();
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
	this->d->lodService->queuedResultCountForDiagnostics() : 0;
    const SbBool havePendingLodResults =
	this->hasPendingLodResults() || queuedLodResults > 0;
    SbBool coalesceLodResults = FALSE;
    const SbBool holdRicherResultsDuringInteraction =
	havePendingLodResults && this->d->lodInteractive &&
	this->d->lodVisibleCoverageComplete &&
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
	    this->d->lodService->inFlightCount() > 0 ||
	    this->d->lodService->pendingTaskCountForDiagnostics() > 0;
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
	 * there preserves the renderer's last prepared generation during camera
	 * input.  Once the quiet debounce ends they publish normally.  Coverage
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
	    this->d->lodService->queuedResultCountForDiagnostics() == 0) {
	    this->d->lodResultsFirstReadyMicroseconds.store(0);
	    /* A result can race the empty observation without producing another
	     * empty->non-empty callback.  Recheck after clearing and install a
	     * fresh age origin if necessary. */
	    if (this->d->lodService &&
		this->d->lodService->queuedResultCountForDiagnostics() > 0) {
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
	if (this->d->lodBudgetRescanAfterFrame)
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
	    this->d->lodService->pendingTaskCountForDiagnostics();
	localStatus.pendingTasks +=
	    this->d->lodService->delayedTaskCountForDiagnostics();
	localStatus.inFlight += this->d->lodService->inFlightCount();
	localStatus.queuedResults +=
	    this->d->lodService->queuedResultCountForDiagnostics();
	localStatus.queuedCacheWrites +=
	    this->d->lodService->queuedCacheWriteCountForDiagnostics();
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
	(this->d->lodService->pendingTaskCountForDiagnostics() > 0 ||
	 this->d->lodService->inFlightCount() > 0 ||
	 this->d->lodService->queuedResultCountForDiagnostics() > 0 ||
	 this->d->lodService->queuedCacheWriteCountForDiagnostics() > 0) ?
	1 : 0;

    /* Refinement and reclamation are separate phases.  A quiet view first
     * reaches its 1 px display target; only after a longer quiet interval and
     * an idle service do we replace this consumer's complete demand snapshot
     * and trim shared CPU prefixes to the aggregate maximum. */
    if (this->d->lodAutoSubmit &&
	this->d->lodResidentCompactionPending &&
	!this->d->lodInteractive &&
	this->d->lodResidentCompactionAfterMicroseconds > 0) {
	if (bu_gettime() <
		this->d->lodResidentCompactionAfterMicroseconds ||
	    pending_lod_realization_work ||
	    this->d->lodSubmissionPending) {
	    localStatus.hasMore = 1;
	} else if (this->d->lodService) {
	    std::vector<BObolLodResidentDemand> demands;
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    const uint64_t demandRevision = viewState ?
		viewState->residentMeshDemandRevision() :
		this->d->lodViewRevision.value();
	    const SbBool continuePlanning =
		this->d->lodResidentCompactionPlanning &&
		this->d->lodResidentCompactionDemandRevision ==
		    demandRevision;
	    if (viewState && !continuePlanning)
		viewState->residentMeshDemands(demands);
	    SbBool planningComplete = FALSE;
	    const size_t queued =
		this->d->lodService->scheduleResidentMeshCompaction(
		    residentConsumerId,
		    demandRevision, demands,
		    &planningComplete);
	    this->d->lodResidentCompactionPending =
		planningComplete ? FALSE : TRUE;
	    this->d->lodResidentCompactionPlanning =
		planningComplete ? FALSE : TRUE;
	    this->d->lodResidentCompactionDemandRevision =
		planningComplete ? 0 : demandRevision;
	    this->d->lodResidentCompactionAfterMicroseconds =
		planningComplete ? 0 : bu_gettime();
	    if (queued || !planningComplete)
		localStatus.hasMore = 1;
	}
    }

    if (this->d->lodUnpresentedResultCount > 0) {
	localStatus.hasMore = 1;
	const int64_t now = bu_gettime();
	if (!this->d->lodResultPublicationAwaitingFrame &&
	    this->d->lodUnpresentedResultsDue(now)) {
	    this->requestRender("lod-result-batch");
	    this->d->lodResultPublicationAwaitingFrame = TRUE;
	    if (this->d->lodRefinementAwaitingFrame)
		this->scheduleLodRefinementFrame("lod-result-batch");
	}
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
	if (!this->isRenderRequested())
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
	 this->d->lodResultPublicationAwaitingFrame))
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
    this->d->reconcilePhase(phasePending);

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
    BObolLodService *UNUSED(service), void *userData)
{
    BObolViewController *controller =
	static_cast<BObolViewController *>(userData);
    if (controller) {
	controller->d->lodResultsPending.store(1);
	int64_t expected = 0;
	(void)controller->d->lodResultsFirstReadyMicroseconds.
	    compare_exchange_strong(expected, bu_gettime());
	controller->markProgressiveWorkPending();
    }
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
    this->d->lodResultSubscriberId = 0;
    this->d->lodResultsPending.store(0);
    this->d->lodResultsFirstReadyMicroseconds.store(0);
    this->d->lodActiveGeneration = 0;
    this->d->lodSubmissionSourceIndex = 0;
    this->d->lodSubmissionEntryOffset = 0;
    this->d->clearLodSubmissionPlan();
    this->d->lodSubmissionPending = FALSE;
    this->d->lodSubmissionRescanPending = FALSE;
    this->d->lodVisibleCoverageComplete = FALSE;
    this->d->lodRetainedRefinementPending = FALSE;
    this->d->lodRetainedRefinementCutAdvanced = FALSE;
    this->d->lodRetainedRefinementBudgetBlocked = FALSE;
    this->d->lodBudgetRescanAfterFrame = FALSE;
    this->d->lodStableCalibrationFramesRemaining = 0;
    this->d->lodPassAdmittedWork = FALSE;
    this->d->lodRefinementFaceBudgetInitialized = FALSE;
    this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
    this->d->lodRetainedAdmissionActive = FALSE;
    this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
    this->d->lodRefinementBudgetRetryMicroseconds = 0;
    this->d->lodSettleAfterRenderSerial = 0;
    this->d->lodLastSubmittedViewRevision.reset();
    this->d->lodLastSubmittedPolicyRevision.reset();
    this->d->lodLastSubmittedSources.clear();
    this->d->lodSubmissionDeltaActive = FALSE;
    this->d->lodSubmissionDeltaSources.clear();
    this->d->lodSubmissionDeltaPlans.clear();
    this->d->lodResidentCompactionPending = service ? TRUE : FALSE;
    this->d->lodResidentCompactionPlanning = FALSE;
    this->d->lodResidentCompactionDemandRevision = 0;
    this->d->lodResidentAdmissionRevision =
	service ? service->residentMeshAdmissionRevision() : 0;
    this->d->lodResidentCompactionAfterMicroseconds =
	service ? bu_gettime() + 750000 : 0;
    this->d->lodRefinementAwaitingFrame = FALSE;
    this->d->lodRefinementResumeAfterRenderSerial = 0;
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodResultPublicationAwaitingFrame = FALSE;
    this->d->lodUnpresentedResultCount = 0;
    this->d->lodUnpresentedResultFirstMicroseconds = 0;

    if (this->d->lodService)
	this->d->lodResultSubscriberId =
	    this->d->lodService->subscribeResultReady(
		BObolViewController::lodResultReadyCB, this);
    this->d->reconcilePhase();
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
    this->d->lodRetainedRefinementPending = FALSE;
    this->d->lodRetainedRefinementCutAdvanced = FALSE;
    this->d->lodRetainedRefinementBudgetBlocked = FALSE;
    this->d->lodBudgetRescanAfterFrame = FALSE;
    this->d->lodStableCalibrationFramesRemaining = 0;
    this->d->lodPassAdmittedWork = FALSE;
    this->d->lodRefinementFaceBudgetInitialized = FALSE;
    this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
    this->d->lodRetainedAdmissionActive = FALSE;
    this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
    this->d->lodRefinementBudgetRetryMicroseconds = 0;
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
    this->d->lodResultPublicationAwaitingFrame = FALSE;
    this->d->lodUnpresentedResultCount = 0;
    this->d->lodUnpresentedResultFirstMicroseconds = 0;
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodInteractiveProgressiveCeiling = -1;
    this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodStablePointProxyRelaxationActive = FALSE;
    this->d->lodStablePresentationHandoffActive = FALSE;
    this->d->lodStablePresentationSnapshotValid = FALSE;
    this->d->lodStablePresentationSnapshotRestored = FALSE;
    this->d->lodInteractionScaleChanged = FALSE;
    this->d->lodStablePresentationProgressiveCeiling = -1;
    this->d->lodStablePresentationPointProxyPixelThreshold = 1.0f;
    this->d->lodStablePresentationActiveFaces = 0;
    this->d->lodStablePresentationCadRevision = 0;
    if (this->d->viewAttachment &&
	this->d->viewAttachment->getViewLodState())
	{
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    viewState->setCadPresentationProgressiveLodCeiling(-1);
	    viewState->setCadPresentationPointProxyPixelThreshold(1.0f);
	    viewState->setCadPresentationCameraMotionFrameReuse(FALSE);
	}
    this->d->reconcilePhase();
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
    if (!this->d->managedLodService)
	this->d->managedLodService.reset(new BObolLodService);
    if (this->d->managedLodService->isRunning() &&
	this->d->managedLodWorkerCount != workerCount)
	this->d->managedLodService->stop();
    if (!this->d->managedLodService->isRunning() &&
	!this->d->managedLodService->start(workerCount, TRUE))
	return NULL;
    this->d->managedLodWorkerCount = workerCount;
    this->setLodService(this->d->managedLodService.get());
    return this->d->managedLodService.get();
}

void
BObolViewController::stopManagedLodService(void)
{
    this->setLodAutoSubmit(FALSE);
    if (this->d->lodService == this->d->managedLodService.get())
	this->setLodService(NULL);
    if (this->d->managedLodService)
	this->d->managedLodService->stop();
    this->d->managedLodWorkerCount = 0;
}

size_t
BObolViewController::getManagedLodWorkerCount(void) const
{
    return this->d->managedLodWorkerCount;
}

void
BObolViewController::setLodAutoSubmit(SbBool enabled)
{
    this->d->lodAutoSubmit = enabled ? TRUE : FALSE;
    if (this->d->lodAutoSubmit) {
	this->requestRender("lod-auto-submit");
    } else {
	this->d->lodInteractiveProgressiveCeiling = -1;
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	this->d->lodStablePointProxyRelaxationActive = FALSE;
	this->d->lodStablePresentationHandoffActive = FALSE;
	this->d->lodStablePresentationSnapshotValid = FALSE;
	this->d->lodStablePresentationSnapshotRestored = FALSE;
	this->d->lodInteractionScaleChanged = FALSE;
	if (this->d->viewAttachment->getViewLodState()) {
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    viewState->setCadPresentationProgressiveLodCeiling(-1);
	    viewState->setCadPresentationPointProxyPixelThreshold(1.0f);
	    viewState->setCadPresentationCameraMotionFrameReuse(FALSE);
	}
    }
    this->d->reconcilePhase();
}

SbBool
BObolViewController::isLodAutoSubmitEnabled(void) const
{
    return this->d->lodAutoSubmit;
}

void
BObolViewController::setLodForcedLevel(int level)
{
    if (level < 0)
	level = 0;

    if (this->d->lodUseForcedLevel && this->d->lodForcedLevel == level)
	return;

    this->d->lodUseForcedLevel = TRUE;
    this->d->lodForcedLevel = level;
    this->advanceLodPolicyRevision();
    this->requestRender("lod-policy");
}

void
BObolViewController::clearLodForcedLevel(void)
{
    if (!this->d->lodUseForcedLevel)
	return;

    this->d->lodUseForcedLevel = FALSE;
    this->d->lodForcedLevel = 0;
    this->advanceLodPolicyRevision();
    this->requestRender("lod-policy");
}

SbBool
BObolViewController::hasLodForcedLevel(void) const
{
    return this->d->lodUseForcedLevel;
}

int
BObolViewController::getLodForcedLevel(void) const
{
    return this->d->lodForcedLevel;
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
	this->d->lodBudgetRescanAfterFrame ||
	this->d->lodStablePresentationHandoffActive ||
	this->d->lodStablePointProxyCalibrationPending;
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
	    (this->d->lodService &&
	     this->d->lodService->queuedResultCountForDiagnostics() > 0))
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
	if (this->d->lodSubmissionPending || this->d->lodResultsPending.load() != 0 ||
	    (this->d->lodService &&
	     this->d->lodService->queuedResultCountForDiagnostics() > 0))
	    this->markProgressiveWorkPending();
    };

    if (!this->hasPendingLodResults() &&
	this->d->lodService->queuedResultCountForDiagnostics() == 0) {
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
	this->d->lodCoveragePassActive = FALSE;
	this->d->lodCoveragePassSawBoundedSource = FALSE;
	this->d->lodCoveragePassVisibleCount = 0;
	this->d->lodCoveragePassCoveredCount = 0;
	this->d->lodLastSubmittedSources.clear();
	return 0;
    }

    if (this->d->lodLastSubmittedViewRevision == this->d->lodViewRevision &&
	this->d->lodLastSubmittedPolicyRevision == this->d->lodPolicyRevision &&
	signatures.sameInventories(this->d->lodLastSubmittedSources)) {
	if (this->d->lodSubmissionPending && this->d->lodActiveGeneration != 0)
	    return this->submitLodRequests(this->d->lodService,
		this->d->lodActiveGeneration, this->d->lodSubmissionRefreshMissing,
		this->d->lodSubmissionReset);
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
	this->d->lodCoveragePassActive = TRUE;
	this->d->lodVisibleCoverageComplete = FALSE;
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	this->d->lodStablePointProxyRelaxationActive = FALSE;
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
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodBudgetRescanAfterFrame = FALSE;
	this->d->lodStableCalibrationFramesRemaining = 0;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodRefinementFaceBudgetInitialized = FALSE;
	this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	this->d->lodRetainedAdmissionActive = FALSE;
	this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	this->d->lodRefinementBudgetRetryMicroseconds = 0;
	if (useSourceDelta) {
	    this->d->lodCoveragePassSawBoundedSource = FALSE;
	    this->d->lodCoveragePassVisibleCount = 0;
	    this->d->lodCoveragePassCoveredCount = 0;
	}
	if (sourceSetChanged || inventoryChanged)
	    this->d->lodStableOverloadRecoveryPerformed = FALSE;
	if (sourceSetChanged || inventoryChanged)
	    this->d->lodStableOverloadRecoveryActiveFaces = 0;
    } else if (!extendedPendingDelta) {
	this->d->lodSubmissionRescanPending = TRUE;
    }
    if (pendingDeltaNeedsFullRescan)
	this->d->lodSubmissionRescanPending = TRUE;
    this->d->lodSubmissionPending = TRUE;
    this->d->lodStableBudgetLimited = FALSE;
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
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodRefinementFaceBudgetInitialized = FALSE;
	this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	this->d->lodRetainedAdmissionActive = FALSE;
	this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	this->d->lodRefinementBudgetRetryMicroseconds = 0;
	this->d->lodCoveragePassSawBoundedSource = FALSE;
	this->d->lodCoveragePassVisibleCount = 0;
	this->d->lodCoveragePassCoveredCount = 0;
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
    if (!this->d->lodRefinementFaceBudgetInitialized) {
	const BObolViewLodState *lodState =
	    this->d->viewAttachment->getViewLodState();
	const size_t activeFaces = lodState ? lodState->activeFaceCount() : 0;
	const float targetFps = this->d->lodInteractive ?
	    this->d->lodInteractiveTargetFps :
	    this->d->lodStableTargetFps;
	const long double calibratedFacesPerSecond =
	    this->d->lodInteractive ?
	    this->d->lodInteractiveCalibratedFacesPerSecond :
	    this->d->lodStableCalibratedFacesPerSecond;
	const long double targetNanoseconds =
	    targetFps > 0.0f ?
	    1000000000.0L / static_cast<long double>(targetFps) : 0.0L;
	const long double calibratedFaceBudget =
	    targetFps > 0.0f && calibratedFacesPerSecond > 0.0L ?
	    calibratedFacesPerSecond /
		static_cast<long double>(targetFps) : 0.0L;
	const long double observedStableNanoseconds = std::max(
	    this->d->lodLastRenderWasPreparedCadReplay ?
		static_cast<long double>(
		    this->d->lastRenderTimeNanoseconds) : 0.0L,
	    static_cast<long double>(
		this->d->lodLastCadGpuTimeNanoseconds));
	const SbBool severeStableOverload =
	    !this->d->lodInteractive &&
	    targetNanoseconds > 0.0L &&
	    observedStableNanoseconds > targetNanoseconds * 2.0L;
	/* Quiet refinement is normally monotonic, but a discrete PoP jump can
	 * overshoot a nonlinear renderer ceiling.  Ordinary misses require
	 * three unchanged probes.  A frame over twice the target is already a
	 * severe user-visible stall and backs off immediately; treating that as
	 * a possible one-frame setup outlier makes renewed input wait behind
	 * repeated known-bad frames.  Recovery changes only retained draw cuts,
	 * never the cache or source geometry. */
	const SbBool stableOverloadRecovery =
	    !this->d->lodInteractive &&
	    !this->d->forceTerminalLodRefinement &&
	    (!this->d->lodStableOverloadRecoveryPerformed ||
	     this->d->lodStableOverloadRecoveryActiveFaces != activeFaces) &&
	    activeFaces > 0 &&
	    (this->d->lodStableCalibrationProbeCount >= 3 ||
	     severeStableOverload) &&
	    targetNanoseconds > 0.0L &&
	    observedStableNanoseconds > targetNanoseconds * 1.20L &&
	    (severeStableOverload ||
	     (calibratedFaceBudget > 0.0L &&
	      calibratedFaceBudget <
		static_cast<long double>(activeFaces) * 0.95L));
	size_t faceBudget = this->d->lodSeedFaceBudget;
	if (this->d->forceTerminalLodRefinement) {
	    faceBudget = SIZE_MAX;
	} else if (targetFps > 0.0f &&
	    calibratedFacesPerSecond > 0.0L) {
	    /* A measured presentation interval may already include a display
	     * cadence cap.  Multiplying every new sample by a headroom factor
	     * would then geometrically collapse the budget even while holding
	     * the target FPS.  The EMA and bounded growth handle noise; use the
	     * calibrated target directly. */
	    const long double affordable =
		calibratedFacesPerSecond /
		static_cast<long double>(targetFps);
	    faceBudget = affordable >= static_cast<long double>(SIZE_MAX) ?
		SIZE_MAX : static_cast<size_t>(
		    std::max<long double>(
			static_cast<long double>(
			    this->d->lodSeedFaceBudget), affordable));
	    /* Learn in bounded steps.  A fast empty/very-coarse hardware frame
	     * must not jump straight into a many-million-face cliff. */
	    if (activeFaces > 0 && faceBudget > activeFaces) {
		/* Quiet refinement uses smaller population steps so a throughput
		 * estimate from a cheap cut cannot overshoot into a much slower
		 * nonlinear rasterization regime.  Conversely, making a 1-3 ms
		 * retained System GL frame grow by only 25% forces thousands of
		 * instances through needless visible staging.  Scale the growth
		 * step from measured headroom: 4x below one quarter of the stable
		 * frame time, 2x below one half, and 1.25x near the limit.  Motion
		 * has no reason to spend multiple frames approaching its coarse
		 * responsive cut. */
		size_t growthMultipleNumerator = 4;
		size_t growthMultipleDenominator = 1;
		if (!this->d->lodInteractive) {
		    const long double growthTargetNanoseconds =
			targetFps > 0.0f ?
			1000000000.0L /
			    static_cast<long double>(targetFps) : 0.0L;
		    /* The faces/second calibration is already smoothed.  Growth
		     * tiering must describe the population currently on screen:
		     * an EMA spanning the earlier blank and minimum-cut frames
		     * can keep a now-cheap cut at 1.25x growth indefinitely. */
		    const uint64_t observedNanoseconds =
			static_cast<uint64_t>(observedStableNanoseconds);
		    if (growthTargetNanoseconds > 0.0L &&
			static_cast<long double>(observedNanoseconds) >=
			    growthTargetNanoseconds * 0.50L) {
			growthMultipleNumerator = 5;
			growthMultipleDenominator = 4;
		    } else if (growthTargetNanoseconds > 0.0L &&
			static_cast<long double>(observedNanoseconds) >=
			    growthTargetNanoseconds * 0.25L) {
			growthMultipleNumerator = 2;
		    }
		}
		/* A budget probe may intentionally render the same population
		 * because the next PoP level is a discrete jump.  Grow from the
		 * already-probed allowance in a quiet view, not solely from the
		 * unchanged active face count; otherwise every probe recomputes
		 * the same too-small budget and one occurrence can remain coarse
		 * forever despite ample calibrated capacity. */
		const size_t growthBase =
		    !this->d->lodInteractive &&
		    this->d->lodCurrentFaceBudget != SIZE_MAX ?
		    std::max(activeFaces, this->d->lodCurrentFaceBudget) :
		    activeFaces;
		const size_t growthIncrement =
		    growthBase > SIZE_MAX / growthMultipleNumerator ?
		    SIZE_MAX :
		    growthBase * growthMultipleNumerator /
			growthMultipleDenominator;
		const size_t growthLimit =
		    std::max(this->d->lodSeedFaceBudget, growthIncrement);
		faceBudget = std::min(faceBudget, growthLimit);
	    }
	    /* Refinement within one quiet camera epoch is monotonic.  A later
	     * slow frame may lower the learned capacity used by the next
	     * interaction, but must not make geometry which was already shown
	     * disappear from an otherwise unchanged view.  Monotonicity is also
	     * the correct noise filter: suppressing a small upward budget change
	     * can contradict the capacity probe's measured headroom and create
	     * an endless no-progress repaint loop.  Camera interaction remains
	     * free to select a coarser responsive cut. */
	    if (!this->d->lodInteractive && !stableOverloadRecovery &&
		this->d->lodCurrentFaceBudget != SIZE_MAX)
		faceBudget = std::max(faceBudget,
		    this->d->lodCurrentFaceBudget);
	    /* The throughput EMA necessarily remembers cheaper earlier cuts.
	     * Near the target frame time it can therefore remain slightly above
	     * the capacity of the population now on screen.  Letting that stale
	     * optimism raise the allowance by a few faces starts another full
	     * many-leaf scan, fills the tiny increment, and repeats forever.
	     *
	     * Ten percent frame-time headroom is the stable convergence margin:
	     * below it, preserve the current coherent cut and report the view as
	     * performance-limited.  Larger headroom still grows normally, and a
	     * later camera epoch recalibrates from the retained cut. */
	    if (!this->d->lodInteractive &&
		this->d->lodCurrentFaceBudget != SIZE_MAX &&
		faceBudget > this->d->lodCurrentFaceBudget &&
		targetNanoseconds > 0.0L) {
		const uint64_t observedNanoseconds =
		    static_cast<uint64_t>(observedStableNanoseconds);
		if (observedNanoseconds > 0 &&
		    static_cast<long double>(observedNanoseconds) >=
			targetNanoseconds * 0.90L)
		    faceBudget = this->d->lodCurrentFaceBudget;
	    }
	    if (stableOverloadRecovery &&
		observedStableNanoseconds > 0.0L) {
		/* Use the measured current population directly for backoff.
		 * The smoothed throughput intentionally lags normal noise, but
		 * that same lag retained a 187 ms frame for dozens of 50k-entry
		 * scans.  Target 80 percent of the frame-time allowance so the
		 * recovered cut leaves room for interaction and bookkeeping. */
		const long double recovered =
		    static_cast<long double>(activeFaces) *
		    targetNanoseconds * 0.80L /
		    observedStableNanoseconds;
		const size_t recoveredBudget =
		    recovered >= static_cast<long double>(SIZE_MAX) ?
			SIZE_MAX : static_cast<size_t>(
			    std::max<long double>(
				static_cast<long double>(
				    this->d->lodSeedFaceBudget),
				recovered));
		faceBudget = std::min(faceBudget, recoveredBudget);
	    }
	}
	/* The first interaction on an already-settled scene has no independent
	 * motion calibration yet.  Do not replace a demonstrably responsive cut
	 * with the nominal 50k seed.  Bootstrap an allowance from the current
	 * retained frame, with 20 percent timing headroom and bounded 4x growth.
	 * Pose-only motion will retain the current cut; zoom may spend the extra
	 * allowance while retargeting from that cut.  A later slow interaction
	 * frame immediately lowers the dedicated motion calibration. */
	if (this->d->lodInteractive && activeFaces > 0 &&
	    targetNanoseconds > 0.0L) {
	    uint64_t observedNanoseconds = std::max(
		this->d->lastRenderTimeNanoseconds,
		this->d->smoothedRenderTimeNanoseconds);
	    if (this->d->lodGestureActive) {
		std::lock_guard<std::mutex> lock(
		    this->d->presentationTimingMutex);
		observedNanoseconds = std::max(observedNanoseconds,
		    this->d->
			smoothedInteractivePresentationIntervalNanoseconds);
	    }
	    if (observedNanoseconds > 0 &&
		static_cast<long double>(observedNanoseconds) <=
		    targetNanoseconds * 1.05L) {
		const long double affordable =
		    static_cast<long double>(activeFaces) *
		    targetNanoseconds * 0.80L /
		    static_cast<long double>(observedNanoseconds);
		const size_t growthLimit =
		    activeFaces > SIZE_MAX / 4 ? SIZE_MAX : activeFaces * 4;
		const size_t responsiveBudget =
		    affordable >= static_cast<long double>(growthLimit) ?
			growthLimit :
			static_cast<size_t>(std::max<long double>(
			    static_cast<long double>(activeFaces),
			    affordable));
		faceBudget = std::max(faceBudget, responsiveBudget);
	    }
	}
	/* Button release does not define a new quality epoch.  Preserve the
	 * complete held-button cut through the debounce even if recomputing the
	 * interactive calibration would return the smaller nominal seed. */
	if (this->d->lodReleaseCutFloorActive &&
	    this->d->lodCurrentFaceBudget != SIZE_MAX)
	    faceBudget = std::max(faceBudget,
		this->d->lodCurrentFaceBudget);
	/*
	 * Translate retained-byte pressure into a conservative occurrence-face
	 * allowance.  Bytes are unique per asset while faces are charged per
	 * occurrence, so this is intentionally an iterative admission bound,
	 * not a claim of exact byte/triangle proportionality.  The following
	 * stable demand compaction measures the actual result and may schedule
	 * another bounded pass.  Targeting 80 percent leaves room for incoming
	 * worker publications and renderer-side copies.
	 */
	this->d->lodCurrentFaceBudget = faceBudget;
	const size_t additionalFaces =
	    faceBudget == SIZE_MAX ? SIZE_MAX :
	    (faceBudget > activeFaces ? faceBudget - activeFaces : 0);
	this->d->lodRefinementFaceBudgetRemaining = additionalFaces;
	/* Quiet cameras remain visually monotonic except for the confirmed
	 * overload case above.  Recovery reuses already resident PoP arrays and
	 * only lowers their draw counts, avoiding box flashes and reload churn. */
	this->d->lodRetainedAdmissionActive =
	    (this->d->lodInteractive || stableOverloadRecovery ||
	     this->d->lodStablePresentationHandoffActive) &&
	    faceBudget != SIZE_MAX && activeFaces > faceBudget ?
	    TRUE : FALSE;
	if (stableOverloadRecovery) {
	    this->d->lodStableOverloadRecoveryPerformed = TRUE;
	    this->d->lodStableOverloadRecoveryActiveFaces = activeFaces;
	}
	this->d->lodRetainedAdmissionBudgetRemaining =
	    this->d->lodRetainedAdmissionActive ? faceBudget : SIZE_MAX;
	this->d->lodRefinementBudgetRetryMicroseconds = 0;
	this->d->lodRefinementFaceBudgetInitialized = TRUE;
	if (getenv("BOBOL_LOD_TRACE_BUDGET"))
	    bu_log("BObol LoD frame budget active_faces=%zu "
		   "last_render_ms=%.3f smooth_render_ms=%.3f target_fps=%.3f "
		   "calibrated_mfaces_s=%.3f total_faces=%zu "
		   "additional_faces=%zu%s\n",
		   activeFaces,
		   this->d->lastRenderTimeNanoseconds / 1000000.0,
		   this->d->smoothedRenderTimeNanoseconds / 1000000.0,
		   targetFps,
		   static_cast<double>(
		       calibratedFacesPerSecond / 1000000.0L),
		   faceBudget, additionalFaces,
		   faceBudget == SIZE_MAX ? " unbounded" : "");
    }

    float scenePixelError = this->d->lodTargetPixelError;
    const BObolViewLodState *sceneLodState =
	this->d->viewAttachment->getViewLodState();
    const size_t sceneActiveFaces = sceneLodState ?
	sceneLodState->activeFaceCount() : 0;
    if (!this->d->lodUseForcedLevel &&
	this->d->lodCurrentFaceBudget != SIZE_MAX &&
	this->d->lodCurrentFaceBudget > 0 &&
	sceneActiveFaces > this->d->lodCurrentFaceBudget) {
	const long double over =
	    static_cast<long double>(sceneActiveFaces) /
	    static_cast<long double>(this->d->lodCurrentFaceBudget);
	scenePixelError *= static_cast<float>(std::sqrt(over));
    }

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_source_roots(this);
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
	const SbBool poseOnlyFullyResident =
	    this->d->lodInteractive &&
	    !this->d->lodViewScaleChanging &&
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
	double aspect = controller_aspect_from_region(this->d->viewportRegion);
	if (aspect <= SMALL_FASTF)
	    aspect = 1.0;
	const SbViewVolume volume = this->d->activeCamera->getViewVolume(
	    static_cast<float>(aspect));
	action.setViewVolume(&volume, scenePixelError);
	action.setGeneration(generation);
	action.setRevisions(this->d->lodViewRevision.value(),
	    this->d->lodPolicyRevision.value());
	action.setRefreshMissing(refreshMissing);
	action.setReset(reset);
	action.setViewLodState(viewLodState);
	/* Coarsening a retained prefix only changes the draw count/snap level;
	 * it does not rebuild or reread the mesh.  Permit it during motion so a
	 * previously settled, expensive cut cannot pin interactive FPS. */
	const SbBool scaleDemandChanged =
	    !this->d->lodInteractive || this->d->lodViewScaleChanging;
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
	    this->d->lodRetainedAdmissionActive &&
	    (!this->d->lodInteractive || this->d->lodViewScaleChanging);
	action.setAllowLevelDowngrade(
	    scaleDemandChanged || applyRetainedAdmission);
	action.setAllowRetainedRefinement(
	    (this->d->forceTerminalLodRefinement || scaleDemandChanged) &&
	    !this->d->lodRefinementAwaitingFrame ? TRUE : FALSE);
	action.setAllowRepresentationRefinement(
	    !this->d->lodInteractive &&
	    !this->d->lodRefinementAwaitingFrame ? TRUE : FALSE);
	/* A calibrated face budget governs richness above the minimum drawable
	 * mesh floor.  Returning a visible resident PoP payload to its box is
	 * both a visual regression and a false economy: the renderer can batch
	 * tiny minimum prefixes into the aggregate point channel without
	 * discarding the useful retained asset. */
	action.setPreserveMeshCoverage(TRUE);
	action.setRefinementFaceBudget(
	    this->d->lodRefinementFaceBudgetRemaining);
	/* Every bounded window owns a disjoint part of the pinned occurrence
	 * plan and consumes from this carried scene-wide remainder.  The action
	 * deliberately skips its full priority recovery when a finite window is
	 * configured, avoiding an O(scene size) input stall. */
	if (applyRetainedAdmission)
	    action.setRetainedSceneFaceBudget(
		this->d->lodRetainedAdmissionBudgetRemaining);
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
	const bool usingUnsatisfiedFrontier =
	    boundedLargeCompact &&
	    !this->d->lodSubmissionDeltaActive &&
	    !this->d->lodCoveragePassActive &&
	    this->d->lodVisibleCoverageComplete &&
	    !this->d->lodRetainedAdmissionActive &&
	    viewLodState;
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
	if (usingUnsatisfiedFrontier &&
	    (!this->d->lodSubmissionPlanValid ||
	     this->d->lodSubmissionPlanSource != source)) {
	    std::vector<SbString> unsatisfiedKeys;
	    viewLodState->unsatisfiedCadOccurrenceKeys(
		source, unsatisfiedKeys);
	    controller_prioritize_lod_frontier(source, viewLodState,
		unsatisfiedKeys, this->d->lodSubmissionPlanEntries);
	    this->d->lodSubmissionPlanSource = source;
	    this->d->lodSubmissionPlanValid = TRUE;
	}
	/* The retained frontier is ordered by the marginal transition which made
	 * it eligible.  Do not let one visit consume every remaining level and
	 * invalidate that ordering.  Coverage/index and exact delta scans keep
	 * their existing bulk behavior. */
	if (usingUnsatisfiedFrontier)
	    action.setTransitionLimitedRefinement(TRUE);
	if (boundedLargeCompact && this->d->lodCoveragePassActive) {
	    /*
	     * A bounded index-order window cannot know whether a later window
	     * still contains an uncovered visible leaf.  Do not let already
	     * covered early entries consume upward-refinement allowance until a
	     * complete projected pass proves coverage.  Missing entries still
	     * bind resident/cache data or submit their minimum provider request.
	     */
	    action.setAllowRetainedRefinement(FALSE);
	    this->d->lodCoveragePassSawBoundedSource = TRUE;
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
	if (this->d->lodUseForcedLevel)
	    action.setForcedLevel(this->d->lodForcedLevel);
	action.apply(source);
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
	if (this->d->lodCoveragePassActive) {
	    const size_t visible = action.getVisibleMeshCount();
	    const size_t covered = action.getCoveredVisibleMeshCount();
	    this->d->lodCoveragePassVisibleCount =
		visible > SIZE_MAX -
			this->d->lodCoveragePassVisibleCount ?
		    SIZE_MAX :
		    this->d->lodCoveragePassVisibleCount + visible;
	    this->d->lodCoveragePassCoveredCount =
		covered > SIZE_MAX -
			this->d->lodCoveragePassCoveredCount ?
		    SIZE_MAX :
		    this->d->lodCoveragePassCoveredCount + covered;
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
		       this->d->lodCoveragePassActive ? 1 : 0,
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
	/* Retained-admission changes are recovery/coarsening (including
	 * retiring an occurrence to its structural proxy), not upward PoP
	 * refinement.  Marking every such bounded slice as a richer cut routed
	 * it through the slow-frame refinement cooldown; a 50k scene could then
	 * spend nearly two seconds between each of 25 recovery windows. */
	if (action.getUpdatedCutCount() > 0 &&
	    !this->d->lodRetainedAdmissionActive)
	    this->d->lodRetainedRefinementCutAdvanced = TRUE;
	if (action.getPendingRetainedRefinementCount() > 0)
	    this->d->lodRetainedRefinementPending = TRUE;
	if (action.getRefinementBudgetBlockedCount() > 0)
	    this->d->lodRetainedRefinementBudgetBlocked = TRUE;
	if (getenv("BOBOL_LOD_TRACE_BUDGET") &&
	    (action.getRefinementFaceBudgetUsed() > 0 ||
	     action.getRefinementBudgetBlockedCount() > 0))
	    bu_log("BObol LoD frame budget source=%s used_faces=%zu "
		   "blocked=%u cuts=%u tasks=%u remaining_before=%zu\n",
		   source->path.getValue().getString(),
		   action.getRefinementFaceBudgetUsed(),
		   action.getRefinementBudgetBlockedCount(),
		   action.getUpdatedCutCount(), action.getSubmittedTaskCount(),
		   this->d->lodRefinementFaceBudgetRemaining);
	if (this->d->lodRefinementFaceBudgetRemaining != SIZE_MAX) {
	    const size_t used = action.getRefinementFaceBudgetUsed();
	    this->d->lodRefinementFaceBudgetRemaining =
		used >= this->d->lodRefinementFaceBudgetRemaining ?
		0 : this->d->lodRefinementFaceBudgetRemaining - used;
	}
	if (applyRetainedAdmission) {
	    const size_t used =
		action.getRetainedSceneFaceBudgetUsed();
	    this->d->lodRetainedAdmissionBudgetRemaining =
		used >= this->d->lodRetainedAdmissionBudgetRemaining ?
		0 : this->d->lodRetainedAdmissionBudgetRemaining - used;
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
    const bool completedPassChangedCut =
	completedPass &&
	this->d->lodRetainedRefinementCutAdvanced;
    if (completedPass && getenv("BOBOL_LOD_TRACE_PASS")) {
	static std::atomic<unsigned int> completedPassTraceCount(0);
	const unsigned int traceIndex =
	    completedPassTraceCount.fetch_add(1);
	if (traceIndex < 512)
	    bu_log("BObol LoD completed pass coverage=%d "
		   "coverage_seen=%d coverage_complete=%d delta=%d "
		   "handoff=%d refinement_pending=%d cut_advanced=%d "
		   "budget_rescan=%d active_faces=%zu face_budget=%zu "
		   "visited=%u cuts=%u\n",
		   this->d->lodCoveragePassActive ? 1 : 0,
		   this->d->lodCoveragePassSawBoundedSource ? 1 : 0,
		   this->d->lodVisibleCoverageComplete ? 1 : 0,
		   this->d->lodSubmissionDeltaActive ? 1 : 0,
		   this->d->lodStablePresentationHandoffActive ? 1 : 0,
		   this->d->lodRetainedRefinementPending ? 1 : 0,
		   this->d->lodRetainedRefinementCutAdvanced ? 1 : 0,
		   this->d->lodBudgetRescanAfterFrame ? 1 : 0,
		   sceneActiveFaces, this->d->lodCurrentFaceBudget,
		   this->d->lastLodVisitedMeshCount,
		   this->d->lastLodUpdatedCutCount);
    }
    if (completedPass && this->d->lodCoveragePassActive &&
	!this->d->lodCoveragePassSawBoundedSource &&
	!this->d->lodSubmissionRescanPending) {
	/* A targeted delta may address a small compact source.  Its single
	 * synchronous action is already a complete coverage pass; do not leave
	 * the global coverage state armed merely because no 2k window was
	 * needed. */
	this->d->lodVisibleCoverageComplete =
	    this->d->lodCoveragePassCoveredCount >=
		this->d->lodCoveragePassVisibleCount ? TRUE : FALSE;
	this->d->setLodConvergenceCompleteVisibleCount(
	    this->d->lodCoveragePassVisibleCount);
	if (!this->d->lodVisibleCoverageComplete) {
	    BObolViewLodState *coverageState =
		this->d->viewAttachment->getViewLodState();
	    if (coverageState)
		coverageState->setCadPresentationCameraMotionFrameReuse(
		    FALSE);
	}
	this->d->lodCoveragePassActive = FALSE;
	this->d->lodCoveragePassVisibleCount = 0;
	this->d->lodCoveragePassCoveredCount = 0;
    }
    const bool completedCoveragePass =
	completedPass && this->d->lodCoveragePassActive &&
	this->d->lodCoveragePassSawBoundedSource &&
	!this->d->lodSubmissionRescanPending;
    if (completedCoveragePass) {
	const bool missingCoverage =
	    this->d->lodCoveragePassCoveredCount <
		this->d->lodCoveragePassVisibleCount;
	if (getenv("BOBOL_LOD_TRACE_PASS"))
	    bu_log("BObol LoD completed coverage pass visible=%zu "
		   "covered=%zu missing=%d admitted=%d "
		   "budget_blocked=%d rescan_pending=%d\n",
		   this->d->lodCoveragePassVisibleCount,
		   this->d->lodCoveragePassCoveredCount,
		   missingCoverage ? 1 : 0,
		   this->d->lodPassAdmittedWork ? 1 : 0,
		   this->d->lodRetainedRefinementBudgetBlocked ? 1 : 0,
		   this->d->lodSubmissionRescanPending ? 1 : 0);
	this->d->setLodConvergenceCompleteVisibleCount(
	    this->d->lodCoveragePassVisibleCount);
	this->d->lodVisibleCoverageComplete =
	    missingCoverage ? FALSE : TRUE;
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodStableBudgetLimited = FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRefinementFaceBudgetInitialized = FALSE;
	this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	this->d->lodRetainedAdmissionActive = FALSE;
	this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	this->d->lodRefinementBudgetRetryMicroseconds = 0;
	this->d->lodCoveragePassSawBoundedSource = FALSE;
	this->d->lodCoveragePassVisibleCount = 0;
	this->d->lodCoveragePassCoveredCount = 0;
	if (missingCoverage) {
	    BObolViewLodState *coverageState =
		this->d->viewAttachment->getViewLodState();
	    if (coverageState)
		coverageState->setCadPresentationCameraMotionFrameReuse(
		    FALSE);
	    /*
	     * Capacity for the next minimum-prefix wave must be learned from a
	     * presented frame, just like ordinary scene-budget admission.
	     * completeRenderTiming() restarts this coverage pass with the new
	     * allowance.  Until then, pending provider work may continue on the
	     * worker pool without turning the GUI loop into a busy rescan.
	     */
	    this->d->lodSubmissionPending = FALSE;
	    this->d->lodBudgetRescanAfterFrame = TRUE;
	    this->requestRender("lod-coverage-admission");
	    this->notifyFrameRequest("lod-coverage-admission");
	} else {
	    /* Every projected leaf was already represented by a mesh.  Begin a
	     * fresh bounded pass which may now spend the remaining scene budget
	     * on screen-value-ordered PoP refinement. */
	    this->d->lodCoveragePassActive = FALSE;
	    this->d->lodSubmissionPending = TRUE;
	    this->markProgressiveWorkPending();
	}
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
	if (this->d->lodCoveragePassActive) {
	    this->d->lodCoveragePassSawBoundedSource = FALSE;
	    this->d->lodCoveragePassVisibleCount = 0;
	    this->d->lodCoveragePassCoveredCount = 0;
	}
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = sources.empty() ? FALSE : TRUE;
	this->d->lodStableBudgetLimited = FALSE;
    } else if (completedPass && this->d->lodRetainedRefinementPending &&
	this->d->lodRetainedRefinementCutAdvanced) {
	/* The just-completed pass selected the next resident cut using the
	 * newest view already.  Present it before any requested rescan; the
	 * post-frame submission is itself a full current-view pass. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = FALSE;
	this->d->lodStableBudgetLimited = FALSE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodStableCalibrationProbeActiveFaces = 0;
	this->d->lodStableCalibrationProbeCount = 0;
	this->d->lodStableCalibrationFramesRemaining = 0;
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
	    this->d->lodStableTargetFps > 0.0f ?
	    static_cast<uint64_t>(1000000000.0L /
		static_cast<long double>(this->d->lodStableTargetFps)) : 0;
	const uint64_t observedStableNanoseconds = std::max(
	    this->d->lodLastRenderWasPreparedCadReplay ?
		(this->d->lastRenderTimeNanoseconds ?
		    this->d->lastRenderTimeNanoseconds :
		    this->d->smoothedRenderTimeNanoseconds) : 0,
	    this->d->lodLastCadGpuTimeNanoseconds);
	long double calibratedStableFaceBudget = 0.0L;
	if (this->d->lodStableTargetFps > 0.0f &&
	    this->d->lodStableCalibratedFacesPerSecond > 0.0L)
	    calibratedStableFaceBudget =
		this->d->lodStableCalibratedFacesPerSecond /
		static_cast<long double>(this->d->lodStableTargetFps);
	const SbBool calibratedBudgetHeadroom =
	    calibratedStableFaceBudget >
		static_cast<long double>(this->d->lodCurrentFaceBudget) *
		    1.05L &&
	    observedStableNanoseconds <
		static_cast<uint64_t>(
		    static_cast<long double>(stableTargetNanoseconds) * 1.20L);
	const SbBool calibrationProbeCandidate =
	    !this->d->lodInteractive &&
	    !this->d->lodStablePresentationHandoffActive &&
	    !this->d->lodPassAdmittedWork &&
	    sceneActiveFaces > 0 &&
	    this->d->lodCurrentFaceBudget != SIZE_MAX &&
	    stableTargetNanoseconds > 0 &&
	    (observedStableNanoseconds <
		static_cast<uint64_t>(
		    static_cast<long double>(stableTargetNanoseconds) * 0.90L) ||
	     calibratedBudgetHeadroom);
	/* A probe is useful only while it can establish a larger affordable
	 * cut.  Bound probes at an unchanged visible population so floating
	 * calibration convergence, a discrete next-level jump, or contradictory
	 * thresholds cannot keep an otherwise stable retained view repainting
	 * forever.  Twelve 1.25x probes can span more than a 14x discrete jump
	 * while costing at most about 0.6 seconds at the default 20 FPS stable
	 * target; faster frames use the existing 2x/4x growth tiers. */
	if (this->d->lodPassAdmittedWork) {
	    this->d->lodStableCalibrationProbeActiveFaces = 0;
	    this->d->lodStableCalibrationProbeCount = 0;
	    this->d->lodStableCalibrationFramesRemaining = 0;
	} else if (this->d->lodStableCalibrationProbeActiveFaces !=
		   sceneActiveFaces) {
	    this->d->lodStableCalibrationProbeActiveFaces = sceneActiveFaces;
	    this->d->lodStableCalibrationProbeCount = 0;
	    this->d->lodStableCalibrationFramesRemaining = 0;
	}
	static constexpr unsigned int minimumStableCalibrationProbes = 3;
	static constexpr unsigned int maxStableCalibrationProbes = 12;
	const SbBool stableCalibrationEligible =
	    !this->d->lodInteractive &&
	    !this->d->lodStablePresentationHandoffActive &&
	    !this->d->lodPassAdmittedWork &&
	    sceneActiveFaces > 0 &&
	    this->d->lodCurrentFaceBudget != SIZE_MAX &&
	    stableTargetNanoseconds > 0;
	/* Never classify a newly reached quiet cut from a single frame.  OSMesa
	 * in particular can pay shader/cache/setup work in the first traversal,
	 * while a compositor stall can perturb System GL.  Three samples remain
	 * a small bounded cost at the stable target and let subsequent ordinary
	 * traversal cost govern the retained budget. */
	const SbBool calibrationProbe =
	    stableCalibrationEligible &&
	    (this->d->lodStableCalibrationProbeCount <
		minimumStableCalibrationProbes ||
	     calibrationProbeCandidate) &&
	    this->d->lodStableCalibrationProbeCount <
		maxStableCalibrationProbes;
	if (getenv("BOBOL_LOD_TRACE_PASS")) {
	    static std::atomic<unsigned int> passTraceCount(0);
	    const unsigned int traceIndex = passTraceCount.fetch_add(1);
	    if (traceIndex < 256)
		bu_log("BObol LoD completed budget pass "
		       "active_faces=%zu face_budget=%zu admitted=%d "
		       "probe_candidate=%d probe_eligible=%d "
		       "probe_count=%u probe_next=%d "
		       "observed_ms=%.3f target_ms=%.3f "
		       "calibrated_mfaces_s=%.3f\n",
		       sceneActiveFaces, this->d->lodCurrentFaceBudget,
		       this->d->lodPassAdmittedWork ? 1 : 0,
		       calibrationProbeCandidate ? 1 : 0,
		       stableCalibrationEligible ? 1 : 0,
		       this->d->lodStableCalibrationProbeCount,
		       calibrationProbe ? 1 : 0,
		       observedStableNanoseconds / 1000000.0,
		       stableTargetNanoseconds / 1000000.0,
		       static_cast<double>(
			   this->d->lodStableCalibratedFacesPerSecond /
			       1000000.0L));
	}
	if (calibrationProbe) {
	    this->d->lodStableCalibrationProbeCount++;
	    const unsigned int targetProbeCount =
		calibrationProbeCandidate ?
		    maxStableCalibrationProbes :
		    minimumStableCalibrationProbes;
	    this->d->lodStableCalibrationFramesRemaining =
		targetProbeCount > this->d->lodStableCalibrationProbeCount ?
		    targetProbeCount -
			this->d->lodStableCalibrationProbeCount : 0;
	} else if (!calibrationProbeCandidate) {
	    this->d->lodStableCalibrationProbeCount = 0;
	    this->d->lodStableCalibrationFramesRemaining = 0;
	}
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = FALSE;
	this->d->lodBudgetRescanAfterFrame =
	    this->d->lodPassAdmittedWork || calibrationProbe;
	this->d->lodStableBudgetLimited =
	    this->d->lodBudgetRescanAfterFrame ? FALSE : TRUE;
	if (this->d->lodBudgetRescanAfterFrame) {
	    const char *frameReason = calibrationProbe ?
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
	if (!this->d->lodBudgetRescanAfterFrame &&
	    this->d->lodStableBudgetLimited &&
	    !this->d->lodInteractive &&
	    sceneActiveFaces > this->d->lodCurrentFaceBudget &&
	    stableTargetNanoseconds > 0 &&
	    observedStableNanoseconds >
		static_cast<long double>(stableTargetNanoseconds) * 1.20L) {
	    const float nextThreshold =
		controller_point_proxy_threshold(
		    this->d->lodPresentationPointProxyPixelThreshold,
		    static_cast<uint64_t>(observedStableNanoseconds), 0,
		    this->d->lodStableTargetFps);
	    if (nextThreshold >
		this->d->lodPresentationPointProxyPixelThreshold + 0.01f) {
		this->d->lodPresentationPointProxyPixelThreshold =
		    nextThreshold;
		BObolViewLodState *presentationState =
		    this->d->viewAttachment->getViewLodState();
		if (presentationState)
		    presentationState->
			setCadPresentationPointProxyPixelThreshold(
			    nextThreshold);
		this->d->lodStablePointProxyCalibrationPending = TRUE;
		/*
		 * The multiplicative pressure correction intentionally lands on
		 * the safe side of the target.  Its next unchanged frame is a
		 * bounded relaxation probe so terminal quality is the finest cut
		 * which meets the stable FPS contract, not the first coarse cut
		 * which happens to do so.
		 */
		this->d->lodStablePointProxyRelaxationActive = TRUE;
		this->requestRender("lod-stable-point-calibration");
		this->notifyFrameRequest("lod-stable-point-calibration");
	    }
	}
	append_controller_lod_diagnostic(this->d->lastLodDiagnostics,
	    SbString(""),
	    this->d->lodStablePointProxyCalibrationPending ?
		"scene LoD calibrating small-part aggregation" :
	    this->d->lodBudgetRescanAfterFrame ?
		(calibrationProbe ?
		    "scene LoD probing stable calibrated capacity" :
		    "scene LoD admission awaiting calibrated frame") :
		"scene LoD reached its calibrated face budget");
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRefinementFaceBudgetInitialized = FALSE;
	this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	this->d->lodRetainedAdmissionActive = FALSE;
	this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	this->d->lodRefinementBudgetRetryMicroseconds = 0;
    } else {
	this->d->lodSubmissionPending = completedPass ? FALSE : TRUE;
	if (completedPass)
	    this->d->lodPassAdmittedWork = FALSE;
	if (completedPass)
	    this->d->lodStableBudgetLimited = FALSE;
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
    if (completedPass &&
	this->d->lodStablePresentationHandoffActive &&
	!this->d->lodSubmissionPending &&
	!this->d->lodBudgetRescanAfterFrame &&
	!completedPassChangedCut) {
	this->d->lodStablePresentationHandoffActive = FALSE;
	this->d->lodInteractiveProgressiveCeiling = -1;
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	if (presentationState) {
	    presentationState->setCadPresentationProgressiveLodCeiling(-1);
	    presentationState->setCadPresentationCameraMotionFrameReuse(FALSE);
	}
	this->d->lodStablePointProxyCalibrationPending = TRUE;
	this->d->lodStablePointProxyRelaxationActive =
	    this->d->lodPresentationPointProxyPixelThreshold > 1.01f ?
		TRUE : FALSE;
	this->requestRender("lod-stable-handoff");
	this->notifyFrameRequest("lod-stable-handoff");
    }
    if (!this->d->lodSubmissionPending &&
	this->d->lodRetainedRefinementCutAdvanced &&
	(this->d->lodRetainedRefinementPending ||
	 this->d->lodStablePresentationHandoffActive)) {
	/* A richer prefix is already in memory, but expose it one cut per
	 * completed frame.  This presentation barrier is also required when
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
	this->d->lodRefinementFaceBudgetInitialized = FALSE;
	this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	this->d->lodRetainedAdmissionActive = FALSE;
	this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	this->d->lodRefinementBudgetRetryMicroseconds = 0;
	this->scheduleLodRefinementFrame("lod-cut");
    }
    if (completedPass && !this->d->lodSubmissionPending &&
	!this->d->lodRetainedRefinementPending) {
	this->d->lodRefinementFaceBudgetInitialized = FALSE;
	this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	this->d->lodRetainedAdmissionActive = FALSE;
	this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	this->d->lodRefinementBudgetRetryMicroseconds = 0;
    }
    if (this->d->lodSubmissionPending)
	this->markProgressiveWorkPending();

    if (this->d->lastLodUpdatedCutCount > 0) {
	this->d->lodResidentCompactionPending = TRUE;
	this->d->lodResidentCompactionAfterMicroseconds =
	    bu_gettime() + 750000;
	this->requestRender("lod-cut");
    }
    if (completedPass && this->d->lodSubmissionDeltaActive &&
	!this->d->lodSubmissionPending &&
	!this->d->lodBudgetRescanAfterFrame &&
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
				       size_t maxEstimatedBytes)
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

    std::vector<BObolLodResult> drained;
    this->d->lastLodResultCount = service->drainResults(
	drained, maxResults, maxEstimatedBytes);
    this->d->lodResultsPending.store(
	service->queuedResultCountForDiagnostics() > 0 ? 1 : 0);
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
	const SbBool traceResult = controller_lod_trace_result(drained[i]);
	if (traceResult) {
	    const BObolViewLodState::CadPayload *resident =
		findResidentCad();
	    bu_log("BObol LoD apply trace object=%s request_level=%d "
		   "loaded_level=%d incoming_view=%llu incoming_policy=%llu "
		   "current_view=%llu current_policy=%llu resident_level=%d "
		   "resident_view=%llu resident_policy=%llu\n",
		   drained[i].request.objectName.getString(),
		   drained[i].request.requestedLevel,
		   drained[i].geometry.activeLevel,
		   static_cast<unsigned long long>(
		       drained[i].request.viewRevision.value()),
		   static_cast<unsigned long long>(
		       drained[i].request.policyRevision.value()),
		   static_cast<unsigned long long>(
		       this->d->lodViewRevision.value()),
		   static_cast<unsigned long long>(
		       this->d->lodPolicyRevision.value()),
		   resident ? resident->activeLevel : -1,
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
			       "occurrence=%s incoming_level=%d requested=%d "
			       "incoming_view=%llu incoming_policy=%llu "
			       "current_view=%llu current_policy=%llu "
			       "resident_level=%d resident_view=%llu "
			       "resident_policy=%llu supersedes=%d\n",
			       drained[i].request.objectName.getString(),
			       drained[i].request.occurrenceKey.getString(),
			       drained[i].geometry.activeLevel,
			       drained[i].request.requestedLevel,
			       static_cast<unsigned long long>(
				   drained[i].request.viewRevision.value()),
			       static_cast<unsigned long long>(
				   drained[i].request.policyRevision.value()),
			       static_cast<unsigned long long>(
				   this->d->lodViewRevision.value()),
			       static_cast<unsigned long long>(
				   this->d->lodPolicyRevision.value()),
			       cadPayload ? cadPayload->activeLevel :
				   (meshPayload ?
				       meshPayload->activeLevel : -1),
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
	const SbBool partialRefinementResult =
	    drained[i].providerStatus == BOBOL_LOD_PROVIDER_READY &&
	    drained[i].resultKind == BOBOL_LOD_RESULT_MESH &&
	    drained[i].geometry.activeLevel >= 0 &&
	    drained[i].request.requestedLevel >
		drained[i].geometry.activeLevel &&
	    !drained[i].terminal &&
	    (drained[i].request.sourceRoutingId == 0 ||
	     !residentCadBefore ||
	     residentCadBefore->activeLevel <
		 drained[i].geometry.activeLevel ||
	     (residentMeshBefore &&
	      residentMeshBefore->activeLevel <
		  drained[i].geometry.activeLevel)) ?
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
				   "object=%s occurrence=%s level=%d "
				   "requested=%d progressive=%d valid=%d "
				   "route=%p diagnostic=%s\n",
				   drained[i].request.objectName.getString(),
				   drained[i].request.occurrenceKey.getString(),
				   drained[i].geometry.activeLevel,
				   drained[i].request.requestedLevel,
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
	service->queuedResultCountForDiagnostics() > 0 ? 1 : 0);

    if (this->d->lastLodAppliedResultCount > 0) {
	const size_t applied =
	    static_cast<size_t>(this->d->lastLodAppliedResultCount);
	if (applied > SIZE_MAX - this->d->lodUnpresentedResultCount)
	    this->d->lodUnpresentedResultCount = SIZE_MAX;
	else
	    this->d->lodUnpresentedResultCount += applied;
	const int64_t now = bu_gettime();
	if (this->d->lodUnpresentedResultFirstMicroseconds <= 0)
	    this->d->lodUnpresentedResultFirstMicroseconds = now;
	const bool firstUsefulMesh =
	    this->getActiveLodMeshPayloadCount() <= applied;
	const bool streamIdle =
	    service->queuedResultCountForDiagnostics() == 0 &&
	    service->pendingTaskCountForDiagnostics() == 0 &&
	    service->inFlightCount() == 0;
	const bool publishNow =
	    firstUsefulMesh || streamIdle ||
	    this->d->lodUnpresentedResultsDue(now);
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
	if (publishNow) {
	    this->requestRender("lod-result");
	    this->d->lodResultPublicationAwaitingFrame = TRUE;
	}
	this->d->lodResidentCompactionPending = TRUE;
	this->d->lodResidentCompactionAfterMicroseconds =
	    bu_gettime() + 750000;
	(void)this->enforceMeshResidencyBudget();
	if (partialRefinementCandidate && publishNow)
	    this->scheduleLodRefinementFrame("lod-result");
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
    this->d->lodPolicyRevision.set(newRevision);
    /*
     * A quality-policy revision changes the requested PoP cut, not camera
     * visibility.  Retain the last complete current-view denominator until a
     * camera or source-inventory revision invalidates it.  Clearing it here
     * let the quiet transition report a 50k scene converged with a zero
     * visible target while its unsatisfied-frontier fast path deliberately
     * skipped another full visibility scan.
     */
    this->d->resetLodConvergenceFraction();
    this->d->lodStableBudgetLimited = FALSE;
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodStablePointProxyRelaxationActive = FALSE;
    this->requestRender("lod-policy");
    this->d->reconcilePhase();
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
    this->d->clearLodConvergenceCandidates();
    this->d->resetLodConvergenceFraction();
    this->d->lodStableBudgetLimited = FALSE;
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodStablePointProxyRelaxationActive = FALSE;
    this->d->lodStablePresentationHandoffActive = FALSE;
    this->d->lodCoveragePassActive = TRUE;
    this->d->lodVisibleCoverageComplete = FALSE;
    this->d->reconcilePhase();
}

void
BObolViewController::advanceLodPolicyRevision(void)
{
    this->d->lodPolicyRevision.advance();
    /* Quality changes preserve the current view's proven visibility
     * denominator.  Source and camera revisions clear it explicitly. */
    this->d->resetLodConvergenceFraction();
    this->d->lodStableBudgetLimited = FALSE;
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodStablePointProxyRelaxationActive = FALSE;
    this->d->lodCoveragePassActive = TRUE;
    this->d->lodVisibleCoverageComplete = FALSE;
    this->d->reconcilePhase();
}

static float
controller_interactive_pixel_error(uint64_t renderNanoseconds,
				   uint64_t presentationNanoseconds,
				   float targetFps)
{
    /* Preserve the already-displayed cut when it meets the configured motion
     * target.  Camera motion is not itself a reason to discard detail.  Once
     * measured render/presentation cost misses the target, increase error
     * with the square root of the overrun: PoP population commonly grows near
     * quadratically.  The two-pixel first pressure tier removes one PoP bit;
     * severe overload may climb to 64 pixels. */
    const double targetNanoseconds =
	std::isfinite(targetFps) && targetFps > 0.0f ?
	1000000000.0 / static_cast<double>(targetFps) :
	1000000000.0 / 60.0;
    const uint64_t observed = std::max(renderNanoseconds,
	presentationNanoseconds);
    if (!observed)
	return 4.0f;
    if (static_cast<double>(observed) <= targetNanoseconds * 1.05)
	return 1.0f;
    const double scale = std::sqrt(
	static_cast<double>(observed) / targetNanoseconds);
    return static_cast<float>(std::max(2.0, std::min(64.0, 2.0 * scale)));
}

static float
controller_point_proxy_threshold(
    float currentThreshold, uint64_t renderNanoseconds,
    uint64_t presentationNanoseconds, float targetFps)
{
    const double targetNanoseconds =
	std::isfinite(targetFps) && targetFps > 0.0f ?
	1000000000.0 / static_cast<double>(targetFps) :
	1000000000.0 / 60.0;
    const uint64_t observed = std::max(
	renderNanoseconds, presentationNanoseconds);
    const double current = std::max(1.0,
	std::min(64.0, static_cast<double>(currentThreshold)));
    if (!observed ||
	static_cast<double>(observed) <= targetNanoseconds * 1.05)
	return static_cast<float>(current);

    /*
     * The screen threshold describes the cut which produced the measured
     * frame.  Recomputing an absolute threshold from that frame reaches a
     * false fixed point when per-object command work dominates: a 4 px cut
     * can remain at 40-60 ms forever.  Apply the measured correction to the
     * current cut, just as the discrete PoP ceiling ratchets downward.
     * Quiet recovery is the only path which lowers this threshold.
     */
    const double correction = std::sqrt(
	static_cast<double>(observed) / targetNanoseconds);
    const double absoluteSeed =
	controller_interactive_pixel_error(
	    renderNanoseconds, presentationNanoseconds, targetFps);
    return static_cast<float>(std::min(64.0,
	std::max(absoluteSeed, current * correction)));
}

static int
controller_interactive_progressive_ceiling(
    BObolViewLodState *viewState, float pixelError, int currentCeiling = -1)
{
    if (!viewState)
	return -1;
    const int activeMaximum =
	viewState->maximumActiveProgressiveLevel();
    if (activeMaximum < 0)
	return -1;
    const double error = std::max(1.0,
	static_cast<double>(pixelError));
    if (error <= 1.0 + std::numeric_limits<double>::epsilon())
	return -1;
    /*
     * PoP levels are discrete.  Rounding the required reduction down kept a
     * cut which was already known to miss the motion target whenever the
     * estimate fell just below the next power of two (3.99 pixels was the
     * common 50k-scene case).  Round toward the conservative coarser cut;
     * quiet recovery removes this render-only ceiling after the debounce and
     * refines from the retained arrays without a reload.
     */
    const int reduction = std::max(1,
	static_cast<int>(std::ceil(std::log2(error))));
    /*
     * A later camera event observes the cost of the cut already constrained
     * by currentCeiling, not the original level-15/active-maximum cut.  Apply
     * each measured correction relative to the cut which produced it.
     * Reapplying every correction to activeMaximum made level 5 a fixed point
     * in the 50k scene even though that retained cut still took 130-160 ms.
     */
    const int measuredMaximum = currentCeiling >= 0 ?
	std::min(activeMaximum, currentCeiling) : activeMaximum;
    return std::max(0, measuredMaximum - reduction);
}

void
BObolViewController::beginLodInteraction(void)
{
    if (!this->d->lodAutoSubmit || this->d->lodGestureActive)
	return;

    const float previousPixelError = this->d->lodTargetPixelError;
    if (!this->d->lodInteractive) {
	this->d->lodStablePresentationSnapshotValid = TRUE;
	this->d->lodStablePresentationSnapshotRestored = FALSE;
	this->d->lodInteractionScaleChanged = FALSE;
	this->d->lodStablePresentationProgressiveCeiling =
	    this->d->lodInteractiveProgressiveCeiling;
	this->d->lodStablePresentationPointProxyPixelThreshold =
	    this->d->lodPresentationPointProxyPixelThreshold;
	const BObolViewLodState *snapshotState =
	    this->d->viewAttachment->getViewLodState();
	this->d->lodStablePresentationActiveFaces = snapshotState ?
	    snapshotState->activeFaceCount() : 0;
	this->d->lodStablePresentationCadRevision = snapshotState ?
	    snapshotState->cadRevision() : 0;
    }
    this->d->lodGestureActive = TRUE;
    this->d->lodViewScaleChanging = FALSE;
    this->d->lodReleaseCutFloorActive = FALSE;
    this->d->lodStablePresentationHandoffActive = FALSE;
    this->d->lodStableCalibrationProbeActiveFaces = 0;
    this->d->lodStableCalibrationProbeCount = 0;
    this->d->lodStableCalibrationFramesRemaining = 0;
    this->d->lodStableOverloadRecoveryPerformed = FALSE;
    this->d->lodStableOverloadRecoveryActiveFaces = 0;
    this->d->lodStablePointProxyCalibrationPending = FALSE;
    this->d->lodStablePointProxyRelaxationActive = FALSE;
    this->d->lodInteractive = TRUE;
    this->d->lodSettleAfterRenderSerial = 0;
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodRefinementFaceBudgetInitialized = FALSE;
    this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
    this->d->lodRetainedAdmissionActive = FALSE;
    this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
    this->d->lodRefinementBudgetRetryMicroseconds = 0;
    this->d->lodLastViewChangeMicroseconds = bu_gettime();
    uint64_t presentationNanoseconds = 0;
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
	presentationNanoseconds = 0;
    }
    this->d->lodTargetPixelError =
	controller_interactive_pixel_error(
	    this->d->smoothedRenderTimeNanoseconds,
	    presentationNanoseconds,
	    this->d->lodInteractiveTargetFps);
    this->d->lodPresentationPointProxyPixelThreshold =
	controller_point_proxy_threshold(
	    this->d->lodPresentationPointProxyPixelThreshold,
	    this->d->lastRenderTimeNanoseconds,
	    presentationNanoseconds,
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
    const int responsiveCeiling =
	controller_interactive_progressive_ceiling(
	    viewState, this->d->lodTargetPixelError,
	    this->d->lodInteractiveProgressiveCeiling);
    if (responsiveCeiling >= 0 &&
	(this->d->lodInteractiveProgressiveCeiling < 0 ||
	 responsiveCeiling < this->d->lodInteractiveProgressiveCeiling))
	this->d->lodInteractiveProgressiveCeiling = responsiveCeiling;
    this->d->lodInteractiveCeilingFeedbackRenderSerial =
	this->d->renderCompletionSerial;
    if (viewState)
	viewState->setCadPresentationProgressiveLodCeiling(
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
    this->d->reconcilePhase();
}

void
BObolViewController::endLodInteraction(void)
{
    if (!this->d->lodGestureActive)
	return;

    this->d->lodGestureActive = FALSE;
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
    const bool orthographicPresentation = this->d->activeCamera &&
	this->d->activeCamera->isOfType(
	    SoOrthographicCamera::getClassTypeId());
    const bool unchangedOccurrencePopulation = viewState ?
	viewState->activeFaceCount() ==
	    this->d->lodStablePresentationActiveFaces &&
	viewState->cadRevision() ==
	    this->d->lodStablePresentationCadRevision :
	this->d->lodStablePresentationActiveFaces == 0;
    if (orthographicPresentation && unchangedOccurrencePopulation &&
	!this->d->lodInteractionScaleChanged &&
	this->d->lodStablePresentationSnapshotValid) {
	this->d->lodInteractiveProgressiveCeiling =
	    this->d->lodStablePresentationProgressiveCeiling;
	this->d->lodPresentationPointProxyPixelThreshold =
	    this->d->lodStablePresentationPointProxyPixelThreshold;
	this->d->lodStablePresentationSnapshotRestored = TRUE;
	this->d->lodReleaseCutFloorActive = FALSE;
	if (viewState) {
	    viewState->setCadPresentationProgressiveLodCeiling(
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
    if (!this->d->lodStablePresentationSnapshotRestored)
	this->d->lodReleaseCutFloorActive = TRUE;
    /* A gesture pass may have exhausted the old, smaller allowance.  Force
     * the release pass to derive its admission state from the coherent floor
     * above rather than reusing that stale remainder. */
    this->d->lodRefinementFaceBudgetInitialized = FALSE;
    this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
    this->d->lodRetainedAdmissionActive = FALSE;
    this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
    this->d->lodRefinementBudgetRetryMicroseconds = 0;
    /* Keep the interaction epoch through the normal quiet-view debounce.
     * Pose-only presentation may already be restored, but new projection and
     * refinement work still waits so release cannot block on a full software
     * planning frame. */
    this->d->lodLastViewChangeMicroseconds = bu_gettime();
    this->d->lodSettleAfterRenderSerial = 0;
    this->markProgressiveWorkPending();
    this->requestRender("lod-interaction-end");
    this->d->reconcilePhase();
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
    return this->d->lodInteractive && this->d->lodViewScaleChanging ?
	TRUE : FALSE;
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
    this->d->lodRefinementFaceBudgetInitialized = FALSE;
    this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
    this->d->lodRetainedAdmissionActive = FALSE;
    this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
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
BObolViewController::getCurrentLodFaceBudget(void) const
{
    return this->d->lodCurrentFaceBudget;
}

size_t
BObolViewController::getActiveLodFaceCount(void) const
{
    const BObolViewLodState *state = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    return state ? state->activeFaceCount() : 0;
}

void
BObolViewController::getLodConvergenceStatus(
	BObolLodConvergenceStatus &status) const
{
    status.clear();
    static_assert(
	static_cast<int>(BObolLodPhaseTracker::Phase::FALLBACK) ==
	    BOBOL_LOD_COORDINATOR_FALLBACK &&
	static_cast<int>(BObolLodPhaseTracker::Phase::COVERAGE) ==
	    BOBOL_LOD_COORDINATOR_COVERAGE &&
	static_cast<int>(BObolLodPhaseTracker::Phase::INTERACTIVE) ==
	    BOBOL_LOD_COORDINATOR_INTERACTIVE &&
	static_cast<int>(BObolLodPhaseTracker::Phase::SETTLING) ==
	    BOBOL_LOD_COORDINATOR_SETTLING &&
	static_cast<int>(BObolLodPhaseTracker::Phase::STABLE) ==
	    BOBOL_LOD_COORDINATOR_STABLE &&
	static_cast<int>(BObolLodPhaseTracker::Phase::COMPACTING) ==
	    BOBOL_LOD_COORDINATOR_COMPACTING,
	"public and private LoD coordinator phases must agree");
    const BObolLodPhaseTracker::Phase coordinatorPhase =
	this->d->currentPhase();
    status.coordinatorPhase = static_cast<int>(coordinatorPhase);
    status.coordinatorTransitionSerial =
	this->d->phaseTransitionSerial();
    status.coordinatorProgressSequence =
	this->d->phaseWitness(coordinatorPhase).sequence;
    status.viewRevision = this->d->lodViewRevision.value();
    status.failedSourceCount = this->getLastFailedSourceCount();
    status.visibleTargetCount =
	this->d->lodConvergenceCandidateCount();
    status.activeFaces = this->getActiveLodFaceCount();
    status.faceBudget = this->getCurrentLodFaceBudget();

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
    if (viewState)
	viewState->convergencePayloadCounts(status.activePayloadCount,
	    status.satisfiedPayloadCount,
	    status.memoryLimitedPayloadCount);

    if (this->d->lodService) {
	status.pendingTasks =
	    this->d->lodService->pendingTaskCountForDiagnostics();
	const size_t delayed =
	    this->d->lodService->delayedTaskCountForDiagnostics();
	status.pendingTasks = delayed > SIZE_MAX - status.pendingTasks ?
	    SIZE_MAX : status.pendingTasks + delayed;
	status.inFlight = this->d->lodService->inFlightCount();
	status.queuedResults =
	    this->d->lodService->queuedResultCountForDiagnostics();
	status.queuedCacheWrites =
	    this->d->lodService->queuedCacheWriteCountForDiagnostics();
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
    const SbBool calibrationPending =
	this->hasPendingLodRefinementFrame();
    status.refinementFramePending =
	this->d->lodRefinementAwaitingFrame;
    status.budgetCalibrationPending =
	this->d->lodBudgetRescanAfterFrame;
    status.stablePresentationHandoffPending =
	this->d->lodStablePresentationHandoffActive;
    status.pointProxyCalibrationPending =
	this->d->lodStablePointProxyCalibrationPending;
    const SbBool visualPending =
	structuralPending || this->d->lodSubmissionPending ||
	resultPending || status.pendingTasks > 0 || status.inFlight > 0 ||
	calibrationPending;

    status.viewReady =
	!visualPending && !this->d->lodInteractive ? TRUE : FALSE;
    status.backgroundPending =
	status.queuedCacheWrites > 0 ||
	this->d->lodResidentCompactionPending ||
	(status.viewReady && this->hasProgressiveWorkPending()) ?
	TRUE : FALSE;
    status.memoryLimited =
	status.memoryLimitedPayloadCount > 0 ||
	(status.residentMeshLimitBytes != SIZE_MAX &&
	 status.stableResidentMeshBytes >
	    status.residentMeshLimitBytes) ?
	TRUE : FALSE;
    status.performanceLimited =
	status.viewReady &&
	(this->d->lodStableBudgetLimited ||
	 status.memoryLimited ||
	 (status.visibleTargetCount > 0 &&
	  (status.activePayloadCount < status.visibleTargetCount ||
	   status.satisfiedPayloadCount < status.visibleTargetCount))) ?
	TRUE : FALSE;

    const bool hasLodState =
	status.expectedLeafCount > 0 || status.availableLeafCount > 0 ||
	status.visibleTargetCount > 0 || status.activePayloadCount > 0 ||
	visualPending || status.backgroundPending;
    if (!hasLodState) {
	status.phase = BOBOL_LOD_CONVERGENCE_IDLE;
	status.fraction = 1.0f;
	return;
    }

    if (structuralPending || structuralDiscovery) {
	status.phase = BOBOL_LOD_CONVERGENCE_DISCOVERING;
	const long double coverage =
	    structuralPending && status.expectedLeafCount > 0 ?
	    static_cast<long double>(status.availableLeafCount) /
		static_cast<long double>(status.expectedLeafCount) : 0.0L;
	status.fraction = static_cast<float>(
	    std::min<long double>(0.40L, 0.40L * coverage));
    } else if (this->d->lodInteractive) {
	status.phase = BOBOL_LOD_CONVERGENCE_INTERACTIVE;
	if (status.visibleTargetCount == 0) {
	    /* A new camera revision clears the projected target plan before
	     * the first bounded scan.  Do not score stale payloads against
	     * themselves and briefly claim 85% before the real denominator
	     * arrives. */
	    status.fraction = 0.40f;
	} else {
	    const size_t target = status.visibleTargetCount;
	    const long double quality =
		static_cast<long double>(
		    std::min(status.satisfiedPayloadCount, target)) /
		    static_cast<long double>(target);
	    status.fraction = static_cast<float>(
		0.40L + std::min<long double>(0.45L, 0.45L * quality));
	}
    } else if (calibrationPending) {
	status.phase = BOBOL_LOD_CONVERGENCE_CALIBRATING;
	status.fraction = 0.95f;
    } else if (visualPending) {
	status.phase = BOBOL_LOD_CONVERGENCE_REFINING;
	size_t target = std::max(
	    status.visibleTargetCount, status.activePayloadCount);
	if (!target) {
	    target = status.pendingTasks > SIZE_MAX - status.inFlight ?
		SIZE_MAX : status.pendingTasks + status.inFlight;
	    target = std::max<size_t>(1, target);
	}
	const long double quality =
	    static_cast<long double>(
		std::min(status.satisfiedPayloadCount, target)) /
	    static_cast<long double>(target);
	status.fraction = static_cast<float>(
	    0.40L + std::min<long double>(0.55L, 0.55L * quality));
    } else if (status.backgroundPending) {
	status.phase = BOBOL_LOD_CONVERGENCE_BACKGROUND;
	status.fraction = 1.0f;
    } else {
	status.phase = BOBOL_LOD_CONVERGENCE_IDLE;
	status.fraction = 1.0f;
    }

    if (status.failedSourceCount > 0)
	status.phase = BOBOL_LOD_CONVERGENCE_ERROR;
    status.fraction = std::max(0.0f, std::min(1.0f, status.fraction));

    /* Planning discovers the denominator in bounded waves and may raise a
     * payload's requested cut after first coverage.  Both can make a raw
     * completed/target ratio fall even though a fixed camera view has only
     * improved.  Keep the user-facing fraction monotonic within one
     * view/policy revision; a camera or policy change explicitly resets it. */
    if (this->d->lodConvergenceFractionViewRevision !=
	    this->d->lodViewRevision ||
	this->d->lodConvergenceFractionPolicyRevision !=
	    this->d->lodPolicyRevision) {
	this->d->lodConvergenceFractionViewRevision =
	    this->d->lodViewRevision;
	this->d->lodConvergenceFractionPolicyRevision =
	    this->d->lodPolicyRevision;
	this->d->lodConvergenceFractionFloor = 0.0f;
    }
    if (status.phase != BOBOL_LOD_CONVERGENCE_ERROR) {
	status.fraction = std::max(status.fraction,
	    this->d->lodConvergenceFractionFloor);
	this->d->lodConvergenceFractionFloor = status.fraction;
    }
}

double
BObolViewController::getCalibratedLodFacesPerSecond(void) const
{
    const long double value = this->d->lodInteractive ?
	this->d->lodInteractiveCalibratedFacesPerSecond :
	this->d->lodStableCalibratedFacesPerSecond;
    return value >= static_cast<long double>(
	std::numeric_limits<double>::max()) ?
	std::numeric_limits<double>::max() : static_cast<double>(value);
}

double
BObolViewController::getInteractiveCalibratedLodFacesPerSecond(void) const
{
    const long double value =
	this->d->lodInteractiveCalibratedFacesPerSecond;
    return value >= static_cast<long double>(
	std::numeric_limits<double>::max()) ?
	std::numeric_limits<double>::max() : static_cast<double>(value);
}

double
BObolViewController::getStableCalibratedLodFacesPerSecond(void) const
{
    const long double value = this->d->lodStableCalibratedFacesPerSecond;
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
	this->d->lodViewScaleChanging = scaleChanged;
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
	    if (!continuingInteractive) {
		this->d->lodStablePresentationSnapshotValid = TRUE;
		this->d->lodStablePresentationSnapshotRestored = FALSE;
		this->d->lodInteractionScaleChanged = FALSE;
		this->d->lodStablePresentationProgressiveCeiling =
		    this->d->lodInteractiveProgressiveCeiling;
		this->d->lodStablePresentationPointProxyPixelThreshold =
		    this->d->lodPresentationPointProxyPixelThreshold;
		const BObolViewLodState *snapshotState =
		    this->d->viewAttachment->getViewLodState();
		this->d->lodStablePresentationActiveFaces = snapshotState ?
		    snapshotState->activeFaceCount() : 0;
		this->d->lodStablePresentationCadRevision = snapshotState ?
		    snapshotState->cadRevision() : 0;
	    }
	    if (scaleChanged)
		this->d->lodInteractionScaleChanged = TRUE;
	    const int64_t now = bu_gettime();
	    this->d->lodLastViewChangeMicroseconds = now;
	    this->d->lodInteractive = TRUE;
	    this->d->lodReleaseCutFloorActive = FALSE;
	    this->d->lodStablePresentationHandoffActive = FALSE;
	    this->d->lodStableCalibrationProbeActiveFaces = 0;
	    this->d->lodStableCalibrationProbeCount = 0;
	    this->d->lodStableCalibrationFramesRemaining = 0;
	    this->d->lodStableOverloadRecoveryPerformed = FALSE;
	    this->d->lodStableOverloadRecoveryActiveFaces = 0;
	    this->d->lodRefinementNotBeforeMicroseconds = 0;
	    this->d->lodRefinementFaceBudgetInitialized = FALSE;
	    this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	    this->d->lodRetainedAdmissionActive = FALSE;
	    this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	    this->d->lodRefinementBudgetRetryMicroseconds = 0;
	    this->d->lodSettleAfterRenderSerial =
		this->d->renderCompletionSerial + 1;
	    if (this->d->lodSettleAfterRenderSerial == 0)
		this->d->lodSettleAfterRenderSerial = 1;
	    this->d->lodResidentCompactionPending = TRUE;
	    this->d->lodResidentCompactionAfterMicroseconds = now + 750000;
	    /* Preserve a pending cut's frame gate across a newer camera
	     * signature.  The action may still coarsen immediately, but it must
	     * not expose another finer prefix before any intervening frame. */
	    uint64_t presentationNanoseconds = 0;
	    if (this->d->lodGestureActive) {
		std::lock_guard<std::mutex> lock(
		    this->d->presentationTimingMutex);
		presentationNanoseconds = this->d->
		    smoothedInteractivePresentationIntervalNanoseconds;
	    }
	    /* A wheel/trackpad camera update is event-driven unless the host
	     * explicitly brackets it as a gesture.  The gap since an unrelated
	     * stable repaint measures user/test cadence, not renderer capacity;
	     * feeding it to the motion policy can turn an ordinary 4-8 px cut
	     * into a destructive 64 px cut.  Render time remains valid for
	     * one-shot changes, while actual drag gestures retain measured
	     * presentation cadence for queued System GL work. */
	    this->d->lodTargetPixelError =
		controller_interactive_pixel_error(
		    this->d->smoothedRenderTimeNanoseconds,
		    presentationNanoseconds,
		    this->d->lodInteractiveTargetFps);
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
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
	     * timings back into controller_point_proxy_threshold() used to grow
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
	    if (hasNewRenderFeedback && !continuingInteractive)
		this->d->lodPresentationPointProxyPixelThreshold =
		    controller_point_proxy_threshold(
			this->d->lodPresentationPointProxyPixelThreshold,
			this->d->lastRenderTimeNanoseconds,
			presentationNanoseconds,
			this->d->lodInteractiveTargetFps);
	    if (hasNewRenderFeedback) {
		responsiveCeiling =
		    controller_interactive_progressive_ceiling(
			viewState, this->d->lodTargetPixelError,
			this->d->lodInteractiveProgressiveCeiling);
		this->d->lodInteractiveCeilingFeedbackRenderSerial =
		    this->d->renderCompletionSerial;
	    }
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
		viewState->setCadPresentationProgressiveLodCeiling(
		    this->d->lodInteractiveProgressiveCeiling);
	    if (viewState)
		viewState->setCadPresentationPointProxyPixelThreshold(
		    this->d->lodPresentationPointProxyPixelThreshold);
	    this->markProgressiveWorkPending();
	}
    }
    this->d->reconcilePhase();
}
