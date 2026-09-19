/*            V I E W _ C O N T R O L L E R _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file view_controller_host.cpp
 *
 * View-controller lifetime, camera synchronization, and host render-request
 * plumbing.  Progressive LoD policy and execution live in their dedicated
 * controller units.
 */

#include "common.h"

#include "bu/log.h"
#include "bv.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BViewController.h"
#include "cad_assembly_private.h"
#include "database_source_realization.h"
#include "scalar_publication_private.h"
#include "view_controller_private.h"
#include "raytrace.h"
#include "rt/view.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <string.h>
#include <vector>

#include <Inventor/SoDB.h>
#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>
#include <Inventor/tools/SbModernUtils.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/gl.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoClipPlane.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSeparator.h>

static const double controller_orthographic_camera_distance_scale = 5.0;
static const double controller_camera_near_scale = 0.001;
static const double controller_camera_far_scale = 100.0;
static const double controller_camera_minimum_distance = 1.0;
static const double controller_camera_minimum_near_distance = 1.0e-6;
static constexpr float controller_default_camera_z = 10.0f;
static constexpr float controller_default_camera_height = 2.0f;
static constexpr float controller_default_camera_near = 1.0f;
static constexpr float controller_default_camera_far = 20.0f;
static constexpr float controller_default_camera_focal = 10.0f;

static SbModernUtils::SoNodeRef
controller_default_camera(void)
{
    SbModernUtils::SoNodeRef cameraOwner(new SoOrthographicCamera);
    auto *camera = static_cast<SoOrthographicCamera *>(cameraOwner.get());
    camera->position.setValue(0.0f, 0.0f, controller_default_camera_z);
    camera->height = controller_default_camera_height;
    camera->nearDistance = controller_default_camera_near;
    camera->farDistance = controller_default_camera_far;
    camera->focalDistance = controller_default_camera_focal;
    return cameraOwner;
}

struct ControllerCallbackDispatchState;
struct ControllerCallbackDispatchFrame {
    ControllerCallbackDispatchState *state;
    ControllerCallbackDispatchFrame *previous;
};
static thread_local ControllerCallbackDispatchFrame *
    controllerActiveCallbackDispatch = NULL;

struct ControllerCallbackDispatchState {
    ControllerCallbackDispatchState(void) :
	dispatches(0),
	closing(false),
	deleteAfterDispatch(false)
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

    bool finishDispatch(void)
    {
	std::lock_guard<std::mutex> lock(this->mutex);
	this->dispatches--;
	if (!this->dispatches)
	    this->cv.notify_all();
	return !this->dispatches && this->deleteAfterDispatch;
    }

    bool close(void)
    {
	std::unique_lock<std::mutex> lock(this->mutex);
	this->closing = true;
	ControllerCallbackDispatchFrame *frame =
	    controllerActiveCallbackDispatch;
	while (frame && frame->state != this)
	    frame = frame->previous;
	if (frame) {
	    /* A callback may unsubscribe itself.  Its dispatch reference remains
	     * live until the callback returns, so waiting here would deadlock and
	     * deleting here would make finishDispatch use freed storage. */
	    this->deleteAfterDispatch = true;
	    return false;
	}
	this->cv.wait(lock, [this]() {
	    return this->dispatches == 0;
	});
	return true;
    }

    std::mutex mutex;
    std::condition_variable cv;
    unsigned int dispatches;
    bool closing;
    bool deleteAfterDispatch;
};

struct ControllerFrameRequestState : ControllerCallbackDispatchState {
    ControllerFrameRequestState(BObolFrameRequestCallback callback_in,
	void *user_data_in) :
	callback(callback_in),
	userData(user_data_in)
    {
    }

    BObolFrameRequestCallback callback;
    void *userData;
};

struct ControllerPresentationSyncState : ControllerCallbackDispatchState {
    ControllerPresentationSyncState(
	BObolPresentationSyncCallback callback_in, void *user_data_in) :
	callback(callback_in),
	userData(user_data_in)
    {
    }

    BObolPresentationSyncCallback callback;
    void *userData;
};

template <typename State>
class ControllerCallbackDispatchScope {
public:
    explicit ControllerCallbackDispatchScope(State *state_in) :
	state(state_in),
	frame {state_in, controllerActiveCallbackDispatch}
    {
	controllerActiveCallbackDispatch = &this->frame;
    }

    ~ControllerCallbackDispatchScope()
    {
	controllerActiveCallbackDispatch = this->frame.previous;
	if (this->state->finishDispatch())
	    delete this->state;
    }

    ControllerCallbackDispatchScope(
	const ControllerCallbackDispatchScope &) = delete;
    ControllerCallbackDispatchScope &operator=(
	const ControllerCallbackDispatchScope &) = delete;

private:
    State *state;
    ControllerCallbackDispatchFrame frame;
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

bool
controller_lod_source_inputs_unsubmitted(
    const std::vector<SoBRLDatabaseSource *> &sources,
    const std::vector<BObolLodSourceSnapshot> &submitted)
{
    for (SoBRLDatabaseSource *source : sources) {
	if (!source || !source->hasDisplayLodTargets())
	    continue;
	const uint64_t routingId = source->getCompactSourceRoutingId();
	const auto found = std::find_if(submitted.begin(), submitted.end(),
	    [routingId](const BObolLodSourceSnapshot &snapshot) {
		return snapshot.routingId.value() == routingId;
	    });
	const uint64_t inventoryRevision =
	    source->getDisplayMeshLodRevision();
	const uint64_t visibilityRevision =
	    source->getDisplayMeshLodVisibilityRevision();
	if (found == submitted.end() ||
	    found->inventoryRevision.value() != inventoryRevision ||
	    found->visibilityRevision != visibilityRevision) {
	    if (getenv("BOBOL_LOD_TRACE_SOURCE_CONTRACT"))
		bu_log("BObol LoD source contract unsubmitted source=%p "
		       "path=%s inventory=%llu/%llu visibility=%llu/%llu\n",
		       static_cast<void *>(source),
		       source->path.getValue().getString(),
		       static_cast<unsigned long long>(
			   found == submitted.end() ? 0 :
			   found->inventoryRevision.value()),
		       static_cast<unsigned long long>(inventoryRevision),
		       static_cast<unsigned long long>(
			   found == submitted.end() ? 0 :
			   found->visibilityRevision),
		       static_cast<unsigned long long>(visibilityRevision));
	    return true;
	}
    }
    return false;
}
/*
 * All normal views in a process share resident mesh assets, worker threads,
 * cache writes, and memory governors.  Each controller still owns an isolated
 * generation and resident-demand consumer.  The weak broker lets the service
 * shut down naturally after the last view releases it.
 */
std::shared_ptr<BObolLodService>
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
    this->initializeControllerState(NULL, NULL, TRUE);
}

BObolViewController::BObolViewController(SoNode *root, SoCamera *camera) :
    d(new Impl(this))
{
    this->initializeControllerState(root, camera, FALSE);
}

void
BObolViewController::initializeControllerState(SoNode *root, SoCamera *camera,
	SbBool createDefaultRoot)
{
    try {
	this->d->viewAttachment->ref();
	SbModernUtils::SoNodeRef batchOwner(new SoBRLCadRenderBatch);
	auto *batch = static_cast<SoBRLCadRenderBatch *>(batchOwner.get());
	batch->addChild(this->d->viewport->getRoot());
	SbModernUtils::SoNodeRef presentationOwner(new SoSeparator);
	auto *presentation = static_cast<SoSeparator *>(presentationOwner.get());
	SbModernUtils::SoNodeRef underlayOwner(new SoGroup);
	auto *underlay = static_cast<SoGroup *>(underlayOwner.get());
	/* Unlike underlay/overlay this root is parented by GED's per-view render
	 * composition, so the controller keeps a reference across rebinds. */
	SbModernUtils::SoNodeRef interlayOwner(new SoGroup);
	SbModernUtils::SoNodeRef overlayOwner(new SoGroup);
	auto *overlay = static_cast<SoGroup *>(overlayOwner.get());
	presentation->addChild(underlay);
	presentation->addChild(batch);
	presentation->addChild(overlay);
	controller_initialize_render_action(this->d->renderManager);
	controller_configure_render_environment(this->d->viewport);
	this->d->renderBatchRoot = batchOwner.release();
	this->d->renderPresentationRoot = presentationOwner.release();
	this->d->framebufferUnderlayRoot = underlay;
	this->d->framebufferInterlayRoot =
	    static_cast<SoGroup *>(interlayOwner.release());
	this->d->framebufferOverlayRoot = overlay;

	SbModernUtils::SoNodeRef defaultRoot(createDefaultRoot ?
	    static_cast<SoNode *>(new SoBRLSceneGroup) : NULL);
	SbModernUtils::SoNodeRef defaultCamera(createDefaultRoot ?
	    controller_default_camera() : SbModernUtils::SoNodeRef(NULL));
	this->setSceneRoot(createDefaultRoot ? defaultRoot.get() : root);
	this->setCamera(createDefaultRoot ?
	    static_cast<SoCamera *>(defaultCamera.get()) : camera);
	if (createDefaultRoot) {
	    /* An empty construction root and private viewport camera make the
	     * default controller immediately renderable without putting camera
	     * state in modeled scene content. Construction creates no visible work
	     * for an unattached endpoint. */
	    this->clearRenderRequest();
	    this->clearProgressiveWorkPending();
	}
    } catch (...) {
	this->releaseControllerStateNoexcept();
	throw;
    }
}

BObolViewController::~BObolViewController(void)
{
    this->releaseControllerStateNoexcept();
}

void
BObolViewController::releaseControllerStateNoexcept(void) noexcept
{
    if (!this->d)
	return;
    try { this->setPresentationSyncCallback(NULL, NULL); }
    catch (...) {}
    ControllerFrameRequestState *frameRequestState = NULL;
    {
	std::lock_guard<std::mutex> lock(this->d->frameRequestMutex);
	frameRequestState = static_cast<ControllerFrameRequestState *>(
	    this->d->frameRequestUserData);
	this->d->frameRequestCallback = NULL;
	this->d->frameRequestUserData = NULL;
    }
    if (frameRequestState && frameRequestState->close())
	delete frameRequestState;
    delete this->d->imageRenderer;
    this->d->imageRenderer = NULL;
    this->d->imageRendererManager = NULL;
    /* Live provider removal republishes convergence state. During destruction
     * there is no successor state to publish, and that calculation allocates. */
    for (const BObolProgressiveProviderRecord &record :
	 this->d->progressiveProviders) {
	if (!record.userDataFree)
	    continue;
	try { (*record.userDataFree)(record.userData); }
	catch (...) {}
    }
    this->d->progressiveProviders.clear();
    /* Unsubscribe before retiring service work.  This is the worker-callback
     * quiescence barrier and must not be preceded by the fallible live
     * setLodService() transition: construction failure can destroy a
     * partially initialized controller while allocation denial is active.
     * Once no callback can address this controller, cancellation and resident
     * demand release are independent best-effort service cleanup. */
    BObolLodService *lodService = this->d->lodService;
    const BObolLodSubscriberId lodSubscriber =
	this->d->lodResultSubscriberId;
    const uint64_t lodGeneration = this->d->lodActiveGeneration;
    if (lodService && lodSubscriber)
	lodService->unsubscribeResultReady(lodSubscriber);
    this->d->lodResultSubscriberId = 0;
    this->d->lodActiveGeneration = 0;
    this->d->lodService = NULL;
    if (lodService && lodGeneration) {
	try { lodService->cancelGeneration(lodGeneration); }
	catch (...) {}
    }
    if (lodService) {
	try {
	    lodService->releaseResidentMeshConsumer(
		this->d->residentMeshConsumerId());
	} catch (...) {}
    }
    this->d->managedLodService.reset();
    this->d->managedLodWorkerCount = 0;
    try { this->clearRtPickCaches(); }
    catch (...) {}
    delete this->d->featureStore;
    this->d->featureStore = NULL;
    delete this->d->polygonStore;
    this->d->polygonStore = NULL;
    delete this->d->selectionStore;
    this->d->selectionStore = NULL;
    /* SoViewport releases its camera during its own teardown. Retire the
     * controller's separate reference directly; setCamera() repairs live
     * lighting and root order and therefore allocates. */
    if (this->d->activeCamera) {
	this->d->activeCamera->unref();
	this->d->activeCamera = NULL;
    }
    /* Destruction retires retained graph references directly.  setSceneRoot()
     * is a live publication operation: it prepares repository membership,
     * rebuilds render composition, synchronizes the manager, and requests a
     * frame.  None of those successor-state effects exists during teardown. */
    if (this->d->viewAttachment)
	this->d->viewAttachment->setSceneRoot(NULL);
    SoBRLCadRenderBatch *batch =
	dynamic_cast<SoBRLCadRenderBatch *>(this->d->renderBatchRoot);
    if (batch)
	batch->setBatchSourceRoot(NULL);
    if (this->d->renderLodRoot) {
	this->d->viewport->setSceneGraph(NULL);
	this->d->renderLodRoot->unref();
	this->d->renderLodRoot = NULL;
    }
    if (this->d->viewAttachment) {
	this->d->viewAttachment->unref();
	this->d->viewAttachment = NULL;
    }
    if (this->d->renderManager) {
	this->d->renderManager->setSceneGraph(NULL);
	delete this->d->renderManager;
	this->d->renderManager = NULL;
    }
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
BObolViewController::setSceneRoot(SoNode *root, SbBool preserveLodState)
{
    if (!preserveLodState)
	this->cancelActiveLodGeneration();
    this->clearRtPickCaches();
    this->d->orthographicDepthReferenceSize = 0.0;
    this->d->orthographicDepthReferenceStructuralRevision = 0;
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
    this->requestLodCapacityRender("scene-root");
}

SoNode *
BObolViewController::getSceneRoot(void) const
{
    return this->d->viewAttachment->getSceneRoot();
}

void
BObolViewController::setRenderSceneRoot(SoNode *root, SbBool preserveLodState)
{
    if (!preserveLodState)
	this->cancelActiveLodGeneration();
    this->clearRtPickCaches();
    if (!preserveLodState)
	this->d->viewAttachment->clearViewLodState();
    this->setViewportSceneGraphWithLod(root);
    this->syncRenderManager();
    this->requestLodCapacityRender("render-scene-root");
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
    this->requestLodCapacityRender("view-attachment");
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
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
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

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);

    const double perspectiveDegrees =
	bv_perspective_get(view);
    const SbBool wantPerspective = perspectiveDegrees > SMALL_FASTF ?
				   TRUE : FALSE;

    SoCamera *const previousCamera = this->d->activeCamera;
    const SbBool cameraReplaced = !previousCamera ||
	(wantPerspective &&
	 !previousCamera->isOfType(SoPerspectiveCamera::getClassTypeId())) ||
	(!wantPerspective &&
	 !previousCamera->isOfType(SoOrthographicCamera::getClassTypeId()));
    if (cameraReplaced && !createCamera)
	return FALSE;

    const int viewWidth = bv_width_get(view);
    const int viewHeight = bv_height_get(view);
    const SbViewportRegion previousRegion = this->d->viewportRegion;
    SbViewportRegion nextRegion = previousRegion;
    SbVec2s window = previousRegion.getWindowSize();
    if (window[0] <= 1 && window[1] <= 1 &&
	viewWidth > 0 && viewHeight > 0) {
	nextRegion = controller_viewport_region_with_size(previousRegion,
	    static_cast<unsigned int>(viewWidth),
	    static_cast<unsigned int>(viewHeight));
    }
    const bool viewportChanged = nextRegion != previousRegion;

    double aspect = controller_aspect_from_region(nextRegion);
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

    const SbBool orthographic = wantPerspective ? FALSE : TRUE;
    double depthReferenceSize = horizontalSize;
    double nextDepthReferenceSize =
	this->d->orthographicDepthReferenceSize;
    uint64_t nextDepthReferenceRevision =
	this->d->orthographicDepthReferenceStructuralRevision;
    if (orthographic) {
	const uint64_t structuralRevision =
	    this->d->sceneController.getStructuralRevision();
	if (nextDepthReferenceRevision != structuralRevision) {
	    nextDepthReferenceRevision = structuralRevision;
	    nextDepthReferenceSize = horizontalSize;
	} else if (horizontalSize > nextDepthReferenceSize) {
	    nextDepthReferenceSize = horizontalSize;
	}
	depthReferenceSize = std::max(horizontalSize, nextDepthReferenceSize);
    }
    double distance = orthographic ?
		      depthReferenceSize *
			  controller_orthographic_camera_distance_scale :
		      (verticalSize * 0.5) / std::tan(heightAngle * 0.5);
    if (distance <= horizontalSize * 0.5)
	distance = horizontalSize * 0.5;
    if (distance <= SMALL_FASTF)
	distance = controller_camera_minimum_distance;

    const SbRotation orientation =
	bobol_camera_orientation_from_brl_rotation(viewRotation);
    const double viewZ[3] = {
	viewRotation[8], viewRotation[9], viewRotation[10]
    };

    const float desiredAspect = static_cast<float>(aspect);
    const SbVec3f desiredPosition(
	static_cast<float>(center[X] + viewZ[X] * distance),
	static_cast<float>(center[Y] + viewZ[Y] * distance),
	static_cast<float>(center[Z] + viewZ[Z] * distance));
    const float desiredFocal = static_cast<float>(distance);
    const float desiredNear = static_cast<float>(std::max(
	depthReferenceSize * controller_camera_near_scale,
	controller_camera_minimum_near_distance));
    const float desiredFar = static_cast<float>(
	std::max(distance + depthReferenceSize * controller_camera_far_scale,
	distance + controller_camera_minimum_distance));
    const auto float_changed = [](float current, float desired) {
	const float scale = std::max(1.0f, std::fabs(desired));
	return std::fabs(current - desired) > 1.0e-6f * scale;
    };

    const float desiredProjectionHeight = static_cast<float>(
	wantPerspective ? heightAngle : verticalSize);
    bool cameraFieldsChanged = cameraReplaced;
    if (!cameraReplaced) {
	cameraFieldsChanged =
	    previousCamera->viewportMapping.getValue() !=
		SoCamera::LEAVE_ALONE ||
	    float_changed(previousCamera->aspectRatio.getValue(),
		desiredAspect) ||
	    previousCamera->position.getValue() != desiredPosition ||
	    previousCamera->orientation.getValue() != orientation ||
	    float_changed(previousCamera->focalDistance.getValue(),
		desiredFocal) ||
	    float_changed(previousCamera->nearDistance.getValue(), desiredNear) ||
	    float_changed(previousCamera->farDistance.getValue(), desiredFar);
	if (wantPerspective) {
	    const auto *perspectiveCamera =
		static_cast<const SoPerspectiveCamera *>(previousCamera);
	    cameraFieldsChanged = cameraFieldsChanged || float_changed(
		perspectiveCamera->heightAngle.getValue(),
		desiredProjectionHeight);
	} else {
	    const auto *orthographicCamera =
		static_cast<const SoOrthographicCamera *>(previousCamera);
	    cameraFieldsChanged = cameraFieldsChanged || float_changed(
		orthographicCamera->height.getValue(), desiredProjectionHeight);
	}
    }

    controller_configure_render_environment(this->d->viewport);
    const auto lights = controller_camera_lights(this->d->viewport);
    std::array<SbVec3f, controller_camera_light_count> lightDirections;
    std::array<bool, controller_camera_light_count> lightChanges{{
	false, false, false}};
    if (this->d->headlightEnabled && this->d->headlightCameraTracked) {
	lightDirections = controller_camera_light_directions(orientation,
	    this->d->headlightOffsetEye);
	for (size_t i = 0; i < lights.size(); ++i)
	    lightChanges[i] = lights[i] &&
		lights[i]->direction.getValue() != lightDirections[i];
    }

    SbPlane minimumPlane;
    SbPlane maximumPlane;
    if (!controller_camera_relative_clip_planes(center, viewZ,
	    horizontalSize, this->d->clipMinimum, this->d->clipMaximum,
	    minimumPlane, maximumPlane))
	return FALSE;
    const std::array<SbPlane, BObolViewController::CLIP_PLANE_CAPACITY>
	clipPlanes{{
	minimumPlane,
	maximumPlane,
	this->d->cuttingPlane
    }};
    const SbBool clipping = bv_zclip_get(view) ? TRUE : FALSE;
    const std::array<SbBool, BObolViewController::CLIP_PLANE_CAPACITY>
	clipEnabled{{
	clipping, clipping, this->d->cuttingPlaneEnabled
    }};
    const std::array<SoClipPlane *, BObolViewController::CLIP_PLANE_CAPACITY>
	clipNodes{{
	controller_clip_plane(this->d->viewport, TRUE),
	controller_clip_plane(this->d->viewport, FALSE),
	controller_cutting_plane(this->d->viewport)
    }};
    std::array<bool, BObolViewController::CLIP_PLANE_CAPACITY> clipChanges{{
	false, false, false}};
    for (size_t i = 0; i < clipNodes.size(); ++i)
	clipChanges[i] = clipNodes[i] &&
	    (clipNodes[i]->plane.getValue() != clipPlanes[i] ||
	    clipNodes[i]->on.getValue() != clipEnabled[i]);

    const auto double_changed = [](double current, double desired) {
	const double scale = std::max(1.0, std::fabs(desired));
	return std::fabs(current - desired) > 1.0e-12 * scale;
    };
    const bool affordanceGeometryChanged = viewportChanged ||
	cameraFieldsChanged || double_changed(
	    this->d->cuttingPlaneAffordanceHorizontalSize, horizontalSize) ||
	double_changed(this->d->cuttingPlaneAffordanceAspect, aspect);
    const bool affordanceChanged =
	controller_cutting_plane_affordance_update_needed(
	    this->d->framebufferOverlayRoot,
	    this->d->cuttingPlaneEnabled, affordanceGeometryChanged);
    const bool lightChanged = std::any_of(lightChanges.begin(),
	lightChanges.end(), [](bool changed) { return changed; });
    const bool clipChanged = std::any_of(clipChanges.begin(),
	clipChanges.end(), [](bool changed) { return changed; });
    const bool changed = viewportChanged || cameraFieldsChanged ||
	lightChanged || clipChanged || affordanceChanged;
    if (!changed)
	return TRUE;

    /* Keep the installed camera identity for ordinary navigation.  The
     * detached camera is the complete scalar successor and also supplies the
     * view volume used to build the section aid before any live field moves. */
    SbModernUtils::SoNodeRef cameraCandidateOwner(
	cameraFieldsChanged ? (wantPerspective ?
	    static_cast<SoNode *>(new SoPerspectiveCamera) :
	    static_cast<SoNode *>(new SoOrthographicCamera)) : NULL);
    SoCamera *cameraCandidate = cameraFieldsChanged ?
	static_cast<SoCamera *>(cameraCandidateOwner.get()) : previousCamera;
    if (cameraFieldsChanged) {
	cameraCandidate->viewportMapping = SoCamera::LEAVE_ALONE;
	cameraCandidate->aspectRatio = desiredAspect;
	cameraCandidate->position = desiredPosition;
	cameraCandidate->orientation = orientation;
	cameraCandidate->focalDistance = desiredFocal;
	cameraCandidate->nearDistance = desiredNear;
	cameraCandidate->farDistance = desiredFar;
	if (wantPerspective)
	    static_cast<SoPerspectiveCamera *>(cameraCandidate)->heightAngle =
		desiredProjectionHeight;
	else
	    static_cast<SoOrthographicCamera *>(cameraCandidate)->height =
		desiredProjectionHeight;
    }

    PreparedScalarNodes scalarFields;
    scalarFields.reserve(1 + controller_camera_light_count +
	BObolViewController::CLIP_PLANE_CAPACITY);
    if (cameraFieldsChanged && !cameraReplaced)
	scalarFields.prepare(*previousCamera, *cameraCandidate);

    std::array<SbModernUtils::SoNodeRef, controller_camera_light_count>
	lightCandidateOwners{{
	    SbModernUtils::SoNodeRef(NULL),
	    SbModernUtils::SoNodeRef(NULL),
	    SbModernUtils::SoNodeRef(NULL)}};
    for (size_t i = 0; i < lights.size(); ++i) {
	if (!lightChanges[i])
	    continue;
	lightCandidateOwners[i] =
	    SbModernUtils::SoNodeRef(new SoDirectionalLight);
	auto *candidate = static_cast<SoDirectionalLight *>(
	    lightCandidateOwners[i].get());
	copy_publication_scalar_fields(*candidate, *lights[i]);
	candidate->direction = lightDirections[i];
	scalarFields.prepare(*lights[i], *candidate);
    }

    std::array<SbModernUtils::SoNodeRef,
	BObolViewController::CLIP_PLANE_CAPACITY> clipCandidateOwners{{
	SbModernUtils::SoNodeRef(NULL),
	SbModernUtils::SoNodeRef(NULL),
	SbModernUtils::SoNodeRef(NULL)}};
    for (size_t i = 0; i < clipNodes.size(); ++i) {
	if (!clipChanges[i])
	    continue;
	clipCandidateOwners[i] = SbModernUtils::SoNodeRef(new SoClipPlane);
	auto *candidate = static_cast<SoClipPlane *>(
	    clipCandidateOwners[i].get());
	copy_publication_scalar_fields(*candidate, *clipNodes[i]);
	candidate->plane = clipPlanes[i];
	candidate->on = clipEnabled[i];
	scalarFields.prepare(*clipNodes[i], *candidate);
    }

    std::unique_ptr<BObolPreparedCuttingPlaneAffordance> affordance;
    if (affordanceChanged)
	affordance =
	    std::make_unique<BObolPreparedCuttingPlaneAffordance>(
		this->d->viewport, this->d->framebufferOverlayRoot,
		cameraCandidate, this->d->cuttingPlane,
		this->d->cuttingPlaneEnabled, horizontalSize, aspect);

    std::unique_ptr<SoViewport::CameraReplacement> cameraReplacement;
    if (cameraReplaced) {
	cameraReplacement = this->d->viewport->prepareCameraReplacement(
	    cameraCandidate, controller_camera_root_index(
		this->d->viewport, previousCamera));
    }

    BObolPreparedRenderRequest renderRequest = this->prepareRenderRequest(
	"rt-view-camera", RenderRequestIntent::LOD_CAPACITY);

    if (cameraFieldsChanged || viewportChanged) {
	/* LoD snapshot preparation reads controller accessors.  Expose the
	 * detached camera and candidate region only while that fallible work is
	 * performed, then restore the public predecessor until commit. */
	this->d->activeCamera = cameraCandidate;
	this->d->viewportRegion = nextRegion;
	this->d->viewport->setViewportRegion(nextRegion);
	this->d->renderManager->setViewportRegion(nextRegion);
	try {
	    this->syncLodViewSignature(TRUE, FALSE);
	} catch (...) {
	    this->d->activeCamera = previousCamera;
	    this->d->viewportRegion = previousRegion;
	    this->d->viewport->setViewportRegion(previousRegion);
	    this->d->renderManager->setViewportRegion(previousRegion);
	    throw;
	}
	this->d->activeCamera = previousCamera;
	this->d->viewportRegion = previousRegion;
	this->d->viewport->setViewportRegion(previousRegion);
	this->d->renderManager->setViewportRegion(previousRegion);
    }

    /* Every operation below is allocation-free.  Commit all public copies and
     * derived nodes before releasing the first observer. */
    if (cameraReplaced)
	cameraCandidate->ref();
    scalarFields.commit();
    if (affordance)
	affordance->commit();
    if (cameraReplacement)
	cameraReplacement->commit();
    this->d->activeCamera = cameraReplaced ? cameraCandidate : previousCamera;
    this->d->viewportRegion = nextRegion;
    this->d->viewport->setViewportRegion(nextRegion);
    this->d->renderManager->setViewportRegion(nextRegion);
    if (cameraReplaced)
	this->d->renderManager->setCamera(cameraCandidate);
    this->d->lastCameraOrientation = orientation;
    this->d->orthographicDepthReferenceSize = nextDepthReferenceSize;
    this->d->orthographicDepthReferenceStructuralRevision =
	nextDepthReferenceRevision;
    this->d->cuttingPlaneAffordanceHorizontalSize = horizontalSize;
    this->d->cuttingPlaneAffordanceAspect = aspect;
    this->commitRenderRequest(renderRequest);

    /* Reentrant observer writes must notify normally, so restore every quiet
     * scalar participant before delivering any graph or field callback. */
    scalarFields.restore();

    std::exception_ptr failure;
    if (cameraReplacement) {
	try { cameraReplacement->notify(); }
	catch (...) { failure = std::current_exception(); }
    }
    if (affordance)
	affordance->notify(failure);
    scalarFields.notify(failure);
    try { this->notifyRenderRequest(renderRequest); }
    catch (...) { if (!failure) failure = std::current_exception(); }
    if (cameraReplaced && previousCamera)
	previousCamera->unref();
    if (changedOut)
	*changedOut = TRUE;
    if (failure)
	std::rethrow_exception(failure);
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
    BObolLodControlTransitionScope controlTransition(this);
    struct RenderEffects : BObolSourceRealizationEffects {
	explicit RenderEffects(BObolViewController &controller) : view(controller) {}
	void prepare(const SoBRLDatabaseSource &, bool realized, const SbString &,
	    const std::vector<SoNode *> &, const std::vector<SoNode *> &) override
	{
	    const bool failed = !realized || view.getLastFailedSourceCount() > 0;
	    request = view.prepareRenderRequest(failed ? "realize-failed" : "realize",
		RenderRequestIntent::LOD_CAPACITY);
	}
	void commit(bool changed) noexcept override
	{
	    if (!changed) return;
	    view.commitRenderRequest(request);
	    sourceCommitted = true;
	}
	void notify() override { view.notifyRenderRequest(request); }
	BObolViewController &view;
	BObolPreparedRenderRequest request;
	bool sourceCommitted = false;
    } effects(*this);
    const SbBool ret = this->d->sceneController.realizePending(&effects);
    /* Keep the explicit API's repaint when no source changed. Otherwise the
     * source commits already published their requests and notified the host. */
    if (!effects.sourceCommitted)
	this->requestLodCapacityRender(ret ? "realize" : "realize-failed");
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
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    const int requested = enabled ? 1 : 0;

    const bool changed = this->d->endpointGraphicalRenderingEnabled.exchange(
	requested, std::memory_order_acq_rel) != requested;
    if (!changed &&
	!this->d->endpointGraphicalRenderingSyncPending.load(
	    std::memory_order_acquire))
	return;
    try {
	this->syncRenderManager();
	this->d->endpointGraphicalRenderingSyncPending.store(
	    0, std::memory_order_release);
    } catch (...) {
	this->d->endpointGraphicalRenderingSyncPending.store(
	    1, std::memory_order_release);
	throw;
    }
    if (!changed)
	return;
    if (enabled)
	this->requestLodCapacityRender("render-engine");
    else
	this->clearRenderRequest();
}

void
BObolViewController::invalidateRendererPerformanceHistory(void)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    PreparedRendererInvalidation publication;
    this->prepareRendererInvalidation(publication);
    this->commitRendererInvalidation(publication, TRUE);
    this->notifyRendererInvalidation(publication);
}

void
BObolViewController::prepareRendererInvalidation(
    PreparedRendererInvalidation &publication) const
{
    this->prepareRendererInvalidation(publication, "renderer-performance");
}

void
BObolViewController::prepareRendererInvalidation(
    PreparedRendererInvalidation &publication, const char *reason) const
{
    publication.automatic = this->automaticLodControlEnabled() != FALSE;
    publication.renderRequest = this->prepareRenderRequest(
	reason, publication.automatic ?
	RenderRequestIntent::LOD_CAPACITY : RenderRequestIntent::PRESENTATION);
}

void
BObolViewController::commitRendererInvalidation(
    PreparedRendererInvalidation &publication, SbBool requestFrame) noexcept
{
    /* Renderer capacity is part of an exact-view quality proof even though it
	* is deliberately not part of camera identity.  Retain immutable scene
	* residency and the coherent currently presented cut, but invalidate every
	* timing-derived proof and arrange one measured-frame rescan. */
    this->d->resetRendererPerformanceEvidence();

    if (publication.automatic) {
	const BObolViewLodState *presentationState =
	    this->d->viewAttachment ?
		this->d->viewAttachment->getViewLodState() : NULL;
	if (presentationState &&
	    presentationState->hasCadPresentationAssemblies())
	    this->d->requireExactPresentationFrame();
	this->publishProgressiveWorkPending(FALSE);
    } else {
	/* resetRendererPerformanceEvidence() creates a fresh capacity-search
	 * certificate.  Policy-off has no consumer for it, so retire that
	 * automatic debt before requesting the ordinary style repaint. */
	this->retireAutomaticLodControl();
    }
    if (requestFrame) {
	this->commitRenderRequest(publication.renderRequest);
	publication.requestCommitted = true;
    }
}

void
BObolViewController::notifyRendererInvalidation(
    const PreparedRendererInvalidation &publication)
{
    if (publication.requestCommitted)
	this->notifyRenderRequest(publication.renderRequest);
}

void
BObolViewController::requestLodCapacityRender(const char *reason)
{
    this->requestRenderImpl(reason, RenderRequestIntent::LOD_CAPACITY);
}

void
BObolViewController::requestPresentationRender(const char *reason)
{
    this->requestRenderImpl(reason, RenderRequestIntent::PRESENTATION);
}

void
BObolViewController::requestLodPresentationRender(const char *reason)
{
    this->requestRenderImpl(reason, RenderRequestIntent::LOD_PLANNING);
}

void
BObolViewController::requestExactCadPresentationRender(const char *reason)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    const BObolViewLodState *presentationState =
	this->d->viewAttachment ?
	    this->d->viewAttachment->getViewLodState() : NULL;
    if (this->automaticLodControlEnabled() && presentationState &&
	presentationState->hasCadPresentationAssemblies())
	this->d->requireExactPresentationFrame();
    this->requestRenderImpl(reason, RenderRequestIntent::PRESENTATION);
}

void
BObolViewController::retireLodCapacityRenderRequest(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    if (!this->d->renderRequest.retireCapacity())
	return;

    /* The repaint remains necessary, but LoD-off has no capacity contract
     * which could consume its timing.  Preserve the existing host wakeup and
     * downgrade only its evidence class. */
    bobol_identity_advance(this->d->hostWorkRevision);
    bobol_identity_advance(this->d->renderRequestSerial);
}

void
BObolViewController::requestRenderImpl(const char *reason,
	RenderRequestIntent intent)
{
    BObolLodControlTransitionScope controlTransition(this);
    auto request = this->prepareRenderRequest(reason, intent);
    this->commitRenderRequest(request);
    this->notifyRenderRequest(request);
}

BObolPreparedRenderRequest
BObolViewController::prepareRenderRequest(const char *reason,
	RenderRequestIntent intent) const
{
    /* The active policy classifies both ordinary requests and requests prepared
     * for source publication. LoD-off must not recreate an automatic owner. */
    const bool lodEnabled = this->lodViewPolicyEnabled();
    return BObolPreparedRenderRequest(reason,
	lodEnabled && intent == RenderRequestIntent::LOD_CAPACITY,
	lodEnabled && intent != RenderRequestIntent::PRESENTATION);
}

void
BObolViewController::commitRenderRequest(BObolPreparedRenderRequest &request) noexcept
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    /* Merge against the current level, so a request arriving after preparation
     * cannot be overwritten by a weaker prepared request. */
    request.decision = this->d->renderRequest.request(std::move(request.reason),
	request.capacityRelevant, request.planningRelevant);
    if (this->d->renderRequest.pending())
	this->d->lodExactPresentationFrame.noteFrameRequested();
    if (request.decision.changed) {
	bobol_identity_advance(this->d->hostWorkRevision);
	bobol_identity_advance(this->d->renderRequestSerial);
	request.serial = this->d->renderRequestSerial;
    }
}

void
BObolViewController::notifyRenderRequest(const BObolPreparedRenderRequest &request)
{
    const char *reason = request.notificationReason.getString();
    std::exception_ptr failure;
    try {
	if (request.decision.changed && controller_lod_trace_enabled(
	    "BOBOL_LOD_TRACE_RENDER_REQUEST", this->d->lodViewRevision.value()))
	    bu_log("BObol LoD render request serial=%llu view=%llu "
		"policy=%llu reason=%s progressive=%d planning=%d capacity=%d\n",
		static_cast<unsigned long long>(request.serial),
		static_cast<unsigned long long>(this->d->lodViewRevision.value()),
		static_cast<unsigned long long>(this->d->lodPolicyRevision.value()),
		reason, this->hasProgressiveWorkPending() ? 1 : 0,
		request.planningRelevant ? 1 : 0, request.capacityRelevant ? 1 : 0);
    } catch (...) { failure = std::current_exception(); }
    /* Diagnostics and source observers cannot prevent the host wakeup. The
     * request level is already complete, even when the host callback throws. */
    try {
	if (request.decision.wakeEndpoint)
	    this->notifyFrameRequest(reason);
    } catch (...) { if (!failure) failure = std::current_exception(); }
    if (failure) std::rethrow_exception(failure);
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
    if (previous && previous->close())
	delete previous;

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
    if (state->close())
	delete state;
}

void
BObolViewController::setPresentationSyncCallback(
    BObolPresentationSyncCallback callback, void *userData)
{
    ControllerPresentationSyncState *replacement = callback ?
	new (std::nothrow) ControllerPresentationSyncState(
	    callback, userData) : NULL;
    if (callback && !replacement)
	return;

    ControllerPresentationSyncState *previous = NULL;
    {
	std::lock_guard<std::mutex> lock(this->d->presentationSyncMutex);
	previous = static_cast<ControllerPresentationSyncState *>(
	    this->d->presentationSyncUserData);
	this->d->presentationSyncCallback = callback;
	this->d->presentationSyncUserData = replacement;
    }
    if (previous && previous->close())
	delete previous;
}

void
BObolViewController::clearPresentationSyncCallback(void *userData)
{
    ControllerPresentationSyncState *state = NULL;
    {
	std::lock_guard<std::mutex> lock(this->d->presentationSyncMutex);
	state = static_cast<ControllerPresentationSyncState *>(
	    this->d->presentationSyncUserData);
	if (!state || state->userData != userData)
	    return;
	this->d->presentationSyncCallback = NULL;
	this->d->presentationSyncUserData = NULL;
    }
    if (state->close())
	delete state;
}

void
BObolViewController::synchronizePresentation(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    const uint64_t started = this->beginRenderTiming();
    BObolPresentationSyncCallback callback = NULL;
    void *userData = NULL;
    ControllerPresentationSyncState *state = NULL;
    {
	std::lock_guard<std::mutex> lock(this->d->presentationSyncMutex);
	state = static_cast<ControllerPresentationSyncState *>(
	    this->d->presentationSyncUserData);
	if (this->d->presentationSyncCallback && state &&
	    state->beginDispatch()) {
	    callback = state->callback;
	    userData = state->userData;
	}
    }
    if (callback) {
	ControllerCallbackDispatchScope<ControllerPresentationSyncState>
	    dispatchScope(state);
	(*callback)(userData);
    }
    controller_synchronize_compact_cad_presentations(this);
    /* Presentation-only hierarchy edits do not enqueue provider work.  They
     * still change the exact LoD visibility denominator, so rendering their
     * retained instance delta must wake the bounded planner which consumes
     * the source journal.  Keep this level-triggered: a frame may coalesce
     * the original GED mutation notification, and an already-consumed source
     * revision is an O(source-count) no-op. */
    if (this->automaticLodControlEnabled() &&
	controller_lod_source_inputs_unsubmitted(
	    controller_render_database_source_roots(this),
	    this->d->lodSourceEvidence.submitted())) {
	(void)this->publishPendingLodSourceRevision();
	this->markProgressiveWorkPending();
    }
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

    ControllerCallbackDispatchScope<ControllerFrameRequestState>
	dispatchScope(state);
    (*callback)(userData, reason ? reason : "");
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
    if (this->d->renderRequest.pending()) {
	snapshot.flags |= BOBOL_HOST_WORK_RENDER;
	if (this->d->renderRequest.capacityRelevant())
	    snapshot.flags |= BOBOL_HOST_WORK_CAPACITY_SAMPLE;
    }
    if (this->d->renderClaimed)
	snapshot.flags |= BOBOL_HOST_WORK_FRAME_CLAIMED;
    if (this->d->capacitySampleClaimed)
	snapshot.flags |= BOBOL_HOST_WORK_CAPACITY_SAMPLE_CLAIMED;
    return snapshot;
}

void
BObolViewController::clearRenderRequest(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    SbBool exactRequestRetired = FALSE;
    {
	std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
	const SbBool changed = this->d->renderRequest.clear() ? TRUE : FALSE;
	if (changed) {
	    bobol_identity_advance(this->d->hostWorkRevision);
	    bobol_identity_advance(this->d->renderRequestSerial);
	    if (this->d->lodExactPresentationFrame.framePending()) {
		this->d->lodExactPresentationFrame.noteRequestRetired();
		exactRequestRetired = TRUE;
	    }
	}
    }
    if (exactRequestRetired)
	this->markProgressiveWorkPending();
}

SbBool
BObolViewController::consumeRenderRequest(SbString *reason,
	SbBool *lodCapacityRelevant, SbBool *lodPlanningRelevant)
{
    BObolLodControlTransitionScope controlTransition(this);
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    const bool requestCapacityRelevant =
	this->d->renderRequest.capacityRelevant();
    const SbBool ret = this->d->renderRequest.consume(
	reason, lodCapacityRelevant, lodPlanningRelevant) ? TRUE : FALSE;
    if (ret) {
	/* The request level disappears when it is claimed, but the capacity
	 * transaction remains foreground work until completed-frame timing
	 * accepts or rejects it.  Keep that classification independently of
	 * whether the caller requested the optional output. */
	this->d->capacitySampleClaimed =
	    this->d->capacitySampleClaimed || requestCapacityRelevant;
	this->d->renderClaimed = TRUE;
	bobol_identity_advance(this->d->hostWorkRevision);
	bobol_identity_advance(this->d->renderRequestSerial);
    }
    return ret;
}

void
BObolViewController::retireClaimedRender(void)
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    if (!this->d->renderClaimed)
	return;
    this->d->renderClaimed = FALSE;
    this->d->capacitySampleClaimed = FALSE;
    bobol_identity_advance(this->d->hostWorkRevision);
}

void
BObolViewController::retireDisplayEndpointWork(void)
{
    try {
	BObolLodControlTransitionScope controlTransition(
	    this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
	this->d->retireDisplayEndpointWorkNoexcept();
    } catch (...) {
	/* Endpoint loss must still cancel every owner when diagnostic transition
	 * setup or an observer rejects the live path. */
	this->d->retireDisplayEndpointWorkNoexcept();
    }

    BU_ASSERT(this->d->automaticLodControlRetired() &&
	this->getHostWorkSnapshot().flags == BOBOL_HOST_WORK_NONE);
}

void
BObolViewController::resumeDisplayEndpointWork(void)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);

    /* Retirement preserves immutable residency and the user's policy.  A new
     * endpoint therefore starts a fresh coverage/capacity transaction while
     * retaining reusable payloads.  Independent providers are included by
     * the shared level projection below. */
    (void)this->synchronizeAutomaticLodControl();
    this->synchronizeProgressiveWorkPending();
}

uint64_t
BObolViewController::renderRequestSerialGet(void) const
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    return this->d->renderRequestSerial;
}
