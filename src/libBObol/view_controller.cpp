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
#include "raytrace.h"
#include "rt/view.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <sstream>
#include <string.h>
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
     * bounded by maxLodApplyMicroseconds, so a larger cap lets a big assembly
     * stream in far faster without risking frame-time spikes (the microsecond
     * budget short-circuits the batch loop). */
    maxLodResults(256),
    maxLodApplyMicroseconds(4000),
    maxProviders(0),
    maxProviderItems(64),
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

struct BObolViewController::Impl {
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
    uint64_t lastRenderTimeNanoseconds = 0;
    uint64_t smoothedRenderTimeNanoseconds = 0;
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
    SbBool lodAutoSubmit = FALSE;
    uint64_t lodActiveGeneration = 0;
    size_t lodSubmissionSourceIndex = 0;
    size_t lodSubmissionEntryOffset = 0;
    SoBRLDatabaseSource *lodSubmissionPlanSource = NULL;
    std::vector<size_t> lodSubmissionPlanEntries;
    SbBool lodSubmissionPlanValid = FALSE;
    void clearLodSubmissionPlan(void)
    {
	lodSubmissionPlanSource = NULL;
	lodSubmissionPlanEntries.clear();
	lodSubmissionPlanValid = FALSE;
    }
    SbBool lodSubmissionPending = FALSE;
    SbBool lodSubmissionRescanPending = FALSE;
    SbBool lodRetainedRefinementPending = FALSE;
    SbBool lodRetainedRefinementCutAdvanced = FALSE;
    SbBool lodRetainedRefinementBudgetBlocked = FALSE;
    SbBool lodBudgetRescanAfterFrame = FALSE;
    SbBool lodPassAdmittedWork = FALSE;
    SbBool forceTerminalLodRefinement = FALSE;
    SbBool lodRefinementFaceBudgetInitialized = FALSE;
    size_t lodRefinementFaceBudgetRemaining = SIZE_MAX;
    SbBool lodRetainedAdmissionActive = FALSE;
    size_t lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
    int64_t lodRefinementBudgetRetryMicroseconds = 0;
    SbBool lodSubmissionRefreshMissing = TRUE;
    int lodSubmissionReset = 0;
    uint64_t lodLastSubmittedViewRevision = 0;
    uint64_t lodLastSubmittedPolicyRevision = 0;
    SbString lodLastSubmittedSourceSignature = SbString("");
    SbString lodViewSignature = SbString("");
    uint64_t lodViewRevision = 1;
    uint64_t lodPolicyRevision = 1;
    int64_t lodLastViewChangeMicroseconds = 0;
    SbBool lodInteractive = FALSE;
    SbBool lodGestureActive = FALSE;
    uint64_t lodSettleAfterRenderSerial = 0;
    SbBool lodResidentCompactionPending = FALSE;
    SbBool lodRefinementAwaitingFrame = FALSE;
    uint64_t lodRefinementResumeAfterRenderSerial = 0;
    int64_t lodRefinementNotBeforeMicroseconds = 0;
    float lodTargetPixelError = 1.0f;
    float lodInteractiveTargetFps = 60.0f;
    float lodStableTargetFps = 20.0f;
    size_t lodSeedFaceBudget = 50000;
    size_t lodCurrentFaceBudget = 50000;
    long double lodInteractiveCalibratedFacesPerSecond = 0.0L;
    long double lodStableCalibratedFacesPerSecond = 0.0L;
    int lodEmergencyProgressiveCeiling = -1;
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
    gl->clear(GL_DEPTH_BUFFER_BIT);
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
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    this->d->renderRequested = TRUE;
    this->d->renderReason = reason ? reason : "";
    this->d->renderRequestSerial++;
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
    this->d->renderManager->render(static_cast<SbBool>(FALSE), static_cast<SbBool>(FALSE));
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
    const size_t calibrationFaces = calibrationState ?
	calibrationState->activeFaceCount() : 0;
    uint64_t calibrationElapsed = elapsed;
    if (this->d->lodGestureActive) {
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
    if (calibrationFaces > 0 && calibrationElapsed > 0) {
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

    /* A budget-limited pass may have admitted a bounded first wave and then
     * scanned the remaining boxes without admitting them.  Re-plan from the
     * highest screen-value occurrence after that wave has actually rendered
     * and supplied a throughput sample.  Replanning before the frame would
     * repeatedly spend the same stale allowance. */
    if (this->d->lodBudgetRescanAfterFrame) {
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
	this->markProgressiveWorkPending();
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
	 * useful cut visible for roughly four render times (bounded at two
	 * seconds) before exposing the next, potentially much larger prefix.
	 * This is feedback from an actually presented frame, so it adapts to
	 * OSMesa, System GL, viewport size, and current scene complexity without
	 * backend-specific policy. */
	const uint64_t responsiveFrame = 33333334ULL;
	int64_t cooldownMicroseconds = 0;
	if (elapsed > responsiveFrame) {
	    const uint64_t elapsedMicroseconds = elapsed / 1000ULL;
	    const uint64_t scaled =
		elapsedMicroseconds > 500000ULL ?
		2000000ULL : elapsedMicroseconds * 4ULL;
	    cooldownMicroseconds = static_cast<int64_t>(
		std::max<uint64_t>(50000ULL,
		    std::min<uint64_t>(2000000ULL, scaled)));
	}
	const int64_t nowMicroseconds = bu_gettime();
	this->d->lodRefinementNotBeforeMicroseconds =
	    cooldownMicroseconds > 0 &&
	    nowMicroseconds <=
		std::numeric_limits<int64_t>::max() - cooldownMicroseconds ?
	    nowMicroseconds + cooldownMicroseconds : nowMicroseconds;
	this->d->lodLastSubmittedViewRevision = 0;
	this->d->lodLastSubmittedPolicyRevision = 0;
	this->markProgressiveWorkPending();
    }
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

	    /* The control EMA above is deliberately quick: it lets the LoD
	     * policy observe a change in renderer capacity within a few frames.
	     * That same response makes a faceplate sampled four times per second
	     * look needlessly nervous.  Give the displayed cadence a 750 ms
	     * elapsed-time constant.  Computing alpha from elapsed time, rather
	     * than using a fixed sample weight, gives the same visual response at
	     * 10, 60, or 240 FPS. */
	    uint64_t &displayed =
		this->d->displayedPresentationIntervalNanoseconds;
	    if (!displayed) {
		displayed = interval;
	    } else {
		constexpr double timeConstantNanoseconds = 750000000.0;
		const double alpha = -std::expm1(
		    -static_cast<double>(interval) / timeConstantNanoseconds);
		const double next = static_cast<double>(displayed) +
		    alpha * (static_cast<double>(interval) -
			static_cast<double>(displayed));
		displayed = static_cast<uint64_t>(
		    std::max(1.0, std::floor(next + 0.5)));
	    }
	} else if (interval > 250000000ULL) {
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
	    }
	}
	this->d->lastInteractivePresentationTimestampNanoseconds = now;
    } else {
	this->d->lastInteractivePresentationTimestampNanoseconds = 0;
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
	    this->d->lodInteractive = FALSE;
	    this->d->lodTargetPixelError = 1.0f;
	    this->d->lodEmergencyProgressiveCeiling = -1;
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    if (viewState) {
		viewState->setCadPresentationProgressiveLodCeiling(-1);
		/* The motion cut is already a measured, coherent presentation.
		 * Do not begin the quiet epoch with a smaller stale face budget:
		 * retained admission would replace a tail of minimum meshes with
		 * boxes on button release, only to restore them moments later.
		 * Stable refinement starts from the complete cut the user just
		 * saw and is monotonic from there. */
		if (this->d->lodCurrentFaceBudget != SIZE_MAX)
		    this->d->lodCurrentFaceBudget = std::max(
			this->d->lodCurrentFaceBudget,
			viewState->activeFaceCount());
	    }
	    this->advanceLodPolicyRevision();
	} else {
	    /* Keep the host frame pump alive through the idle debounce. */
	    localStatus.hasMore = 1;
	}
    }

    if (this->hasPendingLodResults() ||
	(this->d->lodService &&
	 this->d->lodService->queuedResultCountForDiagnostics() > 0)) {
	const uint64_t resultStarted = this->beginRenderTiming();
	(void)this->processPendingLodResults(options->maxLodResults,
	    options->maxLodApplyMicroseconds);
	const uint64_t resultCompleted = this->beginRenderTiming();
	this->d->lastLodResultProcessingTimeNanoseconds =
	    resultCompleted > resultStarted ? resultCompleted - resultStarted : 0;
	localStatus.lodResultsProcessed = this->d->lastLodResultCount;
	localStatus.lodResultsApplied = this->d->lastLodAppliedResultCount;
	if (this->d->lastLodAppliedResultCount > 0)
	    localStatus.changed = 1;
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
    }

    const int pending_service_work =
	(localStatus.pendingTasks > 0 || localStatus.inFlight > 0 ||
	 localStatus.queuedResults > 0 || localStatus.queuedCacheWrites > 0) ?
	1 : 0;
    if (pending_service_work)
	localStatus.hasMore = 1;

    /* Refinement and reclamation are separate phases.  A quiet view first
     * reaches its 1 px display target; only after a longer quiet interval and
     * an idle service do we replace this consumer's complete demand snapshot
     * and trim shared CPU prefixes to the aggregate maximum. */
    if (this->d->lodAutoSubmit &&
	this->d->lodResidentCompactionPending &&
	!this->d->lodInteractive &&
	this->d->lodLastViewChangeMicroseconds > 0) {
	const int64_t elapsed =
	    bu_gettime() - this->d->lodLastViewChangeMicroseconds;
	if (elapsed < 750000 || pending_service_work ||
	    this->d->lodSubmissionPending) {
	    localStatus.hasMore = 1;
	} else if (this->d->lodService) {
	    std::vector<BObolLodResidentDemand> demands;
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    if (viewState)
		viewState->residentMeshDemands(demands);
	    const size_t compacted =
		this->d->lodService->compactResidentMeshes(
		    static_cast<uint64_t>(
			reinterpret_cast<uintptr_t>(this)),
		    this->d->lodViewRevision, demands);
	    this->d->lodResidentCompactionPending = FALSE;
	    if (compacted) {
		if (viewState)
		    viewState->noteResidentMeshesChanged();
		localStatus.changed = 1;
	    }
	}
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
    if (localStatus.changed)
	this->requestRender("progressive-update");

    if (status)
	*status = localStatus;

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

static void
append_signature_string(std::ostringstream &out, const char *value)
{
    size_t len = value ? strlen(value) : 0;
    out << len << ':' << (value ? value : "") << ';';
}

static SbString
controller_lod_view_signature(const struct bv_view_info &view,
			      SbBool haveCamera,
			      SoCamera *camera,
			      const SbViewportRegion &region)
{
    std::ostringstream out;

    out << (haveCamera ? 1 : 0) << ';'
	<< view.width << ';'
	<< view.height << ';'
	<< std::setprecision(17)
	<< view.size << ';'
	<< view.lod.scale << ';'
	<< view.lod.curve_scale << ';'
	<< view.lod.point_scale << ';'
	<< view.lod.bot_threshold << ';';
    if (haveCamera && camera) {
	double aspect = controller_aspect_from_region(region);
	if (aspect <= SMALL_FASTF)
	    aspect = 1.0;
	const SbMatrix matrix = camera->getViewVolume(
	    static_cast<float>(aspect)).getMatrix();
	const float *values = matrix[0];
	for (size_t i = 0; i < 16; i++)
	    out << std::setprecision(9) << values[i] << ';';
    }

    return SbString(out.str().c_str());
}

static SbString
controller_lod_source_signature(const BObolViewController *controller)
{
    std::ostringstream out;

    if (!controller)
	return SbString("");

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_sources(controller);
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (!source || !source->getDatabase())
	    continue;
	if (source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    source->needsRealization() ||
	    (!source->hasRealizedMeshGeometry() &&
	     !source->hasDisplayMeshLodRequests()))
	    continue;

	out << "source;";
	append_signature_string(out, controller_database_id(source->getDatabase()));
	append_signature_string(out, source->path.getValue().getString());
	out << source->drawMode.getValue() << ';'
	    << source->visible.getValue() << ';'
	    << source->lodBotThreshold.getValue() << ';'
	    << source->sourceRevision.getValue() << ';'
	    << source->inputsRevision.getValue() << ';'
	    << source->realizedSourceRevision.getValue() << ';'
	    << source->realizedInputsRevision.getValue() << ';'
	    << source->realizedViewRevision.getValue() << ';'
	    << source->getDisplayMeshLodRevision() << ';';
    }

    return SbString(out.str().c_str());
}

void
BObolViewController::lodResultReadyCB(
    BObolLodService *UNUSED(service), void *userData)
{
    BObolViewController *controller =
	static_cast<BObolViewController *>(userData);
    if (controller) {
	controller->d->lodResultsPending.store(1);
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
    this->d->lodPassAdmittedWork = FALSE;
    this->d->lodRefinementFaceBudgetInitialized = FALSE;
    this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
    this->d->lodRetainedAdmissionActive = FALSE;
    this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
    this->d->lodRefinementBudgetRetryMicroseconds = 0;
    this->d->lodSettleAfterRenderSerial = 0;
    this->d->lodLastSubmittedViewRevision = 0;
    this->d->lodLastSubmittedPolicyRevision = 0;
    this->d->lodLastSubmittedSourceSignature = "";
    this->d->lodResidentCompactionPending = service ? TRUE : FALSE;
    this->d->lodRefinementAwaitingFrame = FALSE;
    this->d->lodRefinementResumeAfterRenderSerial = 0;
    this->d->lodRefinementNotBeforeMicroseconds = 0;

    if (this->d->lodService)
	this->d->lodResultSubscriberId =
	    this->d->lodService->subscribeResultReady(
		BObolViewController::lodResultReadyCB, this);
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
    this->d->lodPassAdmittedWork = FALSE;
    this->d->lodRefinementFaceBudgetInitialized = FALSE;
    this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
    this->d->lodRetainedAdmissionActive = FALSE;
    this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
    this->d->lodRefinementBudgetRetryMicroseconds = 0;
    this->d->lodSettleAfterRenderSerial = 0;
    this->d->lodResultsPending.store(0);
    this->d->lodLastSubmittedViewRevision = 0;
    this->d->lodLastSubmittedPolicyRevision = 0;
    this->d->lodLastSubmittedSourceSignature = "";
    this->d->lodRefinementAwaitingFrame = FALSE;
    this->d->lodRefinementResumeAfterRenderSerial = 0;
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodEmergencyProgressiveCeiling = -1;
    if (this->d->viewAttachment &&
	this->d->viewAttachment->getViewLodState())
	this->d->viewAttachment->getViewLodState()->
	    setCadPresentationProgressiveLodCeiling(-1);
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
	this->d->lodEmergencyProgressiveCeiling = -1;
	if (this->d->viewAttachment->getViewLodState())
	    this->d->viewAttachment->getViewLodState()->
		setCadPresentationProgressiveLodCeiling(-1);
    }
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
    return this->d->lodRefinementAwaitingFrame;
}

size_t
BObolViewController::processPendingLodResults(size_t maxResults,
	uint64_t maxMicroseconds)
{
    if (!this->d->lodService)
	return 0;

    const auto clear_lod_wakeup_if_idle = [this]() {
	if (!this->d->progressiveProviders.empty() ||
	    this->d->lodSubmissionPending ||
	    this->d->lodResultsPending.load() != 0 ||
	    (this->d->lodService &&
	     this->d->lodService->queuedResultCountForDiagnostics() > 0))
	    return;

	/* A result-ready callback may race this drain.  Clear first, then
	 * recheck so a concurrent callback cannot lose its frame wakeup. */
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

    /* Apply drained results in batches rather than one-at-a-time.  Each
     * applyLodResults() call performs a single SoBRLLodUpdateAction traversal
     * of the whole scene, so draining one result per traversal makes the total
     * cost O(results * scene) and stalls the UI while a large assembly streams
     * in.  Batching collapses many results into one traversal; the microsecond
     * budget is still honored by checking the clock between batches. */
    /* The default result queue is 256.  One update action routes all of those
     * results by compact occurrence key in a single scene traversal.  The old
     * 64-result quantum traversed a many-thousand-leaf scene four times per
     * warm-cache wave and made cached box replacement visibly crawl. */
    static const size_t apply_batch = 256;
    const int64_t start = bu_gettime();
    size_t processed = 0;
    unsigned int matched = 0;
    unsigned int applied = 0;
    unsigned int rejected = 0;
    unsigned int unmatched = 0;
    SbString diagnostics;
    while (maxResults == 0 || processed < maxResults) {
	size_t batch = apply_batch;
	if (maxResults != 0) {
	    const size_t remaining = maxResults - processed;
	    if (remaining < batch)
		batch = remaining;
	}
	(void)this->applyLodResults(this->d->lodService, batch);
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
	const int64_t elapsed = bu_gettime() - start;
	if (elapsed >= 0 &&
	    static_cast<uint64_t>(elapsed) >= maxMicroseconds)
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

    SbString signature = controller_lod_source_signature(this);
    if (signature.getLength() == 0)
	return 0;

    if (this->d->lodLastSubmittedViewRevision == this->d->lodViewRevision &&
	this->d->lodLastSubmittedPolicyRevision == this->d->lodPolicyRevision &&
	bu_strcmp(this->d->lodLastSubmittedSourceSignature.getString(),
	       signature.getString()) == 0) {
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
	this->d->lodLastSubmittedSourceSignature.getLength() > 0 &&
	bu_strcmp(this->d->lodLastSubmittedSourceSignature.getString(),
	    signature.getString()) != 0;
    if (sourceSetChanged || !this->d->lodSubmissionPending) {
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodBudgetRescanAfterFrame = FALSE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodRefinementFaceBudgetInitialized = FALSE;
	this->d->lodRefinementFaceBudgetRemaining = SIZE_MAX;
	this->d->lodRetainedAdmissionActive = FALSE;
	this->d->lodRetainedAdmissionBudgetRemaining = SIZE_MAX;
	this->d->lodRefinementBudgetRetryMicroseconds = 0;
    } else {
	this->d->lodSubmissionRescanPending = TRUE;
    }
    this->d->lodSubmissionPending = TRUE;
    this->d->lodSubmissionRefreshMissing = refreshMissing;
    this->d->lodSubmissionReset = reset;
    int submitted = this->submitLodRequests(this->d->lodService, generation,
					    refreshMissing, reset);
    if (submitted >= 0) {
	this->d->lodActiveGeneration = generation;
	this->d->lodLastSubmittedViewRevision = this->d->lodViewRevision;
	this->d->lodLastSubmittedPolicyRevision = this->d->lodPolicyRevision;
	this->d->lodLastSubmittedSourceSignature = signature;
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
		    const long double targetNanoseconds =
			targetFps > 0.0f ?
			1000000000.0L /
			    static_cast<long double>(targetFps) : 0.0L;
		    const uint64_t observedNanoseconds =
			this->d->smoothedRenderTimeNanoseconds ?
			this->d->smoothedRenderTimeNanoseconds :
			this->d->lastRenderTimeNanoseconds;
		    if (targetNanoseconds > 0.0L &&
			static_cast<long double>(observedNanoseconds) >=
			    targetNanoseconds * 0.50L) {
			growthMultipleNumerator = 5;
			growthMultipleDenominator = 4;
		    } else if (targetNanoseconds > 0.0L &&
			static_cast<long double>(observedNanoseconds) >=
			    targetNanoseconds * 0.25L) {
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
	    /* A stable cut is a retained equilibrium, not a frame-by-frame
	     * controller output.  Ignore calibration noise inside a 15% band;
	     * otherwise multi-leaf scenes repeatedly cross minimum-prefix
	     * admission thresholds and visibly flash mesh -> box -> mesh. */
	    if (!this->d->lodInteractive &&
		this->d->lodCurrentFaceBudget > this->d->lodSeedFaceBudget &&
		this->d->lodCurrentFaceBudget != SIZE_MAX) {
		const long double relative =
		    static_cast<long double>(faceBudget) /
		    static_cast<long double>(
			this->d->lodCurrentFaceBudget);
		if (relative >= 0.85L && relative <= 1.15L)
		    faceBudget = this->d->lodCurrentFaceBudget;
	    }
	    /* Refinement within one quiet camera epoch is monotonic.  A later
	     * slow frame may lower the learned capacity used by the next
	     * interaction, but must not make geometry which was already shown
	     * disappear from an otherwise unchanged view.  Camera interaction
	     * remains free to select a coarser responsive cut. */
	    if (!this->d->lodInteractive &&
		this->d->lodCurrentFaceBudget != SIZE_MAX)
		faceBudget = std::max(faceBudget,
		    this->d->lodCurrentFaceBudget);
	}
	this->d->lodCurrentFaceBudget = faceBudget;
	const size_t additionalFaces =
	    faceBudget == SIZE_MAX ? SIZE_MAX :
	    (faceBudget > activeFaces ? faceBudget - activeFaces : 0);
	this->d->lodRefinementFaceBudgetRemaining = additionalFaces;
	this->d->lodRetainedAdmissionActive =
	    faceBudget != SIZE_MAX && activeFaces > faceBudget ?
	    TRUE : FALSE;
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
	action.setRevisions(this->d->lodViewRevision, this->d->lodPolicyRevision);
	action.setRefreshMissing(refreshMissing);
	action.setReset(reset);
	action.setViewLodState(this->d->viewAttachment->getViewLodState());
	/* Coarsening a retained prefix only changes the draw count/snap level;
	 * it does not rebuild or reread the mesh.  Permit it during motion so a
	 * previously settled, expensive cut cannot pin interactive FPS. */
	action.setAllowLevelDowngrade(TRUE);
	action.setAllowRetainedRefinement(
	    this->d->lodRefinementAwaitingFrame ? FALSE : TRUE);
	action.setRefinementFaceBudget(
	    this->d->lodRefinementFaceBudgetRemaining);
	if (this->d->lodRetainedAdmissionActive)
	action.setRetainedSceneFaceBudget(
		this->d->lodRetainedAdmissionBudgetRemaining);
	action.setSubmissionTaskLimit(capacity);
	/* A complete projected priority sort is valuable for a quiet view, but
	 * it is the wrong unit of work for an input turn containing tens of
	 * thousands of leaves.  During interaction, start with deterministic
	 * compact-index order and consume it in bounded waves.  Every occurrence
	 * is still projected against the newest camera and charged to the same
	 * aggregate face budget; only global ordering is deferred until the
	 * stable pass.  This keeps warm-cache activation and emergency
	 * coarsening progressive without a 50-200 ms monolithic UI stall. */
	static const size_t interactiveCompactWave = 2048;
	const int compactCount = source->getCompactInstanceCount();
	const SbBool boundedInteractiveCompact =
	    this->d->lodInteractive &&
	    compactCount > static_cast<int>(interactiveCompactWave) ?
	    TRUE : FALSE;
	if (boundedInteractiveCompact &&
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
	}
	if (this->d->lodSubmissionPlanValid &&
	    this->d->lodSubmissionPlanSource == source)
	    action.setCompactEntryPlan(
		this->d->lodSubmissionPlanEntries);
	action.setCompactEntryRange(this->d->lodSubmissionEntryOffset,
	    boundedInteractiveCompact ? interactiveCompactWave : SIZE_MAX);
	if (this->d->lodUseForcedLevel)
	    action.setForcedLevel(this->d->lodForcedLevel);
	action.apply(source);
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
	if (action.getUpdatedCutCount() > 0)
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
	if (this->d->lodRetainedAdmissionActive) {
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
    if (completedPass && this->d->lodRetainedRefinementPending &&
	this->d->lodRetainedRefinementCutAdvanced) {
	/* The just-completed pass selected the next resident cut using the
	 * newest view already.  Present it before any requested rescan; the
	 * post-frame submission is itself a full current-view pass. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = FALSE;
	this->d->lodPassAdmittedWork = FALSE;
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
	const uint64_t observedStableNanoseconds =
	    this->d->lastRenderTimeNanoseconds ?
	    this->d->lastRenderTimeNanoseconds :
	    this->d->smoothedRenderTimeNanoseconds;
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
	const SbBool calibrationProbe =
	    !this->d->lodInteractive &&
	    !this->d->lodPassAdmittedWork &&
	    sceneActiveFaces > 0 &&
	    this->d->lodCurrentFaceBudget != SIZE_MAX &&
	    stableTargetNanoseconds > 0 &&
	    (observedStableNanoseconds <
		static_cast<uint64_t>(
		    static_cast<long double>(stableTargetNanoseconds) * 0.90L) ||
	     calibratedBudgetHeadroom);
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = FALSE;
	this->d->lodBudgetRescanAfterFrame =
	    this->d->lodPassAdmittedWork || calibrationProbe;
	if (calibrationProbe)
	    this->requestRender("lod-budget-calibration");
	append_controller_lod_diagnostic(this->d->lastLodDiagnostics,
	    SbString(""),
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
    } else if (completedPass && this->d->lodSubmissionRescanPending) {
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionPending = sources.empty() ? FALSE : TRUE;
    } else {
	this->d->lodSubmissionPending = completedPass ? FALSE : TRUE;
	if (completedPass)
	    this->d->lodPassAdmittedWork = FALSE;
    }
    if (!this->d->lodSubmissionPending &&
	this->d->lodRetainedRefinementPending &&
	this->d->lodRetainedRefinementCutAdvanced) {
	/* A richer prefix is already in memory, but expose it one cut per
	 * completed frame.  Retargeting is metadata/draw-count only; no
	 * provider task, cache read, or geometry rebuild is involved. */
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

    if (this->d->lastLodUpdatedCutCount > 0)
	this->requestRender("lod-cut");

    return size_to_int_saturated(
	       static_cast<size_t>(this->d->lastLodSubmittedTaskCount));
}

int
BObolViewController::applyLodResults(BObolLodService *service,
				       size_t maxResults)
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

    SoNode *root = this->getRenderSceneRoot();
    if (!root)
	root = this->getSceneRoot();
    if (!root) {
	this->d->lastLodDiagnostics = "LoD result application requires a scene root";
	return -1;
    }

    std::vector<BObolLodResult> drained;
    this->d->lastLodResultCount = service->drainResults(drained, maxResults);
    this->d->lodResultsPending.store(
	service->queuedResultCountForDiagnostics() > 0 ? 1 : 0);
    if (this->d->lastLodResultCount == 0)
	return 0;

    SoBRLLodUpdateAction update;
    update.setViewLodState(this->d->viewAttachment->getViewLodState());
    SbBool partialRefinementCandidate = FALSE;
    for (size_t i = 0; i < drained.size(); i++) {
	const SbBool traceResult = controller_lod_trace_result(drained[i]);
	if (traceResult) {
	    const BObolViewLodState::CadPayload *resident =
		this->d->viewAttachment->getViewLodState()->
		    findCadForResult(drained[i]);
	    bu_log("BObol LoD apply trace object=%s request_level=%d "
		   "loaded_level=%d incoming_view=%llu incoming_policy=%llu "
		   "current_view=%llu current_policy=%llu resident_level=%d "
		   "resident_view=%llu resident_policy=%llu\n",
		   drained[i].request.objectName.getString(),
		   drained[i].request.requestedLevel,
		   drained[i].geometry.activeLevel,
		   static_cast<unsigned long long>(
		       drained[i].request.viewRevision),
		   static_cast<unsigned long long>(
		       drained[i].request.policyRevision),
		   static_cast<unsigned long long>(this->d->lodViewRevision),
		   static_cast<unsigned long long>(this->d->lodPolicyRevision),
		   resident ? resident->activeLevel : -1,
		   static_cast<unsigned long long>(
		       resident ? resident->viewRevision : 0),
		   static_cast<unsigned long long>(
		       resident ? resident->policyRevision : 0));
	}
	if (drained[i].request.viewRevision != this->d->lodViewRevision ||
	    drained[i].request.policyRevision != this->d->lodPolicyRevision) {
	    /* A completed mesh from the prior camera remains a valid progressive
	     * bootstrap when this occurrence still has only its box or an older
	     * mesh epoch.  Reject it once equal/newer view data is resident so an
	     * out-of-order completion cannot replay a coarse cut. */
	    const BObolViewLodState::CadPayload *cadPayload =
		this->d->viewAttachment->getViewLodState()->
		    findCadForResult(drained[i]);
	    const BObolViewLodState::MeshPayload *meshPayload =
		this->d->viewAttachment->getViewLodState()->
		    findMeshForResult(drained[i]);
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
		this->d->lastLodRejectedResultCount++;
		append_controller_lod_diagnostic(this->d->lastLodDiagnostics,
						 drained[i].request.objectPath,
						 "superseded LoD result epoch rejected");
		continue;
	    }
	}
	if (drained[i].providerStatus == BOBOL_LOD_PROVIDER_READY &&
	    drained[i].resultKind == BOBOL_LOD_RESULT_MESH &&
	    drained[i].geometry.activeLevel >= 0 &&
	    drained[i].request.requestedLevel >
		drained[i].geometry.activeLevel &&
	    !drained[i].terminal)
	    partialRefinementCandidate = TRUE;
	update.addResult(std::move(drained[i]));
    }

    update.apply(root);
    this->d->lastLodMatchedResultCount = update.getMatchedResultCount();
    this->d->lastLodAppliedResultCount = update.getAppliedResultCount();
    this->d->lastLodRejectedResultCount += update.getRejectedResultCount();
    this->d->lastLodUnmatchedResultCount = update.getUnmatchedResultCount();
    if (update.getDiagnostics().getLength() > 0) {
	if (this->d->lastLodDiagnostics.getLength() > 0)
	    this->d->lastLodDiagnostics += "\n";
	this->d->lastLodDiagnostics += update.getDiagnostics();
    }
    this->d->lodResultsPending.store(
	service->queuedResultCountForDiagnostics() > 0 ? 1 : 0);

    if (this->d->lastLodAppliedResultCount > 0) {
	this->requestRender("lod-result");
	(void)this->enforceMeshResidencyBudget();
	if (partialRefinementCandidate)
	    this->scheduleLodRefinementFrame("lod-result");
    }

    return size_to_int_saturated(
	       static_cast<size_t>(this->d->lastLodAppliedResultCount));
}

uint64_t
BObolViewController::getLodViewRevision(void) const
{
    return this->d->lodViewRevision;
}

void
BObolViewController::setLodPolicyRevision(uint64_t revision)
{
    uint64_t newRevision = revision == 0 ? 1 : revision;
    if (this->d->lodPolicyRevision == newRevision)
	return;
    this->d->lodPolicyRevision = newRevision;
    this->requestRender("lod-policy");
}

uint64_t
BObolViewController::getLodPolicyRevision(void) const
{
    return this->d->lodPolicyRevision;
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
    this->d->lodViewRevision++;
    if (this->d->lodViewRevision == 0)
	this->d->lodViewRevision++;
}

void
BObolViewController::advanceLodPolicyRevision(void)
{
    this->d->lodPolicyRevision++;
    if (this->d->lodPolicyRevision == 0)
	this->d->lodPolicyRevision++;
}

static float
controller_interactive_pixel_error(uint64_t renderNanoseconds,
				   uint64_t presentationNanoseconds)
{
    /* Four pixels is the normal motion cut.  If recent frames miss a 60 Hz
     * budget, increase error with the square root of the overrun: PoP level
     * population commonly grows near quadratically.  A 16-pixel ceiling was
     * still a 100-200 ms motion frame for large meshes in OSMesa.  Permit up
     * to 64 pixels during interaction; the quiet-view path remains exactly
     * one pixel and progressively restores terminal visual accuracy. */
    const double targetNanoseconds = 1000000000.0 / 60.0;
    const uint64_t observed = std::max(renderNanoseconds,
	presentationNanoseconds);
    if (!observed ||
	static_cast<double>(observed) <= targetNanoseconds * 1.20)
	return 4.0f;
    const double scale = std::sqrt(
	static_cast<double>(observed) / targetNanoseconds);
    return static_cast<float>(std::max(4.0, std::min(64.0, 4.0 * scale)));
}

static int
controller_interactive_progressive_ceiling(
    BObolViewLodState *viewState, float pixelError)
{
    if (!viewState)
	return -1;
    const int activeMaximum =
	viewState->maximumActiveProgressiveLevel();
    if (activeMaximum < 0)
	return -1;
    const double error = std::max(1.0,
	static_cast<double>(pixelError));
    const int reduction = std::max(1,
	static_cast<int>(std::floor(std::log2(error))));
    return std::max(0, activeMaximum - reduction);
}

void
BObolViewController::beginLodInteraction(void)
{
    if (!this->d->lodAutoSubmit || this->d->lodGestureActive)
	return;

    const float previousPixelError = this->d->lodTargetPixelError;
    this->d->lodGestureActive = TRUE;
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
	presentationNanoseconds =
	    this->d->smoothedInteractivePresentationIntervalNanoseconds;
    }
    this->d->lodTargetPixelError =
	controller_interactive_pixel_error(
	    this->d->smoothedRenderTimeNanoseconds,
	    presentationNanoseconds);
    BObolViewLodState *viewState =
	this->d->viewAttachment->getViewLodState();
    this->d->lodEmergencyProgressiveCeiling =
	controller_interactive_progressive_ceiling(
	    viewState, this->d->lodTargetPixelError);
    if (viewState)
	viewState->setCadPresentationProgressiveLodCeiling(
	    this->d->lodEmergencyProgressiveCeiling);
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
    if (!this->d->lodGestureActive)
	return;

    this->d->lodGestureActive = FALSE;
    {
	std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
	this->d->lastInteractivePresentationTimestampNanoseconds = 0;
    }
    /* Keep the motion cut through the normal quiet-view debounce.  Refining
     * immediately on button release makes the release event itself block on
     * a full software frame and defeats coarse interaction. */
    this->d->lodLastViewChangeMicroseconds = bu_gettime();
    this->d->lodSettleAfterRenderSerial = 0;
    this->markProgressiveWorkPending();
    this->requestRender("lod-interaction-end");
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

float
BObolViewController::getLodTargetPixelError(void) const
{
    return this->d->lodTargetPixelError;
}

int
BObolViewController::getLodEmergencyProgressiveCeiling(void) const
{
    return this->d->lodEmergencyProgressiveCeiling;
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
    SbString signature = controller_lod_view_signature(view, haveCamera,
	this->d->activeCamera, this->d->viewportRegion);

    if (bu_strcmp(this->d->lodViewSignature.getString(), signature.getString()) == 0)
	return;

    this->d->lodViewSignature = signature;
    if (advanceOnChange) {
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
	    this->d->lodLastViewChangeMicroseconds = bu_gettime();
	    this->d->lodInteractive = TRUE;
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
		    presentationNanoseconds);
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    if (this->d->lodEmergencyProgressiveCeiling < 0)
		this->d->lodEmergencyProgressiveCeiling =
		    controller_interactive_progressive_ceiling(
			viewState, this->d->lodTargetPixelError);
	    if (viewState)
		viewState->setCadPresentationProgressiveLodCeiling(
		    this->d->lodEmergencyProgressiveCeiling);
	    this->markProgressiveWorkPending();
	}
    }
}
