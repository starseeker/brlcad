/*          V I E W _ C O N T R O L L E R _ L I G H T I N G . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file view_controller_lighting.cpp
 *
 * View environment, camera rig, scene lights, clipping, and renderer style.
 */

#include "common.h"

#include "bu/str.h"
#include "cad_assembly_private.h"
#include "view_controller_private.h"

#include <Inventor/SbName.h>
#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoClipPlane.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoDepthBuffer.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoEnvironment.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoLight.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoPointLight.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSpotLight.h>
#include <Inventor/SoOffscreenRenderer.h>

#include <algorithm>
#include <cmath>
#include <vector>

/* Camera-rig directions are expressed in eye space (camera looks down -Z,
 * +X right, +Y up).  Directional-light vectors describe photon travel, so a
 * positive X component places the source to the viewer's left.  The studio
 * rig is intentionally asymmetric: an equal ring behaves like ambient light
 * and erases the shape contrast this policy is meant to recover. */
SbVec3f
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

double
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

static const char *
controller_cutting_affordance_name(void)
{
    return "BObolCuttingPlaneAffordance";
}

static SoSeparator *
controller_cutting_affordance(SoGroup *presentationRoot)
{
    if (!presentationRoot)
	return NULL;

    const char *name = controller_cutting_affordance_name();
    for (int i = 0; i < presentationRoot->getNumChildren(); i++) {
	SoNode *child = presentationRoot->getChild(i);
	if (child && child->isOfType(SoSeparator::getClassTypeId()) &&
	    bu_strcmp(child->getName().getString(), name) == 0)
	    return static_cast<SoSeparator *>(child);
    }

    SoSeparator *affordance = new SoSeparator;
    affordance->setName(SbName(name));
    presentationRoot->addChild(affordance);
    return affordance;
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

void
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
	int hasCuttingPlane = 0;
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
		hasCuttingPlane |= bu_strcmp(name,
		    "BObolCuttingPlane") == 0;
	    }
	}
	if (!hasDepthBuffer)
	    renderEnvironment->insertChild(new SoDepthBuffer, 0);
	/* Repair missing owned environment nodes in place.  A plugin may retain
	 * the viewport root while controllers are detached and reattached, so this
	 * is a current scene-ownership invariant rather than a format migration. */
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
	if (!hasCuttingPlane) {
	    SoClipPlane *plane = new SoClipPlane;
	    plane->setName(SbName("BObolCuttingPlane"));
	    plane->on = FALSE;
	    renderEnvironment->addChild(plane);
	}
	const int index = root->findChild(renderEnvironment);
	if (index > 0) {
	    renderEnvironment->ref();
	    root->removeChild(renderEnvironment);
	    root->insertChild(renderEnvironment, 0);
	    renderEnvironment->unref();
	}
	controller_ensure_camera_rig(viewport, NULL);
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

    SoClipPlane *cuttingPlane = new SoClipPlane;
    cuttingPlane->setName(SbName("BObolCuttingPlane"));
    cuttingPlane->on = FALSE;
    renderEnvironment->addChild(cuttingPlane);

    root->insertChild(renderEnvironment, 0);
    controller_ensure_camera_rig(viewport, headlight);
}

SoClipPlane *
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

SoClipPlane *
controller_cutting_plane(SoViewport *viewport)
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
	if (child && child->isOfType(SoClipPlane::getClassTypeId()) &&
	    bu_strcmp(child->getName().getString(),
		"BObolCuttingPlane") == 0)
	    return static_cast<SoClipPlane *>(child);
    }
    return NULL;
}

void
controller_update_cutting_plane_affordance(SoViewport *viewport,
	SoGroup *presentationRoot,
	const SbPlane &plane,
	SbBool enabled,
	double horizontalSize,
	double aspect)
{
    if (!viewport || !presentationRoot)
	return;
    SoSeparator *affordance = controller_cutting_affordance(presentationRoot);
    if (!affordance)
	return;

    affordance->removeAllChildren();
    if (!enabled || !std::isfinite(horizontalSize) || horizontalSize <= 0.0 ||
	!std::isfinite(aspect) || aspect <= 0.0)
	return;

    SbVec3f normal = plane.getNormal();
    if (normal.length() <= 1.0e-6f)
	return;
    normal.normalize();
    const SbVec3f reference = std::fabs(normal[Z]) < 0.8f ?
	SbVec3f(0.0f, 0.0f, 1.0f) : SbVec3f(0.0f, 1.0f, 0.0f);
    SbVec3f axisU = normal.cross(reference);
    if (axisU.length() <= 1.0e-6f)
	return;
    axisU.normalize();
    SbVec3f axisV = normal.cross(axisU);
    axisV.normalize();

    constexpr int gridHalfSteps = 3;
    constexpr float minimumHalfWidth = 1.0e-4f;
    const float halfWidth = std::max(minimumHalfWidth,
	static_cast<float>(horizontalSize * 0.30));
    const float halfHeight = halfWidth / static_cast<float>(aspect);
    const SbVec3f center = normal * plane.getDistanceFromOrigin();
    std::vector<SbVec3f> points;
    std::vector<int32_t> lineCounts;
    points.reserve(static_cast<size_t>(4 + 2 * (gridHalfSteps * 2 + 1)) * 2);
    lineCounts.reserve(4 + 2 * (gridHalfSteps * 2 + 1));
    const auto addLine = [&points, &lineCounts](const SbVec3f &start,
	const SbVec3f &end) {
	points.push_back(start);
	points.push_back(end);
	lineCounts.push_back(2);
    };
    const SbVec3f lowerLeft = center - axisU * halfWidth - axisV * halfHeight;
    const SbVec3f lowerRight = center + axisU * halfWidth - axisV * halfHeight;
    const SbVec3f upperLeft = center - axisU * halfWidth + axisV * halfHeight;
    const SbVec3f upperRight = center + axisU * halfWidth + axisV * halfHeight;
    addLine(lowerLeft, lowerRight);
    addLine(lowerRight, upperRight);
    addLine(upperRight, upperLeft);
    addLine(upperLeft, lowerLeft);
    for (int step = -gridHalfSteps; step <= gridHalfSteps; step++) {
	const float fraction = static_cast<float>(step) /
	    static_cast<float>(gridHalfSteps);
	addLine(center + axisU * (fraction * halfWidth) - axisV * halfHeight,
		center + axisU * (fraction * halfWidth) + axisV * halfHeight);
	addLine(center - axisU * halfWidth + axisV * (fraction * halfHeight),
		center + axisU * halfWidth + axisV * (fraction * halfHeight));
    }

    SoCamera *camera = viewport->getCamera();
    const SbVec2s viewportSize = viewport->getViewportRegion().getViewportSizePixels();
    if (!camera || viewportSize[0] <= 0 || viewportSize[1] <= 0)
	return;
    std::vector<SbVec3f> screenPoints;
    screenPoints.reserve(points.size());
    const SbViewVolume viewVolume = camera->getViewVolume(
	static_cast<float>(aspect));
    for (const SbVec3f &point : points) {
	SbVec3f projected;
	viewVolume.projectToScreen(point, projected);
	if (!std::isfinite(projected[0]) || !std::isfinite(projected[1]) ||
	    !std::isfinite(projected[2]))
	    return;
	screenPoints.push_back(SbVec3f(
	    projected[0] * static_cast<float>(viewportSize[0]),
	    projected[1] * static_cast<float>(viewportSize[1]), 0.0f));
    }

    SoHUDKit *hud = new SoHUDKit;
    SoSeparator *widget = new SoSeparator;
    SoDepthBuffer *depth = new SoDepthBuffer;
    depth->test = FALSE;
    depth->write = FALSE;
    widget->addChild(depth);
    SoLightModel *lighting = new SoLightModel;
    lighting->model = SoLightModel::BASE_COLOR;
    widget->addChild(lighting);
    SoMaterial *material = new SoMaterial;
    material->diffuseColor = SbColor(1.0f, 0.48f, 0.08f);
    material->emissiveColor = SbColor(0.28f, 0.13f, 0.02f);
    material->transparency = 0.22f;
    widget->addChild(material);
    SoDrawStyle *style = new SoDrawStyle;
    style->lineWidth = 1.0f;
    widget->addChild(style);
    SoCoordinate3 *coordinates = new SoCoordinate3;
    coordinates->point.setValues(0, static_cast<int>(screenPoints.size()),
	screenPoints.data());
    widget->addChild(coordinates);
    SoLineSet *lines = new SoLineSet;
    lines->numVertices.setValues(0, static_cast<int>(lineCounts.size()),
	lineCounts.data());
    widget->addChild(lines);
    hud->addWidget(widget);
    affordance->addChild(hud);
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
    this->requestLodCapacityRender("camera");
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
    this->requestLodCapacityRender("viewport");
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
    this->requestLodCapacityRender("viewport-size");
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
    this->requestLodCapacityRender("background");
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
    this->invalidateRendererPerformanceHistory();
    this->requestLodCapacityRender("depth-test");
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
    this->invalidateRendererPerformanceHistory();
    this->requestLodCapacityRender("transparency");
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
    this->invalidateRendererPerformanceHistory();
    this->requestLodCapacityRender("antialiasing");
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
    const double minimumTolerance = std::numeric_limits<double>::epsilon() *
	std::max(1.0, std::max(std::fabs(this->d->clipMinimum),
		std::fabs(minimum)));
    const double maximumTolerance = std::numeric_limits<double>::epsilon() *
	std::max(1.0, std::max(std::fabs(this->d->clipMaximum),
		std::fabs(maximum)));
    if (std::fabs(this->d->clipMinimum - minimum) <= minimumTolerance &&
	std::fabs(this->d->clipMaximum - maximum) <= maximumTolerance)
	return TRUE;
    this->d->clipMinimum = minimum;
    this->d->clipMaximum = maximum;
    this->invalidateRendererPerformanceHistory();
    this->requestLodCapacityRender("clip-bounds");
    return TRUE;
}

void
BObolViewController::getClipBounds(double &minimum, double &maximum) const
{
    minimum = this->d->clipMinimum;
    maximum = this->d->clipMaximum;
}

void
BObolViewController::setCuttingPlaneEnabled(SbBool enabled)
{
    enabled = enabled ? TRUE : FALSE;
    if (this->d->cuttingPlaneEnabled == enabled)
	return;
    this->d->cuttingPlaneEnabled = enabled;
    SoClipPlane *plane = controller_cutting_plane(this->d->viewport);
    if (plane) {
	plane->plane = this->d->cuttingPlane;
	plane->on = enabled;
    }
	controller_update_cutting_plane_affordance(this->d->viewport,
	this->d->framebufferOverlayRoot,
	this->d->cuttingPlane, enabled,
	this->d->cuttingPlaneAffordanceHorizontalSize,
	this->d->cuttingPlaneAffordanceAspect);
    this->requestLodCapacityRender("cutting-plane");
}

SbBool
BObolViewController::isCuttingPlaneEnabled(void) const
{
    return this->d->cuttingPlaneEnabled;
}

SbBool
BObolViewController::setCuttingPlane(const SbPlane &plane)
{
    const SbVec3f normal = plane.getNormal();
    const float distance = plane.getDistanceFromOrigin();
    if (normal.length() <= 1.0e-6f || !std::isfinite(distance))
	return FALSE;
    this->d->cuttingPlane = plane;
    SoClipPlane *node = controller_cutting_plane(this->d->viewport);
    if (node) {
	node->plane = plane;
	node->on = this->d->cuttingPlaneEnabled;
    }
	controller_update_cutting_plane_affordance(this->d->viewport,
	this->d->framebufferOverlayRoot, plane,
	this->d->cuttingPlaneEnabled,
	this->d->cuttingPlaneAffordanceHorizontalSize,
	this->d->cuttingPlaneAffordanceAspect);
    this->requestLodCapacityRender("cutting-plane");
    return TRUE;
}

SbPlane
BObolViewController::getCuttingPlane(void) const
{
    return this->d->cuttingPlane;
}

size_t
BObolViewController::getActiveClipPlanes(
    SbPlane planes[CLIP_PLANE_CAPACITY]) const
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
    SoClipPlane *cutting = controller_cutting_plane(this->d->viewport);
    if (cutting && cutting->on.getValue())
	planes[count++] = cutting->plane.getValue();
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
    if (changed) {
	this->invalidateRendererPerformanceHistory();
	this->requestLodCapacityRender("lighting");
    }
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
    this->requestLodCapacityRender("lighting-profile");
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
	std::fabs(viewState->getNormalCreaseAngle() - beforeAngle) > 1.0e-6f) {
	this->invalidateRendererPerformanceHistory();
	this->requestLodCapacityRender("normal-style");
    }
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
    this->requestLodCapacityRender("lighting");
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
    this->requestLodCapacityRender("lighting");
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
    this->requestLodCapacityRender("lighting");
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
    this->requestLodCapacityRender("lighting");
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
    this->requestLodCapacityRender("lighting");
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
    this->requestLodCapacityRender("lighting");
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
    this->invalidateRendererPerformanceHistory();
    this->requestLodCapacityRender("depth-cue");
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
    this->invalidateRendererPerformanceHistory();
    this->requestLodCapacityRender("software-wire-mode");
}

BObolViewController::SoftwareWireMode
BObolViewController::getSoftwareWireMode(void) const
{
    return this->d->softwareWireMode;
}
