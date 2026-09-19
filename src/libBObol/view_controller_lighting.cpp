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
#include "scalar_publication_private.h"
#include "view_controller_private.h"

#include <Inventor/SbName.h>
#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/misc/SoChildList.h>
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
#include <Inventor/tools/SbModernUtils.h>

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

static constexpr size_t controller_environment_child_count = 6;
static constexpr float controller_default_ambient_intensity = 0.18f;
static constexpr float controller_mged_ambient_intensity = 0.30f;
static constexpr float controller_default_headlight_intensity = 0.68f;
static constexpr float controller_default_fill_intensity = 0.22f;
static constexpr float controller_default_rim_intensity = 0.18f;
static constexpr float controller_minimum_light_intensity = 0.0f;
static constexpr float controller_maximum_light_intensity = 1.0f;
static constexpr float controller_light_intensity_tolerance = 1.0e-6f;
static constexpr int controller_antialiasing_pass_count = 1;
static constexpr float controller_clip_plane_tolerance =
    32.0f * std::numeric_limits<float>::epsilon();
static constexpr float controller_minimum_direction_length = 1.0e-6f;
static constexpr float controller_unit_direction_tolerance = 1.0e-6f;
static constexpr double controller_full_to_half_angle = 0.5;
static constexpr unsigned int controller_minimum_viewport_extent = 1;
static constexpr unsigned int controller_maximum_viewport_extent =
    static_cast<unsigned int>(std::numeric_limits<short>::max());

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

static bool
controller_normalize_direction(SbVec3f &direction)
{
    if (!std::isfinite(direction[0]) || !std::isfinite(direction[1]) ||
	!std::isfinite(direction[2]))
	return false;
    const double x = direction[0];
    const double y = direction[1];
    const double z = direction[2];
    const double length = std::sqrt(x * x + y * y + z * z);
    if (!std::isfinite(length) ||
	length <= static_cast<double>(controller_minimum_direction_length))
	return false;
    direction *= static_cast<float>(1.0 / length);
    return true;
}

double
controller_aspect_from_region(const SbViewportRegion &region)
{
    SbVec2s window = region.getWindowSize();
    if (window[0] <= 0 || window[1] <= 0)
	return 0.0;

    return static_cast<double>(window[0]) / static_cast<double>(window[1]);
}

SbViewportRegion
controller_viewport_region_with_size(const SbViewportRegion &region,
	unsigned int width, unsigned int height)
{
    width = std::min(std::max(width, controller_minimum_viewport_extent),
	controller_maximum_viewport_extent);
    height = std::min(std::max(height, controller_minimum_viewport_extent),
	controller_maximum_viewport_extent);
    SbViewportRegion sized = region;
    sized.setWindowSize(static_cast<short>(width), static_cast<short>(height));
    sized.setViewportPixels(0, 0, static_cast<short>(width),
	static_cast<short>(height));
    return sized;
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
controller_find_cutting_affordance(SoGroup *presentationRoot)
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
    return NULL;
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

void
controller_configure_render_environment(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return;

    SoSeparator *root = viewport->getRoot();
    SoGroup *renderEnvironment = controller_find_render_environment(root);
    SoDepthBuffer *depthBuffer = NULL;
    SoEnvironment *environment = NULL;
    SoLightModel *lightModel = NULL;
    SoClipPlane *clipMinimum = NULL;
    SoClipPlane *clipMaximum = NULL;
    SoClipPlane *cuttingPlane = NULL;
    if (renderEnvironment) {
	for (int i = 0; i < renderEnvironment->getNumChildren(); i++) {
	    SoNode *child = renderEnvironment->getChild(i);
	    if (!depthBuffer && child &&
		child->isOfType(SoDepthBuffer::getClassTypeId()))
		depthBuffer = static_cast<SoDepthBuffer *>(child);
	    if (!environment && child &&
		child->isOfType(SoEnvironment::getClassTypeId()))
		environment = static_cast<SoEnvironment *>(child);
	    if (!lightModel && child &&
		child->isOfType(SoLightModel::getClassTypeId()))
		lightModel = static_cast<SoLightModel *>(child);
	    if (!child || !child->isOfType(SoClipPlane::getClassTypeId()))
		continue;
	    const char *name = child->getName().getString();
	    if (!clipMinimum && bu_strcmp(name,
		    controller_clip_plane_name(TRUE)) == 0)
		clipMinimum = static_cast<SoClipPlane *>(child);
	    if (!clipMaximum && bu_strcmp(name,
		    controller_clip_plane_name(FALSE)) == 0)
		clipMaximum = static_cast<SoClipPlane *>(child);
	    if (!cuttingPlane && bu_strcmp(name, "BObolCuttingPlane") == 0)
		cuttingPlane = static_cast<SoClipPlane *>(child);
	}
    }

    SoDirectionalLight *headlight = controller_find_camera_light(root,
	controller_headlight_name());
    SoDirectionalLight *fill = controller_find_camera_light(root,
	controller_studio_fill_name());
    SoDirectionalLight *rim = controller_find_camera_light(root,
	controller_studio_rim_name());
    int cameraIndex = -1;
    int renderEnvironmentOccurrences = 0;
    int headlightOccurrences = 0;
    int fillOccurrences = 0;
    int rimOccurrences = 0;
    for (int i = 0; i < root->getNumChildren(); ++i) {
	SoNode *child = root->getChild(i);
	renderEnvironmentOccurrences += child == renderEnvironment;
	headlightOccurrences += child == headlight;
	fillOccurrences += child == fill;
	rimOccurrences += child == rim;
	if (cameraIndex < 0 && child &&
	    child->isOfType(SoCamera::getClassTypeId())) {
	    cameraIndex = i;
	}
    }
    const int rigAnchor = cameraIndex >= 0 ? cameraIndex : 0;
    if (renderEnvironment && depthBuffer && environment && lightModel &&
	clipMinimum && clipMaximum && cuttingPlane && headlight && fill && rim &&
	renderEnvironmentOccurrences == 1 && headlightOccurrences == 1 &&
	fillOccurrences == 1 && rimOccurrences == 1 &&
	root->findChild(renderEnvironment) == 0 &&
	root->findChild(headlight) == rigAnchor + 1 &&
	root->findChild(fill) == rigAnchor + 2 &&
	root->findChild(rim) == rigAnchor + 3)
	return;

    /* Build every missing owned node under a temporary reference.  A plugin
     * may retain this root between controller attachments, so repair is a
     * live publication and must leave its extension children intact. */
    SbModernUtils::SoNodeRef renderEnvironmentOwner(NULL);
    SbModernUtils::SoNodeRef depthBufferOwner(NULL);
    SbModernUtils::SoNodeRef environmentOwner(NULL);
    SbModernUtils::SoNodeRef lightModelOwner(NULL);
    SbModernUtils::SoNodeRef clipMinimumOwner(NULL);
    SbModernUtils::SoNodeRef clipMaximumOwner(NULL);
    SbModernUtils::SoNodeRef cuttingPlaneOwner(NULL);
    SbModernUtils::SoNodeRef headlightOwner(NULL);
    SbModernUtils::SoNodeRef fillOwner(NULL);
    SbModernUtils::SoNodeRef rimOwner(NULL);
    const bool newRenderEnvironment = renderEnvironment == NULL;
    if (!renderEnvironment) {
	renderEnvironmentOwner = SbModernUtils::SoNodeRef(new SoGroup);
	renderEnvironment = static_cast<SoGroup *>(renderEnvironmentOwner.get());
	renderEnvironment->setName(SbName(controller_render_environment_name()));
    }
    if (!depthBuffer) {
	depthBufferOwner = SbModernUtils::SoNodeRef(new SoDepthBuffer);
	depthBuffer = static_cast<SoDepthBuffer *>(depthBufferOwner.get());
	depthBuffer->test = TRUE;
	depthBuffer->write = TRUE;
    }
    if (!environment) {
	environmentOwner = SbModernUtils::SoNodeRef(new SoEnvironment);
	environment = static_cast<SoEnvironment *>(environmentOwner.get());
	environment->ambientIntensity = controller_default_ambient_intensity;
	environment->ambientColor = SbColor(1.0f, 1.0f, 1.0f);
    }
    if (!lightModel) {
	lightModelOwner = SbModernUtils::SoNodeRef(new SoLightModel);
	lightModel = static_cast<SoLightModel *>(lightModelOwner.get());
	lightModel->model = SoLightModel::PHONG;
    }
    if (!clipMinimum) {
	clipMinimumOwner = SbModernUtils::SoNodeRef(new SoClipPlane);
	clipMinimum = static_cast<SoClipPlane *>(clipMinimumOwner.get());
	clipMinimum->setName(SbName(controller_clip_plane_name(TRUE)));
	clipMinimum->on = FALSE;
    }
    if (!clipMaximum) {
	clipMaximumOwner = SbModernUtils::SoNodeRef(new SoClipPlane);
	clipMaximum = static_cast<SoClipPlane *>(clipMaximumOwner.get());
	clipMaximum->setName(SbName(controller_clip_plane_name(FALSE)));
	clipMaximum->on = FALSE;
    }
    if (!cuttingPlane) {
	cuttingPlaneOwner = SbModernUtils::SoNodeRef(new SoClipPlane);
	cuttingPlane = static_cast<SoClipPlane *>(cuttingPlaneOwner.get());
	cuttingPlane->setName(SbName("BObolCuttingPlane"));
	cuttingPlane->on = FALSE;
    }
    if (!headlight) {
	headlightOwner = SbModernUtils::SoNodeRef(new SoDirectionalLight);
	headlight = static_cast<SoDirectionalLight *>(headlightOwner.get());
	headlight->setName(SbName(controller_headlight_name()));
	headlight->color = SbColor(1.0f, 1.0f, 1.0f);
	headlight->intensity = controller_default_headlight_intensity;
	headlight->direction = bobol_headlight_default_offset();
    }
    if (!fill) {
	fillOwner = SbModernUtils::SoNodeRef(new SoDirectionalLight);
	fill = static_cast<SoDirectionalLight *>(fillOwner.get());
	fill->setName(SbName(controller_studio_fill_name()));
	fill->color = SbColor(1.0f, 1.0f, 1.0f);
	fill->intensity = controller_default_fill_intensity;
	fill->direction = bobol_studio_fill_offset();
    }
    if (!rim) {
	rimOwner = SbModernUtils::SoNodeRef(new SoDirectionalLight);
	rim = static_cast<SoDirectionalLight *>(rimOwner.get());
	rim->setName(SbName(controller_studio_rim_name()));
	rim->color = SbColor(1.0f, 1.0f, 1.0f);
	rim->intensity = controller_default_rim_intensity;
	rim->direction = bobol_studio_rim_offset();
    }

    std::vector<SoNode *> environmentOrder;
    environmentOrder.reserve(
	static_cast<size_t>(renderEnvironment->getNumChildren()) +
	controller_environment_child_count);
    if (!depthBufferOwner)
	for (int i = 0; i < renderEnvironment->getNumChildren(); ++i)
	    environmentOrder.push_back(renderEnvironment->getChild(i));
    else {
	environmentOrder.push_back(depthBuffer);
	for (int i = 0; i < renderEnvironment->getNumChildren(); ++i)
	    environmentOrder.push_back(renderEnvironment->getChild(i));
    }
    if (environmentOwner) environmentOrder.push_back(environment);
    if (lightModelOwner) environmentOrder.push_back(lightModel);
    if (clipMinimumOwner) environmentOrder.push_back(clipMinimum);
    if (clipMaximumOwner) environmentOrder.push_back(clipMaximum);
    if (cuttingPlaneOwner) environmentOrder.push_back(cuttingPlane);

    std::unique_ptr<SoChildList::Replacement> environmentReplacement;
    if (newRenderEnvironment) {
	for (SoNode *child : environmentOrder)
	    renderEnvironment->addChild(child);
    } else if (environmentOrder.size() !=
	    static_cast<size_t>(renderEnvironment->getNumChildren())) {
	environmentReplacement =
	    renderEnvironment->getChildren()->prepareReplacement(environmentOrder);
    }

    /* A fixed-function light is transformed when submitted.  Keep the root in
     * environment/camera/key/fill/rim order, or put the rig after the
     * environment while no camera is attached. */
    std::vector<SoNode *> rootOrder;
    rootOrder.reserve(static_cast<size_t>(root->getNumChildren()) + 4);
    rootOrder.push_back(renderEnvironment);
    cameraIndex = -1;
    for (int i = 0; i < root->getNumChildren(); ++i) {
	SoNode *child = root->getChild(i);
	if (child == renderEnvironment || child == headlight || child == fill ||
	    child == rim)
	    continue;
	rootOrder.push_back(child);
	if (cameraIndex < 0 && child &&
	    child->isOfType(SoCamera::getClassTypeId()))
	    cameraIndex = static_cast<int>(rootOrder.size()) - 1;
    }
    SoNode *cameraLights[controller_camera_light_count] = {headlight, fill, rim};
    rootOrder.insert(rootOrder.begin() +
	(cameraIndex >= 0 ? cameraIndex + 1 : 1), cameraLights,
	cameraLights + controller_camera_light_count);

    bool rootChanged = rootOrder.size() !=
	static_cast<size_t>(root->getNumChildren());
    for (size_t i = 0; !rootChanged && i < rootOrder.size(); ++i)
	rootChanged = root->getChild(static_cast<int>(i)) != rootOrder[i];
    std::unique_ptr<SoChildList::Replacement> rootReplacement;
    if (rootChanged)
	rootReplacement = root->getChildren()->prepareReplacement(rootOrder);

    if (environmentReplacement) environmentReplacement->commit();
    if (rootReplacement) rootReplacement->commit();

    std::exception_ptr failure;
    const auto notify = [&failure](auto &replacement) {
	if (!replacement) return;
	try { replacement->notify(); }
	catch (...) { if (!failure) failure = std::current_exception(); }
    };
    notify(environmentReplacement);
    notify(rootReplacement);
    if (failure) std::rethrow_exception(failure);
}

int
controller_camera_root_index(SoViewport *viewport, SoCamera *camera)
{
    if (!viewport || !viewport->getRoot())
	return 0;
    SoSeparator *root = viewport->getRoot();
    const int currentIndex = camera ? root->findChild(camera) : -1;
    if (currentIndex >= 0)
	return currentIndex;
    SoGroup *environment = controller_find_render_environment(root);
    return environment ? root->findChild(environment) + 1 : 0;
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

bool
controller_camera_relative_clip_planes(const double center[3],
	const double viewZ[3], double horizontalSize, double minimum,
	double maximum, SbPlane &minimumPlane, SbPlane &maximumPlane)
{
    if (!center || !viewZ || !std::isfinite(horizontalSize) ||
	horizontalSize <= 0.0 || !std::isfinite(minimum) ||
	!std::isfinite(maximum) || minimum > maximum)
	return false;
    for (size_t i = 0; i < 3; ++i)
	if (!std::isfinite(center[i]) || !std::isfinite(viewZ[i]))
	    return false;

    const double viewScale = horizontalSize * 0.5;
    const SbVec3f minimumNormal(
	static_cast<float>(viewZ[0]),
	static_cast<float>(viewZ[1]),
	static_cast<float>(viewZ[2]));
    if (!std::isfinite(viewScale) ||
	minimumNormal.sqrLength() <=
	    controller_minimum_direction_length *
	    controller_minimum_direction_length)
	return false;
    const SbVec3f maximumNormal = -minimumNormal;
    SbVec3f minimumPoint;
    SbVec3f maximumPoint;
    for (size_t i = 0; i < 3; ++i) {
	minimumPoint[int(i)] = static_cast<float>(
	    center[i] + viewZ[i] * minimum * viewScale);
	maximumPoint[int(i)] = static_cast<float>(
	    center[i] + viewZ[i] * maximum * viewScale);
	if (!std::isfinite(minimumPoint[int(i)]) ||
	    !std::isfinite(maximumPoint[int(i)]))
	    return false;
    }
    minimumPlane = SbPlane(minimumNormal, minimumPoint);
    maximumPlane = SbPlane(maximumNormal, maximumPoint);
    return std::isfinite(minimumPlane.getDistanceFromOrigin()) &&
	std::isfinite(maximumPlane.getDistanceFromOrigin());
}

static SbModernUtils::SoNodeRef
controller_build_cutting_plane_affordance(SoViewport *viewport,
	SoCamera *camera,
	const SbPlane &plane,
	SbBool enabled,
	double horizontalSize,
	double aspect)
{
    if (!viewport || !camera || !enabled || !std::isfinite(horizontalSize) ||
	horizontalSize <= 0.0 ||
	!std::isfinite(aspect) || aspect <= 0.0)
	return SbModernUtils::SoNodeRef(NULL);

    SbVec3f normal = plane.getNormal();
    if (normal.length() <= controller_minimum_direction_length)
	return SbModernUtils::SoNodeRef(NULL);
    normal.normalize();
    const SbVec3f reference = std::fabs(normal[Z]) < 0.8f ?
	SbVec3f(0.0f, 0.0f, 1.0f) : SbVec3f(0.0f, 1.0f, 0.0f);
    SbVec3f axisU = normal.cross(reference);
    if (axisU.length() <= controller_minimum_direction_length)
	return SbModernUtils::SoNodeRef(NULL);
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

    const SbVec2s viewportSize = viewport->getViewportRegion().getViewportSizePixels();
    if (viewportSize[0] <= 0 || viewportSize[1] <= 0)
	return SbModernUtils::SoNodeRef(NULL);
    std::vector<SbVec3f> screenPoints;
    screenPoints.reserve(points.size());
    const SbViewVolume viewVolume = camera->getViewVolume(
	static_cast<float>(aspect));
    for (const SbVec3f &point : points) {
	SbVec3f projected;
	viewVolume.projectToScreen(point, projected);
	if (!std::isfinite(projected[0]) || !std::isfinite(projected[1]) ||
	    !std::isfinite(projected[2]))
	    return SbModernUtils::SoNodeRef(NULL);
	screenPoints.push_back(SbVec3f(
	    projected[0] * static_cast<float>(viewportSize[0]),
	    projected[1] * static_cast<float>(viewportSize[1]), 0.0f));
    }
    SbModernUtils::SoNodeRef hudOwner(new SoHUDKit);
    auto *hud = static_cast<SoHUDKit *>(hudOwner.get());
    SbModernUtils::SoNodeRef widgetOwner(new SoSeparator);
    auto *widget = static_cast<SoSeparator *>(widgetOwner.get());
    SbModernUtils::SoNodeRef depthOwner(new SoDepthBuffer);
    auto *depth = static_cast<SoDepthBuffer *>(depthOwner.get());
    depth->test = FALSE;
    depth->write = FALSE;
    widget->addChild(depth);
    SbModernUtils::SoNodeRef lightingOwner(new SoLightModel);
    auto *lighting = static_cast<SoLightModel *>(lightingOwner.get());
    lighting->model = SoLightModel::BASE_COLOR;
    widget->addChild(lighting);
    SbModernUtils::SoNodeRef materialOwner(new SoMaterial);
    auto *material = static_cast<SoMaterial *>(materialOwner.get());
    material->diffuseColor = SbColor(1.0f, 0.48f, 0.08f);
    material->emissiveColor = SbColor(0.28f, 0.13f, 0.02f);
    material->transparency = 0.22f;
    widget->addChild(material);
    SbModernUtils::SoNodeRef styleOwner(new SoDrawStyle);
    auto *style = static_cast<SoDrawStyle *>(styleOwner.get());
    style->lineWidth = 1.0f;
    widget->addChild(style);
    SbModernUtils::SoNodeRef coordinatesOwner(new SoCoordinate3);
    auto *coordinates = static_cast<SoCoordinate3 *>(coordinatesOwner.get());
    coordinates->point.setValues(0, static_cast<int>(screenPoints.size()),
	screenPoints.data());
    widget->addChild(coordinates);
    SbModernUtils::SoNodeRef linesOwner(new SoLineSet);
    auto *lines = static_cast<SoLineSet *>(linesOwner.get());
    lines->numVertices.setValues(0, static_cast<int>(lineCounts.size()),
	lineCounts.data());
    widget->addChild(lines);
    hud->addWidget(widget);
    return hudOwner;
}

bool
controller_cutting_plane_affordance_update_needed(SoGroup *presentationRoot,
	SbBool enabled, bool cameraGeometryChanged)
{
    SoSeparator *affordance =
	controller_find_cutting_affordance(presentationRoot);
    if (!enabled)
	return affordance && affordance->getNumChildren() != 0;
    return cameraGeometryChanged || !affordance ||
	affordance->getNumChildren() != 1;
}

class BObolPreparedCuttingPlaneAffordance::Impl {
public:
    Impl(SoViewport *viewport, SoGroup *presentationRoot, SoCamera *camera,
	const SbPlane &plane, SbBool enabled, double horizontalSize,
	double aspect) :
	presentation(presentationRoot),
	hudOwner(controller_build_cutting_plane_affordance(viewport, camera,
	    plane, enabled, horizontalSize, aspect))
    {
	if (!viewport || !presentation)
	    return;

	SoSeparator *affordance =
	    controller_find_cutting_affordance(presentation);
	if (!hudOwner) {
	    if (affordance && affordance->getNumChildren())
		children = affordance->getChildren()->prepareReplacement({});
	    return;
	}

	if (affordance) {
	    children = affordance->getChildren()->prepareReplacement(
		{hudOwner.get()});
	    return;
	}

	affordanceOwner = SbModernUtils::SoNodeRef(new SoSeparator);
	affordance = static_cast<SoSeparator *>(affordanceOwner.get());
	affordance->setName(SbName(controller_cutting_affordance_name()));
	affordance->addChild(hudOwner.get());
	std::vector<SoNode *> next;
	next.reserve(static_cast<size_t>(presentation->getNumChildren()) + 1);
	for (int i = 0; i < presentation->getNumChildren(); ++i)
	    next.push_back(presentation->getChild(i));
	next.push_back(affordance);
	presentationChildren =
	    presentation->getChildren()->prepareReplacement(next);
    }

    void commit()
    {
	if (committed)
	    return;
	if (children)
	    children->commit();
	if (presentationChildren)
	    presentationChildren->commit();
	committed = true;
    }

    void notify(std::exception_ptr &failure)
    {
	if (!committed || notified)
	    return;
	notified = true;
	const auto publish = [&failure](auto &replacement) {
	    if (!replacement)
		return;
	    try { replacement->notify(); }
	    catch (...) { if (!failure) failure = std::current_exception(); }
	};
	publish(children);
	publish(presentationChildren);
    }

private:
    SoGroup *presentation = NULL;
    SbModernUtils::SoNodeRef hudOwner{NULL};
    SbModernUtils::SoNodeRef affordanceOwner{NULL};
    std::unique_ptr<SoChildList::Replacement> children;
    std::unique_ptr<SoChildList::Replacement> presentationChildren;
    bool committed = false;
    bool notified = false;
};

BObolPreparedCuttingPlaneAffordance::BObolPreparedCuttingPlaneAffordance(
    SoViewport *viewport, SoGroup *presentationRoot, SoCamera *camera,
    const SbPlane &plane, SbBool enabled, double horizontalSize,
    double aspect) :
    impl(std::make_unique<Impl>(viewport, presentationRoot, camera, plane,
	enabled, horizontalSize, aspect))
{
}

BObolPreparedCuttingPlaneAffordance::~BObolPreparedCuttingPlaneAffordance() =
    default;

void
BObolPreparedCuttingPlaneAffordance::commit()
{
    this->impl->commit();
}

void
BObolPreparedCuttingPlaneAffordance::notify(std::exception_ptr &failure)
{
    this->impl->notify(failure);
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
    BObolPreparedCuttingPlaneAffordance publication(viewport,
	presentationRoot, viewport->getCamera(), plane, enabled,
	horizontalSize, aspect);
    publication.commit();
    std::exception_ptr failure;
    publication.notify(failure);
    if (failure)
	std::rethrow_exception(failure);
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

std::array<SoDirectionalLight *, controller_camera_light_count>
controller_camera_lights(SoViewport *viewport)
{
    return {{
	controller_headlight(viewport),
	controller_studio_fill(viewport),
	controller_studio_rim(viewport)
    }};
}

std::array<SbVec3f, controller_camera_light_count>
controller_camera_light_directions(const SbRotation &orientation,
	const SbVec3f &headlightOffset)
{
    const std::array<SbVec3f, controller_camera_light_count> eyeDirections{{
	headlightOffset,
	bobol_studio_fill_offset(),
	bobol_studio_rim_offset()
    }};
    std::array<SbVec3f, controller_camera_light_count> worldDirections;
    for (size_t i = 0; i < worldDirections.size(); ++i) {
	orientation.multVec(eyeDirections[i], worldDirections[i]);
	(void)worldDirections[i].normalize();
    }
    return worldDirections;
}

namespace {

static bool
controller_color_is_finite(const SbColor &color)
{
    return std::isfinite(color[0]) && std::isfinite(color[1]) &&
	std::isfinite(color[2]);
}

static SbColor
controller_clamp_unit_color(const SbColor &color)
{
    return SbColor(
	std::max(0.0f, std::min(1.0f, color[0])),
	std::max(0.0f, std::min(1.0f, color[1])),
	std::max(0.0f, std::min(1.0f, color[2])));
}

static bool
controller_clip_plane_matches(const SbPlane &current,
	const SbPlane &desired)
{
    const float distanceScale = std::max(1.0f, std::max(
	std::fabs(current.getDistanceFromOrigin()),
	std::fabs(desired.getDistanceFromOrigin())));
    return (current.getNormal() - desired.getNormal()).length() <=
	controller_clip_plane_tolerance &&
	std::fabs(current.getDistanceFromOrigin() -
	    desired.getDistanceFromOrigin()) <=
	controller_clip_plane_tolerance * distanceScale;
}

template <typename Node, typename Update>
SbModernUtils::SoNodeRef
controller_scalar_candidate(Node &source, Update update)
{
    SbModernUtils::SoNodeRef owner(new Node);
    auto *candidate = static_cast<Node *>(owner.get());
    copy_publication_scalar_fields(*candidate, source);
    update(*candidate);
    return owner;
}

} // namespace

struct BObolViewController::PreparedScalarAppearancePublication {
    enum class Impact {
	PRESENTATION,
	RENDERER_CAPACITY
    };

    PreparedScalarAppearancePublication(BObolViewController &target,
	SoNode &source, const SoNode &candidate, const char *reason,
	Impact requestedImpact) : controller(target), fields(source, candidate),
	impact(requestedImpact)
    {
	if (this->impact == Impact::RENDERER_CAPACITY) {
	    this->controller.prepareRendererInvalidation(
		this->rendererInvalidation, reason);
	} else {
	    this->renderRequest = this->controller.prepareRenderRequest(reason,
		RenderRequestIntent::PRESENTATION);
	}
    }

    void commitFields() { this->fields.commit(); }

    void commitRequest() noexcept
    {
	if (this->impact == Impact::RENDERER_CAPACITY) {
	    this->controller.commitRendererInvalidation(
		this->rendererInvalidation, TRUE);
	} else {
	    this->controller.commitRenderRequest(this->renderRequest);
	}
    }

    void commit()
    {
	this->commitFields();
	this->commitRequest();
    }

    void notify()
    {
	this->fields.restore();
	std::exception_ptr failure;
	this->fields.notify(failure);
	try {
	    if (this->impact == Impact::RENDERER_CAPACITY)
		this->controller.notifyRendererInvalidation(
		    this->rendererInvalidation);
	    else
		this->controller.notifyRenderRequest(this->renderRequest);
	} catch (...) {
	    if (!failure)
		failure = std::current_exception();
	}
	if (failure)
	    std::rethrow_exception(failure);
    }

private:
    BObolViewController &controller;
    PreparedScalarFields fields;
    PreparedRendererInvalidation rendererInvalidation;
    BObolPreparedRenderRequest renderRequest;
    Impact impact;
};

struct BObolViewController::PreparedRendererAppearancePublication {
    PreparedRendererAppearancePublication(BObolViewController &target,
	const char *reason) : controller(target)
    {
	this->controller.prepareRendererInvalidation(
	    this->rendererInvalidation, reason);
    }

    void commit() noexcept
    {
	this->controller.commitRendererInvalidation(
	    this->rendererInvalidation, TRUE);
    }

    void notify(std::exception_ptr &failure)
    {
	try {
	    this->controller.notifyRendererInvalidation(
		this->rendererInvalidation);
	} catch (...) {
	    if (!failure)
		failure = std::current_exception();
	}
    }

    void notify()
    {
	std::exception_ptr failure;
	this->notify(failure);
	if (failure)
	    std::rethrow_exception(failure);
    }

private:
    BObolViewController &controller;
    PreparedRendererInvalidation rendererInvalidation;
};

static const char *
controller_scene_lights_group_name(void)
{
    return "BObolSceneLights";
}

/* Locate (creating if needed) the in-scene lights group, always positioned
 * after the camera rig in the viewport root so fixed-function GL transforms
 * their world-space positions/directions into eye space. Publish a changed
 * order once so paths and observers never see the group temporarily detached. */
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
    SbModernUtils::SoNodeRef groupOwner(NULL);
    if (!group) {
	group = new SoGroup;
	groupOwner = SbModernUtils::SoNodeRef(group);
	group->setName(SbName(controller_scene_lights_group_name()));
    }

    SoDirectionalLight *cameraLights[controller_camera_light_count] = {
	controller_find_camera_light(root, controller_headlight_name()),
	controller_find_camera_light(root, controller_studio_fill_name()),
	controller_find_camera_light(root, controller_studio_rim_name())
    };

    int groupIndex = -1;
    int groupOccurrences = 0;
    int anchorIndex = -1;
    for (int i = 0; i < root->getNumChildren(); ++i) {
	SoNode *child = root->getChild(i);
	if (child == group) {
	    if (groupIndex < 0)
		groupIndex = i;
	    ++groupOccurrences;
	}
	if ((child && child->isOfType(SoCamera::getClassTypeId())) ||
	    child == cameraLights[0] || child == cameraLights[1] ||
	    child == cameraLights[2])
	    anchorIndex = i;
    }
    if (group && groupOccurrences == 1 && groupIndex == anchorIndex + 1)
	return group;

    std::vector<SoNode *> ordered;
    ordered.reserve(static_cast<size_t>(root->getNumChildren()) + 1);
    size_t insertIndex = 0;
    bool foundAnchor = false;
    for (int i = 0; i < root->getNumChildren(); ++i) {
	SoNode *child = root->getChild(i);
	if (child == group)
	    continue;
	ordered.push_back(child);
	const bool camera = child &&
	    child->isOfType(SoCamera::getClassTypeId());
	const bool cameraLight = child == cameraLights[0] ||
	    child == cameraLights[1] || child == cameraLights[2];
	if (camera || cameraLight) {
	    insertIndex = ordered.size();
	    foundAnchor = true;
	}
    }
    if (!foundAnchor)
	insertIndex = ordered.size();
    ordered.insert(ordered.begin() + static_cast<std::ptrdiff_t>(insertIndex),
	group);

    bool unchanged = ordered.size() ==
	static_cast<size_t>(root->getNumChildren());
    for (size_t i = 0; unchanged && i < ordered.size(); ++i)
	unchanged = root->getChild(static_cast<int>(i)) == ordered[i];
    if (!unchanged) {
	auto replacement = root->getChildren()->prepareReplacement(ordered);
	replacement->commit();
	replacement->notify();
    }
    return group;
}

void
BObolViewController::setCamera(SoCamera *camera)
{
    if (camera == this->d->activeCamera)
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    controller_configure_render_environment(this->d->viewport);
    (void)controller_scene_lights_group(this->d->viewport);

    const int cameraIndex = controller_camera_root_index(
	this->d->viewport, this->d->activeCamera);

    BObolPreparedRenderRequest renderRequest = this->prepareRenderRequest(
	"camera", RenderRequestIntent::LOD_CAPACITY);
    std::unique_ptr<SoViewport::CameraReplacement> viewportReplacement =
	this->d->viewport->prepareCameraReplacement(camera, cameraIndex);
    SoCamera *previousCamera = this->d->activeCamera;
    if (camera)
	camera->ref();

    /* LoD synchronization performs its fallible snapshot work before it
     * installs a new signature.  Give it the candidate camera without waking
     * the endpoint, then publish every public camera copy and the root in one
     * non-throwing commit. */
    this->d->activeCamera = camera;
    try {
	this->syncLodViewSignature(TRUE, FALSE);
    } catch (...) {
	this->d->activeCamera = previousCamera;
	if (camera)
	    camera->unref();
	throw;
    }
    this->d->activeCamera = previousCamera;

    viewportReplacement->commit();
    this->d->activeCamera = camera;
    this->d->renderManager->setCamera(camera);
    this->commitRenderRequest(renderRequest);

    std::exception_ptr failure;
    try {
	viewportReplacement->notify();
    } catch (...) {
	failure = std::current_exception();
    }
    try {
	this->notifyRenderRequest(renderRequest);
    } catch (...) {
	if (!failure)
	    failure = std::current_exception();
    }
    if (previousCamera)
	previousCamera->unref();
    if (failure)
	std::rethrow_exception(failure);
}

SoCamera *
BObolViewController::getCamera(void) const
{
    return this->d->activeCamera;
}

void
BObolViewController::setViewportRegion(const SbViewportRegion &region)
{
    this->publishViewportRegion(region, "viewport");
}

const SbViewportRegion &
BObolViewController::getViewportRegion(void) const
{
    return this->d->viewportRegion;
}

void
BObolViewController::setViewportSize(unsigned int width, unsigned int height)
{
    PreparedViewportPublication publication;
    this->prepareViewportSizePublication(width, height, "viewport-size",
	publication);
    this->finishViewportRegionPublication(publication);
}

void
BObolViewController::prepareViewportSizePublication(unsigned int width,
	unsigned int height, const char *reason,
	PreparedViewportPublication &publication)
{
    const SbViewportRegion region = controller_viewport_region_with_size(
	this->d->viewportRegion, width, height);
    this->prepareViewportRegionPublication(region, reason, publication);
}

void
BObolViewController::publishViewportRegion(const SbViewportRegion &region,
	const char *reason)
{
    PreparedViewportPublication publication;
    this->prepareViewportRegionPublication(region, reason, publication);
    this->finishViewportRegionPublication(publication);
}

void
BObolViewController::prepareViewportRegionPublication(
    const SbViewportRegion &region, const char *reason,
    PreparedViewportPublication &publication)
{
    if (region == this->d->viewportRegion)
	return;
    publication.renderRequest = this->prepareRenderRequest(reason,
	RenderRequestIntent::LOD_CAPACITY);
    const SbViewportRegion previousRegion = this->d->viewportRegion;
    this->d->viewportRegion = region;
    this->d->viewport->setViewportRegion(region);
    this->d->renderManager->setViewportRegion(region);
    try {
	/* LoD synchronization completes its fallible convergence snapshot before
	 * installing the new signature and revision.  Restore the three public
	 * viewport copies if that preparation cannot complete. */
	this->syncLodViewSignature(TRUE, FALSE);
    } catch (...) {
	this->d->viewportRegion = previousRegion;
	this->d->viewport->setViewportRegion(previousRegion);
	this->d->renderManager->setViewportRegion(previousRegion);
	throw;
    }
    publication.changed = true;
}

void
BObolViewController::finishViewportRegionPublication(
    PreparedViewportPublication &publication)
{
    if (!publication.changed)
	return;
    this->commitViewportRegionPublication(publication);
    this->notifyViewportRegionPublication(publication);
}

void
BObolViewController::commitViewportRegionPublication(
    PreparedViewportPublication &publication) noexcept
{
    if (!publication.changed)
	return;
    this->commitRenderRequest(publication.renderRequest);
}

void
BObolViewController::notifyViewportRegionPublication(
    const PreparedViewportPublication &publication)
{
    if (!publication.changed)
	return;
    this->notifyRenderRequest(publication.renderRequest);
}

void
BObolViewController::setBackgroundColors(const SbColor &bottom,
	const SbColor &top)
{
    if (!controller_color_is_finite(bottom) ||
	!controller_color_is_finite(top))
	return;
    SoEnvironment *environment = controller_environment(this->d->viewport);
    if (!environment)
	throw std::logic_error("controller environment is unavailable");
    if (this->d->backgroundBottom == bottom &&
	this->d->backgroundTop == top &&
	environment->fogColor.getValue() == top)
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    SbModernUtils::SoNodeRef candidateOwner = controller_scalar_candidate(
	*environment, [&](SoEnvironment &candidate) {
	    candidate.fogColor = top;
	});
    /* Clear/gradient pixels do not change the CAD population or provide a
     * meaningful retained-geometry capacity sample. */
    PreparedScalarAppearancePublication publication(*this, *environment,
	*candidateOwner.get(), "background",
	PreparedScalarAppearancePublication::Impact::PRESENTATION);

    publication.commitFields();
    this->d->backgroundBottom = bottom;
    this->d->backgroundTop = top;
    publication.commitRequest();
    publication.notify();
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
    enabled = enabled ? TRUE : FALSE;
    SoDepthBuffer *depth = controller_depth_buffer(this->d->viewport);
    if (!depth || (depth->test.getValue() == enabled &&
	depth->write.getValue() == enabled))
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    SbModernUtils::SoNodeRef candidateOwner = controller_scalar_candidate(
	*depth, [=](SoDepthBuffer &candidate) {
	    candidate.test = enabled;
	    candidate.write = enabled;
	});
    PreparedScalarAppearancePublication publication(*this, *depth,
	*candidateOwner.get(), "depth-test",
	PreparedScalarAppearancePublication::Impact::RENDERER_CAPACITY);
    publication.commit();
    publication.notify();
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
    SoGLRenderAction *onscreen =
	this->d->renderManager->getGLRenderAction();
    SoGLRenderAction *offscreen = this->d->imageRenderer ?
	this->d->imageRenderer->getGLRenderAction() : NULL;
    if (!onscreen)
	throw std::logic_error("controller render action is unavailable");
    const SoGLRenderAction::TransparencyType type = enabled ?
	SoGLRenderAction::BLEND : SoGLRenderAction::NONE;
    if (this->d->transparencyEnabled == enabled &&
	onscreen->getTransparencyType() == type &&
	(!offscreen || offscreen->getTransparencyType() == type))
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    PreparedRendererAppearancePublication publication(*this, "transparency");
    this->d->transparencyEnabled = enabled;
    onscreen->setTransparencyType(type);
    if (offscreen)
	offscreen->setTransparencyType(type);
    publication.commit();
    publication.notify();
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
    SoGLRenderAction *onscreen =
	this->d->renderManager->getGLRenderAction();
    SoGLRenderAction *offscreen = this->d->imageRenderer ?
	this->d->imageRenderer->getGLRenderAction() : NULL;
    if (!onscreen)
	throw std::logic_error("controller render action is unavailable");
    const auto matches = [=](const SoGLRenderAction *action) {
	return !action || (action->isSmoothing() == enabled &&
	    action->getNumPasses() == controller_antialiasing_pass_count);
    };
    if (this->d->antialiasingEnabled == enabled && matches(onscreen) &&
	matches(offscreen))
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    PreparedRendererAppearancePublication publication(*this, "antialiasing");
    this->d->antialiasingEnabled = enabled;
    for (SoGLRenderAction *action : {onscreen, offscreen}) {
	if (!action)
	    continue;
	action->setSmoothing(enabled);
	action->setNumPasses(controller_antialiasing_pass_count);
    }
    publication.commit();
    publication.notify();
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
    SoCamera *camera = this->d->activeCamera;
    SoClipPlane *minimumNode = controller_clip_plane(this->d->viewport, TRUE);
    SoClipPlane *maximumNode = controller_clip_plane(this->d->viewport, FALSE);
    if (!camera || !minimumNode || !maximumNode)
	throw std::logic_error("controller clipping state is unavailable");

    SbVec3f viewZ;
    camera->orientation.getValue().multVec(SbVec3f(0.0f, 0.0f, 1.0f),
	viewZ);
    const SbVec3f position = camera->position.getValue();
    const double focalDistance = camera->focalDistance.getValue();
    const double center[3] = {
	static_cast<double>(position[0]) -
	    static_cast<double>(viewZ[0]) * focalDistance,
	static_cast<double>(position[1]) -
	    static_cast<double>(viewZ[1]) * focalDistance,
	static_cast<double>(position[2]) -
	    static_cast<double>(viewZ[2]) * focalDistance
    };
    const double direction[3] = {
	static_cast<double>(viewZ[0]),
	static_cast<double>(viewZ[1]),
	static_cast<double>(viewZ[2])
    };
    SbPlane minimumPlane;
    SbPlane maximumPlane;
    if (!controller_camera_relative_clip_planes(center, direction,
	    this->d->cuttingPlaneAffordanceHorizontalSize, minimum, maximum,
	    minimumPlane, maximumPlane))
	return FALSE;

    const double minimumTolerance = std::numeric_limits<double>::epsilon() *
	std::max(1.0, std::max(std::fabs(this->d->clipMinimum),
		std::fabs(minimum)));
    const double maximumTolerance = std::numeric_limits<double>::epsilon() *
	std::max(1.0, std::max(std::fabs(this->d->clipMaximum),
		std::fabs(maximum)));
    if (std::fabs(this->d->clipMinimum - minimum) <= minimumTolerance &&
	std::fabs(this->d->clipMaximum - maximum) <= maximumTolerance &&
	controller_clip_plane_matches(minimumNode->plane.getValue(),
	    minimumPlane) &&
	controller_clip_plane_matches(maximumNode->plane.getValue(),
	    maximumPlane))
	return TRUE;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    const std::array<SoClipPlane *, 2> nodes{{minimumNode, maximumNode}};
    const std::array<SbPlane, 2> planes{{minimumPlane, maximumPlane}};
    std::array<SbModernUtils::SoNodeRef, 2> candidateOwners{{
	SbModernUtils::SoNodeRef(NULL), SbModernUtils::SoNodeRef(NULL)}};
    PreparedScalarNodes scalarFields;
    scalarFields.reserve(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
	if (controller_clip_plane_matches(nodes[i]->plane.getValue(),
		planes[i]))
	    continue;
	candidateOwners[i] = controller_scalar_candidate(
	    *nodes[i], [&](SoClipPlane &candidate) {
		candidate.plane = planes[i];
	    });
	scalarFields.prepare(*nodes[i], *candidateOwners[i].get());
    }
    PreparedRendererAppearancePublication publication(*this, "clip-bounds");

    scalarFields.commit();
    this->d->clipMinimum = minimum;
    this->d->clipMaximum = maximum;
    publication.commit();

    scalarFields.restore();
    std::exception_ptr failure;
    scalarFields.notify(failure);
    publication.notify(failure);
    if (failure)
	std::rethrow_exception(failure);
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
    this->publishCuttingPlaneState(this->d->cuttingPlane, enabled);
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
    if (!std::isfinite(normal[0]) || !std::isfinite(normal[1]) ||
	!std::isfinite(normal[2]) ||
	normal.length() <= controller_minimum_direction_length ||
	!std::isfinite(distance))
	return FALSE;
    if (this->d->cuttingPlane == plane)
	return TRUE;
    this->publishCuttingPlaneState(plane, this->d->cuttingPlaneEnabled);
    return TRUE;
}

void
BObolViewController::publishCuttingPlaneState(const SbPlane &plane,
	SbBool enabled)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    controller_configure_render_environment(this->d->viewport);
    SoClipPlane *node = controller_cutting_plane(this->d->viewport);
    if (!node)
	throw std::logic_error("controller cutting plane is unavailable");

    SbModernUtils::SoNodeRef candidateOwner(new SoClipPlane);
    auto *candidate = static_cast<SoClipPlane *>(candidateOwner.get());
    copy_publication_scalar_fields(*candidate, *node);
    candidate->plane = plane;
    candidate->on = enabled;
    PreparedScalarFields fields(*node, *candidate);

    const bool geometryChanged = this->d->cuttingPlane != plane ||
	this->d->cuttingPlaneEnabled != enabled;
    std::unique_ptr<BObolPreparedCuttingPlaneAffordance> affordance;
    if (controller_cutting_plane_affordance_update_needed(
	    this->d->framebufferOverlayRoot, enabled, geometryChanged)) {
	affordance = std::make_unique<BObolPreparedCuttingPlaneAffordance>(
	    this->d->viewport, this->d->framebufferOverlayRoot,
	    this->d->activeCamera, plane, enabled,
	    this->d->cuttingPlaneAffordanceHorizontalSize,
	    this->d->cuttingPlaneAffordanceAspect);
    }
    BObolPreparedRenderRequest renderRequest = this->prepareRenderRequest(
	"cutting-plane", RenderRequestIntent::LOD_CAPACITY);

    fields.commit();
    if (affordance)
	affordance->commit();
    this->d->cuttingPlane = plane;
    this->d->cuttingPlaneEnabled = enabled;
    this->commitRenderRequest(renderRequest);

    fields.restore();
    std::exception_ptr failure;
    if (affordance)
	affordance->notify(failure);
    fields.notify(failure);
    try { this->notifyRenderRequest(renderRequest); }
    catch (...) { if (!failure) failure = std::current_exception(); }
    if (failure)
	std::rethrow_exception(failure);
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
    const auto lights = controller_camera_lights(this->d->viewport);
    if (std::any_of(lights.begin(), lights.end(),
	[](SoDirectionalLight *light) { return light == NULL; }))
	return;
    const auto worldDirections = controller_camera_light_directions(
	this->d->lastCameraOrientation, this->d->headlightOffsetEye);
    for (size_t i = 0; i < controller_camera_light_count; i++) {
	if (worldDirections[i].length() > 0.0f &&
	    lights[i]->direction.getValue() != worldDirections[i])
	    lights[i]->direction = worldDirections[i];
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
    SbVec3f offset = profile == LIGHTING_STUDIO ?
	bobol_headlight_default_offset() : bobol_mged_headlight_offset();
    (void)offset.normalize();
    this->publishLightingState(profile, offset,
	this->d->headlightCameraTracked, this->d->headlightEnabled,
	this->d->sceneLights, this->d->sceneLightsEnabled,
	"lighting-profile");
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
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    BObolViewLodState *viewState = this->getViewLodState();
    if (!viewState)
	return;
    const BObolViewLodState::NormalStyle beforeStyle =
	viewState->getNormalStyle();
    const float beforeAngle = viewState->getNormalCreaseAngle();
    viewState->setNormalStyle(style, creaseAngleDegrees);
    if (viewState->getNormalStyle() != beforeStyle ||
	std::fabs(viewState->getNormalCreaseAngle() - beforeAngle) > 1.0e-6f) {
	/* Spatial-page normals are prepared by the LoD workers.  Publish new
	 * presentation work without changing the geometry-policy or renderer-
	 * capacity epochs: normal selection preserves the submitted topology,
	 * resident payloads, and active cuts. */
	if (this->automaticLodControlEnabled()) {
	    this->d->requireExactPresentationFrame();
	    this->beginCadPresentationRepairPass();
	    this->requestLodPresentationRender("normal-style");
	} else {
	    this->requestPresentationRender("normal-style");
	}
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
    this->publishLightingState(this->d->lightingProfile,
	this->d->headlightOffsetEye, this->d->headlightCameraTracked,
	enabled ? TRUE : FALSE, this->d->sceneLights,
	this->d->sceneLightsEnabled, "lighting");
}

SbBool
BObolViewController::isHeadlightEnabled(void) const
{
    return this->d->headlightEnabled;
}

void
BObolViewController::setHeadlightCameraTracked(SbBool tracked)
{
    this->publishLightingState(this->d->lightingProfile,
	this->d->headlightOffsetEye, tracked ? TRUE : FALSE,
	this->d->headlightEnabled, this->d->sceneLights,
	this->d->sceneLightsEnabled, "lighting");
}

SbBool
BObolViewController::isHeadlightCameraTracked(void) const
{
    return this->d->headlightCameraTracked;
}

void
BObolViewController::setHeadlightOffset(const SbVec3f &eyeDir)
{
    this->publishLightingState(this->d->lightingProfile, eyeDir,
	this->d->headlightCameraTracked, this->d->headlightEnabled,
	this->d->sceneLights, this->d->sceneLightsEnabled, "lighting");
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
    SoDirectionalLight *nodes[controller_camera_light_count] = {
	controller_headlight(this->d->viewport),
	controller_studio_fill(this->d->viewport),
	controller_studio_rim(this->d->viewport)
    };
    const char *names[controller_camera_light_count] = {
	"camera-key", "camera-fill", "camera-rim"};
    for (size_t i = 0; i < controller_camera_light_count; i++) {
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
    this->publishLightingState(this->d->lightingProfile,
	this->d->headlightOffsetEye, this->d->headlightCameraTracked,
	this->d->headlightEnabled, this->d->sceneLights,
	enabled ? TRUE : FALSE, "lighting");
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

static bool
controller_scene_light_float_equal(float left, float right)
{
    return std::memcmp(&left, &right, sizeof(float)) == 0;
}

static bool
controller_scene_light_equal(const BObolSceneLightRealization &left,
	const BObolSceneLightRealization &right)
{
    return left.kind == right.kind && left.position == right.position &&
	left.direction == right.direction && left.color == right.color &&
	controller_scene_light_float_equal(left.intensity, right.intensity) &&
	controller_scene_light_float_equal(left.coneAngleDeg,
	    right.coneAngleDeg) && left.name == right.name;
}

static bool
controller_scene_lights_equal(
	const std::vector<BObolSceneLightRealization> &left,
	const std::vector<BObolSceneLightRealization> &right)
{
    return left.size() == right.size() &&
	std::equal(left.begin(), left.end(), right.begin(),
	    controller_scene_light_equal);
}

namespace {

class ControllerPreparedSceneLightChildren {
public:
    ControllerPreparedSceneLightChildren(SoGroup *group,
	const std::vector<BObolSceneLightRealization> &lights, SbBool enabled)
    {
	if (!group && !lights.empty())
	    throw std::logic_error("controller scene-light group is unavailable");
	if (!group)
	    return;

	this->owners.reserve(lights.size());
	std::vector<SoNode *> publishedLights;
	publishedLights.reserve(lights.size());
	for (const BObolSceneLightRealization &light : lights) {
	    SbModernUtils::SoNodeRef nodeOwner(NULL);
	    SoLight *node = NULL;
	    if (light.kind == BOBOL_SCENE_LIGHT_DIRECTIONAL) {
		nodeOwner = SbModernUtils::SoNodeRef(new SoDirectionalLight);
		auto *directional = static_cast<SoDirectionalLight *>(
		    nodeOwner.get());
		directional->direction = light.direction;
		node = directional;
	    } else if (light.kind == BOBOL_SCENE_LIGHT_SPOT) {
		nodeOwner = SbModernUtils::SoNodeRef(new SoSpotLight);
		auto *spot = static_cast<SoSpotLight *>(nodeOwner.get());
		spot->location = light.position;
		spot->direction = light.direction;
		/* Database angle is full beam dispersion; Coin stores the
		 * half-angle from the axis and accepts at most 90 degrees. */
		const float maximumCutoff = static_cast<float>(M_PI_2);
		spot->cutOffAngle = std::min(maximumCutoff,
		    static_cast<float>(light.coneAngleDeg *
			controller_full_to_half_angle * DEG2RAD));
		node = spot;
	    } else {
		nodeOwner = SbModernUtils::SoNodeRef(new SoPointLight);
		auto *point = static_cast<SoPointLight *>(nodeOwner.get());
		point->location = light.position;
		node = point;
	    }
	    node->color = light.color;
	    node->intensity = light.intensity;
	    node->on = enabled;
	    publishedLights.push_back(node);
	    this->owners.push_back(std::move(nodeOwner));
	}
	this->replacement = group->getChildren()->prepareReplacement(
	    publishedLights);
    }

    void commit()
    {
	if (this->replacement)
	    this->replacement->commit();
    }

    void notify(std::exception_ptr &failure)
    {
	if (!this->replacement)
	    return;
	try { this->replacement->notify(); }
	catch (...) { if (!failure) failure = std::current_exception(); }
    }

private:
    std::vector<SbModernUtils::SoNodeRef> owners;
    std::unique_ptr<SoChildList::Replacement> replacement;
};

class ControllerPreparedSceneLightEnablement {
public:
    ControllerPreparedSceneLightEnablement(SoGroup *group, SbBool enabled)
    {
	if (!group)
	    return;
	const size_t count = static_cast<size_t>(group->getNumChildren());
	this->candidates.reserve(count);
	this->fields.reserve(count);
	for (int i = 0; i < group->getNumChildren(); ++i) {
	    SoNode *child = group->getChild(i);
	    if (!child || !child->isOfType(SoLight::getClassTypeId()))
		continue;
	    SbModernUtils::SoNodeRef candidateOwner(static_cast<SoNode *>(
		child->getTypeId().createInstance()));
	    SoNode *candidate = candidateOwner.get();
	    if (!candidate || !candidate->isOfType(SoLight::getClassTypeId()))
		throw std::logic_error("controller has unsupported scene light");
	    copy_publication_scalar_fields(*candidate, *child);
	    static_cast<SoLight *>(candidate)->on = enabled;
	    this->candidates.push_back(std::move(candidateOwner));
	    this->fields.push_back(std::make_unique<PreparedScalarFields>(
		*child, *candidate));
	}
    }

    void commit()
    {
	for (const auto &field : this->fields)
	    field->commit();
    }

    void restore()
    {
	for (const auto &field : this->fields)
	    field->restore();
    }

    void notify(std::exception_ptr &failure)
    {
	this->restore();
	for (const auto &field : this->fields)
	    field->notify(failure);
    }

private:
    std::vector<SbModernUtils::SoNodeRef> candidates;
    std::vector<std::unique_ptr<PreparedScalarFields>> fields;
};

class ControllerPreparedLightingFields {
public:
    ControllerPreparedLightingFields(SoEnvironment *environment,
	const std::array<SoDirectionalLight *, controller_camera_light_count> &lights,
	BObolViewController::LightingProfile profile, bool profileChanged,
	const std::array<SbVec3f, controller_camera_light_count> &directions,
	bool aimDirections, SbBool headlightEnabled, SbBool masterLightingEnabled)
    {
	if (!environment || std::any_of(lights.begin(), lights.end(),
		[](SoDirectionalLight *light) { return light == NULL; }))
	    throw std::logic_error("controller camera-light rig is unavailable");

	this->candidates.reserve(controller_camera_light_count + 1);
	this->fields.reserve(controller_camera_light_count + 1);
	this->prepareCandidate(environment, [=](SoEnvironment &candidate) {
	    if (!profileChanged)
		return;
	    candidate.ambientColor = SbColor(1.0f, 1.0f, 1.0f);
	    candidate.ambientIntensity = profile ==
		BObolViewController::LIGHTING_STUDIO ?
		controller_default_ambient_intensity :
		controller_mged_ambient_intensity;
	});
	const SbBool studio = profile == BObolViewController::LIGHTING_STUDIO;
	const SbBool lightOn = headlightEnabled && masterLightingEnabled;
	const std::array<SbBool, controller_camera_light_count> enabled{{
	    lightOn, lightOn && studio, lightOn && studio}};
	const std::array<float, controller_camera_light_count> intensity{{
	    studio ? controller_default_headlight_intensity : 1.0f,
	    controller_default_fill_intensity,
	    controller_default_rim_intensity}};
	for (size_t i = 0; i < lights.size(); ++i) {
	    this->prepareCandidate(lights[i], [=](SoDirectionalLight &candidate) {
		if (profileChanged) {
		    candidate.color = SbColor(1.0f, 1.0f, 1.0f);
		    candidate.intensity = intensity[i];
		}
		candidate.on = enabled[i];
		if (aimDirections)
		    candidate.direction = directions[i];
	    });
	}
    }

    ControllerPreparedLightingFields(SoLightModel *model,
	const std::array<SoDirectionalLight *, controller_camera_light_count> &lights,
	BObolViewController::LightingProfile profile, SbBool headlightEnabled,
	SbBool masterLightingEnabled)
    {
	if (!model || std::any_of(lights.begin(), lights.end(),
		[](SoDirectionalLight *light) { return light == NULL; }))
	    throw std::logic_error("controller master-lighting rig is unavailable");

	this->candidates.reserve(controller_camera_light_count + 1);
	this->fields.reserve(controller_camera_light_count + 1);
	this->prepareCandidate(model, [=](SoLightModel &candidate) {
	    candidate.model = masterLightingEnabled ? SoLightModel::PHONG :
		SoLightModel::BASE_COLOR;
	});

	const SbBool lightOn = masterLightingEnabled && headlightEnabled;
	const SbBool studio = profile == BObolViewController::LIGHTING_STUDIO;
	const std::array<SbBool, controller_camera_light_count> enabled{{
	    lightOn, lightOn && studio, lightOn && studio}};
	for (size_t i = 0; i < lights.size(); ++i) {
	    this->prepareCandidate(lights[i], [=](SoDirectionalLight &candidate) {
		candidate.on = enabled[i];
	    });
	}
    }

    void commit()
    {
	for (const auto &field : this->fields)
	    field->commit();
    }

    void restore()
    {
	for (const auto &field : this->fields)
	    field->restore();
    }

    void notify(std::exception_ptr &failure)
    {
	this->restore();
	for (const auto &field : this->fields)
	    field->notify(failure);
    }

private:
    template <typename Node, typename Update>
    void prepareCandidate(Node *source, Update update)
    {
	SbModernUtils::SoNodeRef owner = controller_scalar_candidate(
	    *source, update);
	SoNode *candidate = owner.get();
	this->candidates.push_back(std::move(owner));
	this->fields.push_back(std::make_unique<PreparedScalarFields>(
	    *source, *candidate));
    }

    std::vector<SbModernUtils::SoNodeRef> candidates;
    std::vector<std::unique_ptr<PreparedScalarFields>> fields;
};

} // namespace

void
BObolViewController::setLightingEnabled(SbBool enabled)
{
    enabled = enabled ? TRUE : FALSE;
    SoLightModel *model = controller_light_model(this->d->viewport);
    const auto lights = controller_camera_lights(this->d->viewport);
    if (!model || std::any_of(lights.begin(), lights.end(),
	    [](SoDirectionalLight *light) { return light == NULL; }))
	return;

    const int requestedModel = enabled ? SoLightModel::PHONG :
	SoLightModel::BASE_COLOR;
    const SbBool lightOn = enabled && this->d->headlightEnabled;
    const SbBool studio = this->d->lightingProfile == LIGHTING_STUDIO;
    const std::array<SbBool, controller_camera_light_count> requestedLights{{
	lightOn, lightOn && studio, lightOn && studio}};
    if (model->model.getValue() == requestedModel &&
	std::equal(lights.begin(), lights.end(), requestedLights.begin(),
	    [](const SoDirectionalLight *light, SbBool requested) {
		return light->on.getValue() == requested;
	    }))
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    ControllerPreparedLightingFields lighting(model, lights,
	this->d->lightingProfile, this->d->headlightEnabled, enabled);
    PreparedRendererInvalidation rendererInvalidation;
    this->prepareRendererInvalidation(rendererInvalidation, "lighting");

    lighting.commit();
    this->commitRendererInvalidation(rendererInvalidation, TRUE);

    lighting.restore();
    std::exception_ptr failure;
    lighting.notify(failure);
    try { this->notifyRendererInvalidation(rendererInvalidation); }
    catch (...) { if (!failure) failure = std::current_exception(); }
    if (failure)
	std::rethrow_exception(failure);
}

void
BObolViewController::setSceneLights(
	const std::vector<BObolSceneLightRealization> &lights)
{
    this->publishLightingState(this->d->lightingProfile,
	this->d->headlightOffsetEye, this->d->headlightCameraTracked,
	this->d->headlightEnabled, lights, this->d->sceneLightsEnabled,
	"lighting");
}

void
BObolViewController::setSceneLights(
	const std::vector<BObolSceneLightRealization> &lights, SbBool enabled)
{
    this->publishLightingState(this->d->lightingProfile,
	this->d->headlightOffsetEye, this->d->headlightCameraTracked,
	this->d->headlightEnabled, lights, enabled ? TRUE : FALSE,
	"lighting");
}

void
BObolViewController::setLightingState(LightingProfile profile,
	const SbVec3f &headlightOffset, SbBool headlightCameraTracked,
	SbBool headlightEnabled,
	const std::vector<BObolSceneLightRealization> &sceneLights,
	SbBool sceneLightsEnabled)
{
    this->publishLightingState(profile, headlightOffset,
	headlightCameraTracked ? TRUE : FALSE,
	headlightEnabled ? TRUE : FALSE, sceneLights,
	sceneLightsEnabled ? TRUE : FALSE, "lighting");
}

void
BObolViewController::publishLightingState(LightingProfile profile,
	const SbVec3f &headlightOffset, SbBool headlightCameraTracked,
	SbBool headlightEnabled,
	const std::vector<BObolSceneLightRealization> &sceneLights,
	SbBool sceneLightsEnabled, const char *reason)
{
    if (profile != LIGHTING_STUDIO && profile != LIGHTING_MGED)
	return;
    SbVec3f offset = headlightOffset;
    if (!controller_normalize_direction(offset))
	return;
    if ((offset - this->d->headlightOffsetEye).length() <=
	controller_unit_direction_tolerance)
	offset = this->d->headlightOffsetEye;
    headlightCameraTracked = headlightCameraTracked ? TRUE : FALSE;
    headlightEnabled = headlightEnabled ? TRUE : FALSE;
    sceneLightsEnabled = sceneLightsEnabled ? TRUE : FALSE;

    const bool profileChanged = this->d->lightingProfile != profile;
    const bool offsetChanged = this->d->headlightOffsetEye != offset;
    const bool trackingChanged =
	this->d->headlightCameraTracked != headlightCameraTracked;
    const bool headlightChanged =
	this->d->headlightEnabled != headlightEnabled;
    const bool cameraStateChanged = profileChanged || offsetChanged ||
	trackingChanged || headlightChanged;
    const bool cameraFieldsChanged = profileChanged || offsetChanged ||
	headlightChanged || (trackingChanged && headlightCameraTracked);
    const bool sceneLightsChanged = !controller_scene_lights_equal(
	this->d->sceneLights, sceneLights);
    const bool sceneEnablementChanged =
	this->d->sceneLightsEnabled != sceneLightsEnabled;
    if (!cameraStateChanged && !sceneLightsChanged &&
	!sceneEnablementChanged)
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    std::vector<BObolSceneLightRealization> candidateSceneLights;
    if (sceneLightsChanged)
	candidateSceneLights = sceneLights;

    std::unique_ptr<ControllerPreparedLightingFields> cameraLighting;
    if (cameraFieldsChanged) {
	const bool aimDirections =
	    (this->d->headlightEnabled || headlightEnabled) &&
	    (profileChanged || offsetChanged ||
		(trackingChanged && headlightCameraTracked) ||
		(!this->d->headlightEnabled && headlightEnabled));
	cameraLighting = std::make_unique<ControllerPreparedLightingFields>(
	    controller_environment(this->d->viewport),
	    controller_camera_lights(this->d->viewport), profile,
	    profileChanged,
	    controller_camera_light_directions(
		this->d->lastCameraOrientation, offset),
	    aimDirections, headlightEnabled, this->isLightingEnabled());
    }

    SoGroup *sceneGroup = static_cast<SoGroup *>(this->getSceneLightsRoot());
    if (sceneLightsChanged && !sceneLights.empty())
	sceneGroup = controller_scene_lights_group(this->d->viewport);
    std::unique_ptr<ControllerPreparedSceneLightChildren> sceneChildren;
    std::unique_ptr<ControllerPreparedSceneLightEnablement> sceneFields;
    if (sceneLightsChanged)
	sceneChildren = std::make_unique<ControllerPreparedSceneLightChildren>(
	    sceneGroup, candidateSceneLights, sceneLightsEnabled);
    else if (sceneEnablementChanged)
	sceneFields = std::make_unique<ControllerPreparedSceneLightEnablement>(
	    sceneGroup, sceneLightsEnabled);
    BObolPreparedRenderRequest renderRequest = this->prepareRenderRequest(
	reason, RenderRequestIntent::LOD_CAPACITY);

    if (cameraLighting)
	cameraLighting->commit();
    if (sceneChildren)
	sceneChildren->commit();
    if (sceneFields)
	sceneFields->commit();
    if (sceneLightsChanged)
	this->d->sceneLights.swap(candidateSceneLights);
    this->d->lightingProfile = profile;
    this->d->headlightOffsetEye = offset;
    this->d->headlightCameraTracked = headlightCameraTracked;
    this->d->headlightEnabled = headlightEnabled;
    this->d->sceneLightsEnabled = sceneLightsEnabled;
    this->commitRenderRequest(renderRequest);

    if (cameraLighting)
	cameraLighting->restore();
    if (sceneFields)
	sceneFields->restore();
    std::exception_ptr failure;
    if (cameraLighting)
	cameraLighting->notify(failure);
    if (sceneChildren)
	sceneChildren->notify(failure);
    if (sceneFields)
	sceneFields->notify(failure);
    try { this->notifyRenderRequest(renderRequest); }
    catch (...) { if (!failure) failure = std::current_exception(); }
    if (failure)
	std::rethrow_exception(failure);
}

void
BObolViewController::rebuildSceneLights(void)
{
    SoGroup *group = static_cast<SoGroup *>(this->getSceneLightsRoot());
    if (this->d->sceneLights.empty()) {
	if (!group || !group->getNumChildren())
	    return;
    } else {
	group = controller_scene_lights_group(this->d->viewport);
	if (!group)
	    return;
    }

    ControllerPreparedSceneLightChildren children(group,
	this->d->sceneLights, this->d->sceneLightsEnabled);
    children.commit();
    std::exception_ptr failure;
    children.notify(failure);
    if (failure)
	std::rethrow_exception(failure);
}

void
BObolViewController::setHeadlightColor(const SbColor &color)
{
    if (!controller_color_is_finite(color))
	return;
    SoDirectionalLight *light = controller_headlight(this->d->viewport);
    if (!light)
	return;
    const SbColor clamped = controller_clamp_unit_color(color);
    if (light->color.getValue() == clamped)
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    SbModernUtils::SoNodeRef candidateOwner = controller_scalar_candidate(
	*light, [&](SoDirectionalLight &candidate) {
	    candidate.color = clamped;
	});
    PreparedScalarAppearancePublication publication(*this, *light,
	*candidateOwner.get(), "lighting",
	PreparedScalarAppearancePublication::Impact::PRESENTATION);
    publication.commit();
    publication.notify();
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
    if (!std::isfinite(intensity))
	return;
    SoDirectionalLight *light = controller_headlight(this->d->viewport);
    if (!light)
	return;
    const float clamped = std::max(controller_minimum_light_intensity,
	std::min(controller_maximum_light_intensity, intensity));
    if (std::fabs(light->intensity.getValue() - clamped) <=
	controller_light_intensity_tolerance)
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    SbModernUtils::SoNodeRef candidateOwner = controller_scalar_candidate(
	*light, [=](SoDirectionalLight &candidate) {
	    candidate.intensity = clamped;
	});
    PreparedScalarAppearancePublication publication(*this, *light,
	*candidateOwner.get(), "lighting",
	PreparedScalarAppearancePublication::Impact::PRESENTATION);
    publication.commit();
    publication.notify();
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
    enabled = enabled ? TRUE : FALSE;
    SoEnvironment *environment = controller_environment(this->d->viewport);
    if (!environment)
	return;
    const int requested = enabled ? SoEnvironment::HAZE :
	SoEnvironment::NONE;
    if (environment->fogType.getValue() == requested &&
	environment->fogColor.getValue() == this->d->backgroundTop &&
	std::fpclassify(environment->fogVisibility.getValue()) == FP_ZERO)
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    SbModernUtils::SoNodeRef candidateOwner = controller_scalar_candidate(
	*environment, [&](SoEnvironment &candidate) {
	    candidate.fogType = requested;
	    candidate.fogColor = this->d->backgroundTop;
	    /* Zero delegates visibility distance to the active camera volume. */
	    candidate.fogVisibility = 0.0f;
	});
    PreparedScalarAppearancePublication publication(*this, *environment,
	*candidateOwner.get(), "depth-cue",
	PreparedScalarAppearancePublication::Impact::RENDERER_CAPACITY);
    publication.commit();
    publication.notify();
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
    SoBRLViewLodGroup *lodRoot = this->d->renderLodRoot;
    SoBRLCadRenderBatch *batch =
	dynamic_cast<SoBRLCadRenderBatch *>(this->d->renderBatchRoot);
    if (this->d->softwareWireMode == mode &&
	(!lodRoot || lodRoot->getSoftwareWireMode() == mode) &&
	(!batch || batch->getSoftwareWireMode() == mode))
	return;

    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    PreparedRendererAppearancePublication publication(
	*this, "software-wire-mode");
    this->d->softwareWireMode = mode;
    if (lodRoot)
	lodRoot->setSoftwareWireMode(mode);
    if (batch)
	batch->setSoftwareWireMode(mode);
    publication.commit();
    publication.notify();
}

BObolViewController::SoftwareWireMode
BObolViewController::getSoftwareWireMode(void) const
{
    return this->d->softwareWireMode;
}
