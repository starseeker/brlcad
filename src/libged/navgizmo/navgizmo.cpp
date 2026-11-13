/*                    N A V G I Z M O . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file libged/navgizmo/navgizmo.cpp
 *
 * Optional BRL-CAD navigation faceplate and the reference implementation for
 * a GED plugin which deliberately publishes and interacts with a custom Obol
 * node.  All scene access goes through the installed ged/plugin/obol.h API;
 * no libged drawing internals are used.
 */

#include "common.h"

#include "BObol/BNavigationGizmo.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "bu/str.h"
#include "bv.h"
#include "ged/plugin/obol.h"
#include "ged/view.h"

#include "../include/plugin.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/sensors/SoFieldSensor.h>
#include <Inventor/sensors/SoNodeSensor.h>

#include <cmath>
#include <new>

namespace {

static const char *navgizmo_feature_name = "faceplate::navigation_gizmo";
static const char *navgizmo_owner_role = "libged.navigation-gizmo";
static constexpr int navgizmo_overlay_order = 500;
static constexpr int navgizmo_drag_threshold_pixels = 2;
static constexpr int navgizmo_input_priority = 1000;

enum NavgizmoAction {
    NAVGIZMO_ACTION_PRESS = 0x4e470001u,
    NAVGIZMO_ACTION_RELEASE = 0x4e470002u,
    NAVGIZMO_ACTION_MOTION = 0x4e470003u
};

struct NavgizmoState;

static void navgizmo_camera_changed(void *userData, SoSensor *sensor);
static void navgizmo_camera_replaced(void *userData, SoSensor *sensor);

struct NavgizmoSnapshot {
    SoBRLNavigationGizmo::Style style;
    SoBRLNavigationGizmo::Part hoverPart;
    SoBRLNavigationGizmo::Part activePart;
    SoBRLNavigationGizmo::Part pressedPart;
    SbRotation cameraOrientation;
    int pressX;
    int pressY;
    int dragging;
    int moved;
    SbBool visible;

    NavgizmoSnapshot(void) :
	style(SoBRLNavigationGizmo::CUBE),
	hoverPart(SoBRLNavigationGizmo::PART_NONE),
	activePart(SoBRLNavigationGizmo::PART_NONE),
	pressedPart(SoBRLNavigationGizmo::PART_NONE),
	cameraOrientation(SbRotation::identity()), pressX(0), pressY(0),
	dragging(0), moved(0), visible(FALSE)
    {
    }
};

struct NavgizmoState {
    struct ged_view_context *viewContext;
    bobol_display_endpoint_t *endpoint;
    BObolViewController *controller;
    SoCamera *observedCamera;
    SoFieldSensor cameraSensor;
    SoNodeSensor cameraRootSensor;
    NavgizmoSnapshot snapshot;
    int suppressCameraPublication;
    int inputInstalled;

    NavgizmoState(void) :
	viewContext(NULL), endpoint(NULL), controller(NULL), observedCamera(NULL),
	cameraSensor(navgizmo_camera_changed, this),
	cameraRootSensor(navgizmo_camera_replaced, this),
	suppressCameraPublication(0), inputInstalled(0)
    {
	cameraSensor.setPriority(0);
	cameraRootSensor.setPriority(0);
    }

    ~NavgizmoState(void)
    {
	cameraRootSensor.detach();
	cameraSensor.detach();
	if (observedCamera)
	    observedCamera->unref();
    }
};

static SbBool
navgizmo_presentations_equal(const NavgizmoSnapshot &a,
    const NavgizmoSnapshot &b)
{
    if (a.style != b.style || a.hoverPart != b.hoverPart ||
	a.activePart != b.activePart || a.visible != b.visible)
	return FALSE;
    return !a.visible ||
	a.cameraOrientation.equals(b.cameraOrientation, 1.0e-6f);
}

static SbBool
navgizmo_track_camera(NavgizmoState *state)
{
    if (!state || !state->controller)
	return FALSE;
    SoCamera *camera = state->controller->getCamera();
    if (camera == state->observedCamera)
	return FALSE;
    state->cameraSensor.detach();
    if (camera)
	camera->ref();
    if (state->observedCamera)
	state->observedCamera->unref();
    state->observedCamera = camera;
    if (state->observedCamera)
	state->cameraSensor.attach(&state->observedCamera->orientation);
    return TRUE;
}

static SbRotation
navgizmo_camera_orientation(NavgizmoState *state)
{
    navgizmo_track_camera(state);
    return state && state->observedCamera ?
	state->observedCamera->orientation.getValue() : SbRotation::identity();
}

static void navgizmo_feature_result(const BObolCommandResult &result,
    void *userData);

static NavgizmoState *
navgizmo_state(BObolViewController *controller,
    BObolFeatureHandle *handleOut = NULL)
{
    if (handleOut)
	*handleOut = BObolFeatureHandle();
    if (!controller)
	return NULL;
    BObolFeatureHandle handle = controller->features().find(
	navgizmo_feature_name, BOBOL_FEATURE_SCOPE_LOCAL);
    if (!handle.isValid())
	return NULL;
    BObolFeatureRecord record;
    if (!controller->features().record(handle, record) ||
	bu_strcmp(record.owner.ownerRole.getString(), navgizmo_owner_role) != 0 ||
	!record.owner.ownerToken ||
	record.owner.callbackUserData != record.owner.ownerToken)
	return NULL;
    if (handleOut)
	*handleOut = handle;
    return static_cast<NavgizmoState *>(
	const_cast<void *>(record.owner.ownerToken));
}

static SoBRLNavigationGizmo *
navgizmo_current_gizmo(NavgizmoState *state,
    BObolFeatureHandle *handleOut = NULL)
{
    BObolFeatureHandle handle;
    if (!state || navgizmo_state(state->controller, &handle) != state)
	return NULL;
    SoNode *node = state->controller->features().node(handle);
    if (!node || !node->isOfType(SoBRLNavigationGizmo::getClassTypeId()))
	return NULL;
    if (handleOut)
	*handleOut = handle;
    return static_cast<SoBRLNavigationGizmo *>(node);
}

static BObolFeatureOwner
navgizmo_owner(NavgizmoState *state)
{
    BObolFeatureOwner owner;
    owner.ownerToken = state;
    owner.ownerId = "libged::navgizmo";
    owner.ownerRole = navgizmo_owner_role;
    owner.resultCallback = navgizmo_feature_result;
    owner.callbackUserData = state;
    return owner;
}

static BObolOverlayInfo
navgizmo_overlay(NavgizmoState *state)
{
    BObolOverlayInfo overlay;
    overlay.isOverlay = TRUE;
    overlay.ownerToken = state;
    overlay.role = BObolOverlayRole::Screen;
    overlay.overlayClass = BObolOverlayClass::Faceplate;
    overlay.lifecycle = BObolOverlayLifecycle::PerView;
    overlay.order = BObolOverlayOrder::Screen;
    overlay.sortOrder = navgizmo_overlay_order;
    return overlay;
}

struct NavgizmoPublicationCommit {
    NavgizmoState *state;
    NavgizmoSnapshot snapshot;
    SbBool committed;
};

static void
navgizmo_commit_snapshot(void *context) noexcept
{
    NavgizmoPublicationCommit *commit =
	static_cast<NavgizmoPublicationCommit *>(context);
    if (commit && commit->state) {
	commit->state->snapshot = commit->snapshot;
	commit->committed = TRUE;
    }
}

static int
navgizmo_publish_impl(NavgizmoState *state, NavgizmoSnapshot next,
    SbBool forceSuccessor)
{
    if (!state || !state->controller)
	return 0;
    next.cameraOrientation = navgizmo_camera_orientation(state);

    const SbBool hasCurrent = navgizmo_current_gizmo(state) != NULL;
    if (!hasCurrent && state->controller->features().exists(
	    navgizmo_feature_name, BOBOL_FEATURE_SCOPE_LOCAL))
	return 0;
    if (!forceSuccessor && hasCurrent &&
	navgizmo_presentations_equal(state->snapshot, next)) {
	state->snapshot = next;
	return 1;
    }

    SoBRLNavigationGizmo *gizmo = NULL;
    try {
	gizmo = new (std::nothrow) SoBRLNavigationGizmo;
    } catch (...) {
	return 0;
    }
    if (!gizmo)
	return 0;
    gizmo->ref();
    NavgizmoPublicationCommit commit = {state, next, FALSE};
    try {
	gizmo->visible = next.visible;
	gizmo->style = next.style;
	gizmo->hoverPart = next.hoverPart;
	gizmo->activePart = next.activePart;
	gizmo->setCameraOrientationSnapshot(next.cameraOrientation);
	(void)gizmo->rebuildGeometry();

	BObolFeatureStyle style;
	style.hasVisible = TRUE;
	style.visible = next.visible;
	style.hasSelectable = TRUE;
	style.selectable = FALSE;
	style.hud = TRUE;

	BObolFeaturePublication feature;
	feature.name = navgizmo_feature_name;
	feature.kind = BObolFeatureKind::CustomNode;
	feature.scope = BObolFeatureScope::Local;
	feature.style = style;
	feature.owner = navgizmo_owner(state);
	feature.overlay = navgizmo_overlay(state);
	feature.customNode = gizmo;

	BObolFeatureStorePublication storePublication;
	storePublication.store = &state->controller->features();
	storePublication.features.push_back(feature);
	std::vector<BObolFeatureStorePublication> publications;
	publications.push_back(storePublication);
	const SbBool published =
	    BObolFeatureStore::applyCoordinatedPublications(publications,
		navgizmo_commit_snapshot, &commit);
	gizmo->unref();
	return published ? 1 : 0;
    } catch (...) {
	gizmo->unref();
	/* Observer failure follows commit.  Preparation failure leaves the
	 * preceding node and authoritative interaction state in place. */
	return commit.committed ? 1 : 0;
    }
}

static int
navgizmo_publish(NavgizmoState *state, const NavgizmoSnapshot &next,
    SbBool forceSuccessor = FALSE)
{
    try {
	return navgizmo_publish_impl(state, next, forceSuccessor);
    } catch (...) {
	/* Coin input and sensor callbacks must not leak C++ exceptions. */
	return 0;
    }
}

static void
navgizmo_camera_changed(void *userData, SoSensor *UNUSED(sensor))
{
    NavgizmoState *state = static_cast<NavgizmoState *>(userData);
    if (!state || state->suppressCameraPublication)
	return;
    /* A hidden node has no orientation-dependent geometry.  Its next visible
	 * successor snapshots the then-current camera without idle frame churn. */
    if (!state->snapshot.visible)
	return;
    (void)navgizmo_publish(state, state->snapshot);
}

static void
navgizmo_camera_replaced(void *userData, SoSensor *UNUSED(sensor))
{
    NavgizmoState *state = static_cast<NavgizmoState *>(userData);
    if (!navgizmo_track_camera(state) || state->suppressCameraPublication ||
	!state->snapshot.visible)
	return;
    /* Camera identity is part of the retained observation generation even
     * when the replacement happens to have the same orientation. */
    (void)navgizmo_publish(state, state->snapshot, TRUE);
}

static int
navgizmo_dimensions(NavgizmoState *state, int &width, int &height)
{
    width = 0;
    height = 0;
    if (!state || !state->viewContext)
	return 0;
    const struct bv_context *context =
	reinterpret_cast<const struct bv_context *>(state->viewContext);
    width = bv_context_width_get(context);
    height = bv_context_height_get(context);
    if ((width <= 0 || height <= 0) && state->controller) {
	const SbVec2s viewport =
	    state->controller->getViewportRegion().getViewportSizePixels();
	width = static_cast<int>(viewport[0]);
	height = static_cast<int>(viewport[1]);
    }
    return width > 0 && height > 0 ? 1 : 0;
}

static SoBRLNavigationGizmo::Part
navgizmo_hit(NavgizmoState *state, const BObolInputEvent *event)
{
    if (!state || !event)
	return SoBRLNavigationGizmo::PART_NONE;
    if (!navgizmo_publish(state, state->snapshot))
	return SoBRLNavigationGizmo::PART_NONE;
    SoBRLNavigationGizmo *gizmo = navgizmo_current_gizmo(state);
    if (!gizmo)
	return SoBRLNavigationGizmo::PART_NONE;
    int width = 0;
    int height = 0;
    if (!navgizmo_dimensions(state, width, height))
	return SoBRLNavigationGizmo::PART_NONE;
    double pixelRatio = 1.0;
    struct bv_display_property_value value = BV_DISPLAY_PROPERTY_VALUE_INIT;
    if (state->endpoint && bobol_display_endpoint_property_get(state->endpoint,
	    "endpoint.device_pixel_ratio", &value) == BV_DISPLAY_PROPERTY_OK &&
	value.type == BV_DISPLAY_PROPERTY_DOUBLE && value.double_value > 0.0)
	pixelRatio = value.double_value;
    return gizmo->hitTest(
	static_cast<int>(std::lround(event->x * pixelRatio)),
	static_cast<int>(std::lround(event->y * pixelRatio)), width, height);
}

static int
navgizmo_view_update(NavgizmoState *state)
{
    if (!state || !state->viewContext)
	return 0;
    state->suppressCameraPublication = 1;
    const int updated = ged_view_context_update(state->viewContext);
    /* A GED host callback may synchronize this already, but headless and
	 * embedding clients are permitted to have no callback. */
    const int synchronized = updated && (!state->endpoint ||
	bobol_display_endpoint_view_sync(state->endpoint,
	    state->viewContext));
    state->suppressCameraPublication = 0;
    return synchronized;
}

static int
navgizmo_rotate(NavgizmoState *state, const BObolInputEvent *event,
    const NavgizmoSnapshot &next)
{
    if (!state || !state->viewContext || !event)
	return 0;
    struct bv_context *context =
	reinterpret_cast<struct bv_context *>(state->viewContext);
    struct bv *view = bv_context_view(context);
    point_t center = VINIT_ZERO;
    if (!view || !bv_center_get(center, view) ||
	!bv_mouse_delta_adjust(view, event->dx, event->dy, center,
	    BV_ADJUST_ROT))
	return 0;
    return navgizmo_view_update(state) && navgizmo_publish(state, next);
}

static int
navgizmo_orient(NavgizmoState *state, SoBRLNavigationGizmo::Part part)
{
    if (!state || !state->viewContext)
	return 0;
    float azimuth = 0.0f;
    float elevation = 0.0f;
    if (!SoBRLNavigationGizmo::partAet(part, azimuth, elevation))
	return 0;

    struct bv_context *context =
	reinterpret_cast<struct bv_context *>(state->viewContext);
    struct bv *view = bv_context_view(context);
    vect_t aet;
    VSET(aet, static_cast<fastf_t>(azimuth),
	static_cast<fastf_t>(elevation), 0.0);
    if (!view || !bv_aet_set(view, aet))
	return 0;

    return navgizmo_view_update(state);
}

static int
navgizmo_input(void *userData, BObolInputAction action,
    const BObolInputEvent *event)
{
    NavgizmoState *state = static_cast<NavgizmoState *>(userData);
    if (!state || !navgizmo_current_gizmo(state) || !event)
	return BOBOL_INPUT_RESULT_UNHANDLED;

    if (action == NAVGIZMO_ACTION_PRESS) {
	const SoBRLNavigationGizmo::Part part = navgizmo_hit(state, event);
	if (part == SoBRLNavigationGizmo::PART_NONE)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	NavgizmoSnapshot next = state->snapshot;
	next.pressedPart = part;
	next.pressX = event->x;
	next.pressY = event->y;
	next.dragging = 1;
	next.moved = 0;
	next.hoverPart = part;
	next.activePart = part;
	(void)navgizmo_publish(state, next);
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (action == NAVGIZMO_ACTION_MOTION) {
	const SoBRLNavigationGizmo::Part part = navgizmo_hit(state, event);
	if (!state->snapshot.dragging) {
	    NavgizmoSnapshot next = state->snapshot;
	    next.hoverPart = part;
	    (void)navgizmo_publish(state, next);
	    return part == SoBRLNavigationGizmo::PART_NONE ?
		BOBOL_INPUT_RESULT_UNHANDLED : BOBOL_INPUT_RESULT_HANDLED;
	}

	NavgizmoSnapshot next = state->snapshot;
	if (std::abs(event->x - next.pressX) >
		navgizmo_drag_threshold_pixels ||
	    std::abs(event->y - next.pressY) >
		navgizmo_drag_threshold_pixels)
	    next.moved = 1;
	(void)navgizmo_rotate(state, event, next);
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (action == NAVGIZMO_ACTION_RELEASE) {
	if (!state->snapshot.dragging)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	const SoBRLNavigationGizmo::Part pressed = state->snapshot.pressedPart;
	const SoBRLNavigationGizmo::Part hover = navgizmo_hit(state, event);
	const int click = !state->snapshot.moved && hover == pressed;
	if (click)
	    (void)navgizmo_orient(state, pressed);
	NavgizmoSnapshot next = state->snapshot;
	next.dragging = 0;
	next.moved = 0;
	next.pressedPart = SoBRLNavigationGizmo::PART_NONE;
	next.activePart = SoBRLNavigationGizmo::PART_NONE;
	next.hoverPart = hover;
	(void)navgizmo_publish(state, next);
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    return BOBOL_INPUT_RESULT_UNHANDLED;
}

static int
navgizmo_input_install(NavgizmoState *state)
{
    if (!state || !state->endpoint)
	return 0;
    if (state->inputInstalled)
	return 1;
    static const unsigned int allModifiers = BOBOL_INPUT_MOD_SHIFT |
	BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META;
    static const BObolInputBinding bindings[] = {
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0, 0,
	 allModifiers, navgizmo_input_priority, NAVGIZMO_ACTION_PRESS},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0, 0,
	 0, navgizmo_input_priority, NAVGIZMO_ACTION_RELEASE},
	{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY, 0,
	 0, navgizmo_input_priority, NAVGIZMO_ACTION_MOTION}
    };
    static const BObolInputActionLayer layer = {
	"libged-navigation-gizmo", bindings,
	sizeof(bindings) / sizeof(bindings[0]), navgizmo_input
    };
    state->inputInstalled = bobol_display_endpoint_input_action_layer_set(
	state->endpoint, &layer, state, state) ? 1 : 0;
    return state->inputInstalled;
}

static void
navgizmo_input_clear(NavgizmoState *state)
{
    if (!state || !state->endpoint || !state->inputInstalled)
	return;
    (void)bobol_display_endpoint_input_action_layer_clear_if(
	state->endpoint, state);
    state->inputInstalled = 0;
}

static void
navgizmo_feature_result(const BObolCommandResult &result, void *userData)
{
    if (result.status != BObolCommandResultStatus::Removed)
	return;
    NavgizmoState *state = static_cast<NavgizmoState *>(userData);
    if (!state)
	return;
    navgizmo_input_clear(state);
    state->controller = NULL;
    state->endpoint = NULL;
    state->viewContext = NULL;
    delete state;
}

static int
navgizmo_create(struct ged *gedp, struct ged_view_context *viewContext,
    bobol_display_endpoint_t *endpoint, BObolViewController *controller,
    SbBool visible, SoBRLNavigationGizmo::Style style)
{
    if (!controller || controller->features().exists(navgizmo_feature_name,
	BOBOL_FEATURE_SCOPE_LOCAL))
	return BRLCAD_ERROR;

    NavgizmoState *state = NULL;
    try {
	state = new NavgizmoState;
    } catch (...) {
	return BRLCAD_ERROR;
    }
    state->viewContext = viewContext;
    state->endpoint = endpoint;
    state->controller = controller;
    /* The field sensor observes orientation within one camera generation.
     * Root publication is the completion edge for replacing that generation. */
    state->cameraRootSensor.attach(controller->getViewport()->getRoot());

    NavgizmoSnapshot next = state->snapshot;
    next.visible = visible;
    next.style = style;
    if (visible && !navgizmo_input_install(state)) {
	delete state;
	bu_vls_printf(gedp->ged_result_str,
	    "unable to install the navigation gizmo input layer");
	return BRLCAD_ERROR;
    }
    if (!navgizmo_publish(state, next)) {
	navgizmo_input_clear(state);
	delete state;
	bu_vls_printf(gedp->ged_result_str,
	    "unable to publish the navigation gizmo");
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
navgizmo_enable(struct ged *gedp, struct ged_view_context *viewContext,
    bobol_display_endpoint_t *endpoint, BObolViewController *controller)
{
    NavgizmoState *existing = navgizmo_state(controller);
    if (existing) {
	if (existing->snapshot.visible) {
	    if (navgizmo_publish(existing, existing->snapshot))
		return BRLCAD_OK;
	    bu_vls_printf(gedp->ged_result_str,
		"unable to synchronize the navigation gizmo");
	    return BRLCAD_ERROR;
	}
	if (!navgizmo_input_install(existing)) {
	    bu_vls_printf(gedp->ged_result_str,
		"unable to install the navigation gizmo input layer");
	    return BRLCAD_ERROR;
	}
	NavgizmoSnapshot next = existing->snapshot;
	next.visible = TRUE;
	next.hoverPart = SoBRLNavigationGizmo::PART_NONE;
	next.activePart = SoBRLNavigationGizmo::PART_NONE;
	next.pressedPart = SoBRLNavigationGizmo::PART_NONE;
	next.dragging = 0;
	next.moved = 0;
	if (!navgizmo_publish(existing, next)) {
	    navgizmo_input_clear(existing);
	    bu_vls_printf(gedp->ged_result_str,
		"unable to publish the navigation gizmo");
	    return BRLCAD_ERROR;
	}
	return BRLCAD_OK;
    }
    if (controller->features().exists(navgizmo_feature_name,
	BOBOL_FEATURE_SCOPE_LOCAL)) {
	bu_vls_printf(gedp->ged_result_str,
	    "a different owner already published %s", navgizmo_feature_name);
	return BRLCAD_ERROR;
    }
    return navgizmo_create(gedp, viewContext, endpoint, controller, TRUE,
	SoBRLNavigationGizmo::CUBE);
}

static int
navgizmo_disable(BObolViewController *controller)
{
    NavgizmoState *state = navgizmo_state(controller);
    if (!state)
	return BRLCAD_OK;
    if (!state->snapshot.visible)
	return BRLCAD_OK;
    NavgizmoSnapshot next = state->snapshot;
    next.visible = FALSE;
    next.hoverPart = SoBRLNavigationGizmo::PART_NONE;
    next.activePart = SoBRLNavigationGizmo::PART_NONE;
    next.pressedPart = SoBRLNavigationGizmo::PART_NONE;
    next.dragging = 0;
    next.moved = 0;
    if (!navgizmo_publish(state, next))
	return BRLCAD_ERROR;
    navgizmo_input_clear(state);
    return BRLCAD_OK;
}

static void
navgizmo_usage(struct ged *gedp, const char *command)
{
    (void)command;
    bu_vls_printf(gedp->ged_result_str,
	"Usage: view faceplate navgizmo [0|1|off|on|toggle|style [cube|circles]]\n");
}

} // namespace

extern "C" int
ged_navgizmo_core(struct ged *gedp, int argc, const char *argv[])
{
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);
    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct ged_view_context *viewContext = ged_view_active_ctx(gedp);
    bobol_display_endpoint_t *endpoint =
	ged_plugin_obol_endpoint_get(viewContext);
    BObolViewController *controller =
	ged_plugin_obol_view_controller(viewContext);
    if (!endpoint || !controller) {
	bu_vls_printf(gedp->ged_result_str,
	    "%s requires an active Obol display endpoint", argv[0]);
	return BRLCAD_ERROR;
    }

    NavgizmoState *state = navgizmo_state(controller);
    const int enabled = state && state->snapshot.visible ? 1 : 0;
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "%d", enabled);
	return BRLCAD_OK;
    }
    if (argc < 2 || argc > 3) {
	navgizmo_usage(gedp, argv[0]);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(argv[1], "help") || BU_STR_EQUAL(argv[1], "-h") ||
	BU_STR_EQUAL(argv[1], "--help")) {
	navgizmo_usage(gedp, argv[0]);
	return GED_HELP;
    }
    if (BU_STR_EQUAL(argv[1], "style")) {
	if (argc == 2) {
	    const int current = state ? state->snapshot.style :
		SoBRLNavigationGizmo::CUBE;
	    bu_vls_printf(gedp->ged_result_str, "%s",
		current == SoBRLNavigationGizmo::CIRCLES ? "circles" : "cube");
	    return BRLCAD_OK;
	}
	SoBRLNavigationGizmo::Style requested;
	if (BU_STR_EQUAL(argv[2], "cube"))
	    requested = SoBRLNavigationGizmo::CUBE;
	else if (BU_STR_EQUAL(argv[2], "circles"))
	    requested = SoBRLNavigationGizmo::CIRCLES;
	else {
	    navgizmo_usage(gedp, argv[0]);
	    return BRLCAD_ERROR;
	}
	if (!state) {
	    if (controller->features().exists(navgizmo_feature_name,
		    BOBOL_FEATURE_SCOPE_LOCAL)) {
		bu_vls_printf(gedp->ged_result_str,
		    "a different owner already published %s",
		    navgizmo_feature_name);
		return BRLCAD_ERROR;
	    }
	    return navgizmo_create(gedp, viewContext, endpoint, controller,
		FALSE, requested);
	}
	NavgizmoSnapshot next = state->snapshot;
	next.style = requested;
	if (!navgizmo_publish(state, next)) {
	    bu_vls_printf(gedp->ged_result_str,
		"unable to publish the navigation gizmo style");
	    return BRLCAD_ERROR;
	}
	return BRLCAD_OK;
    }
    if (argc != 2) {
	navgizmo_usage(gedp, argv[0]);
	return BRLCAD_ERROR;
    }
    if (BU_STR_EQUAL(argv[1], "toggle"))
	return enabled ? navgizmo_disable(controller) :
	    navgizmo_enable(gedp, viewContext, endpoint, controller);
    if (BU_STR_EQUAL(argv[1], "1") || BU_STR_EQUAL(argv[1], "on"))
	return navgizmo_enable(gedp, viewContext, endpoint, controller);
    if (BU_STR_EQUAL(argv[1], "0") || BU_STR_EQUAL(argv[1], "off"))
	return navgizmo_disable(controller);

    navgizmo_usage(gedp, argv[0]);
    return BRLCAD_ERROR;
}

#define GED_NAVGIZMO_COMMANDS(X, XID) \
    XID(navgizmo, "view.faceplate.navgizmo", ged_navgizmo_core, \
	GED_CMD_UPDATE_VIEW)

GED_DECLARE_COMMAND_SET(GED_NAVGIZMO_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_navgizmo", 1, GED_NAVGIZMO_COMMANDS)

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
