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

#include <Inventor/nodes/SoCamera.h>

#include <cmath>
#include <new>

namespace {

static const char *navgizmo_feature_name = "faceplate::navigation_gizmo";
static const char *navgizmo_owner_role = "libged.navigation-gizmo";

enum NavgizmoAction {
    NAVGIZMO_ACTION_PRESS = 0x4e470001u,
    NAVGIZMO_ACTION_RELEASE = 0x4e470002u,
    NAVGIZMO_ACTION_MOTION = 0x4e470003u
};

struct NavgizmoState {
    struct ged_view_context *viewContext;
    bobol_display_endpoint_t *endpoint;
    BObolViewController *controller;
    SoBRLNavigationGizmo *gizmo;
    SoBRLNavigationGizmo::Part pressedPart;
    int pressX;
    int pressY;
    int dragging;
    int moved;

    NavgizmoState(void) :
	viewContext(NULL), endpoint(NULL), controller(NULL), gizmo(NULL),
	pressedPart(SoBRLNavigationGizmo::PART_NONE), pressX(0), pressY(0),
	dragging(0), moved(0)
    {
    }
};

static void
navgizmo_request_frame(NavgizmoState *state, const char *reason)
{
    if (!state)
	return;
    if (state->endpoint) {
	(void)bobol_display_endpoint_request_presentation_frame(
	    state->endpoint, reason);
	return;
    }
    if (state->controller)
	state->controller->requestPresentationRender(reason);
}

static void
navgizmo_camera_refresh(NavgizmoState *state)
{
    if (!state || !state->controller || !state->gizmo)
	return;
    SoCamera *camera = state->controller->getCamera();
    if (camera != state->gizmo->getCamera())
	state->gizmo->setCamera(camera);
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
    if (!state || !state->gizmo || !event)
	return SoBRLNavigationGizmo::PART_NONE;
    int width = 0;
    int height = 0;
    if (!navgizmo_dimensions(state, width, height))
	return SoBRLNavigationGizmo::PART_NONE;
    navgizmo_camera_refresh(state);
    double pixelRatio = 1.0;
    struct bv_display_property_value value = BV_DISPLAY_PROPERTY_VALUE_INIT;
    if (state->endpoint && bobol_display_endpoint_property_get(state->endpoint,
	    "endpoint.device_pixel_ratio", &value) == BV_DISPLAY_PROPERTY_OK &&
	value.type == BV_DISPLAY_PROPERTY_DOUBLE && value.double_value > 0.0)
	pixelRatio = value.double_value;
    return state->gizmo->hitTest(
	static_cast<int>(std::lround(event->x * pixelRatio)),
	static_cast<int>(std::lround(event->y * pixelRatio)), width, height);
}

static int
navgizmo_rotate(NavgizmoState *state, const BObolInputEvent *event)
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

    (void)ged_view_context_update(state->viewContext);
    navgizmo_camera_refresh(state);
    navgizmo_request_frame(state, "navigation-gizmo-drag");
    return 1;
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

    (void)ged_view_context_update(state->viewContext);
    navgizmo_camera_refresh(state);
    navgizmo_request_frame(state, "navigation-gizmo-orient");
    return 1;
}

static int
navgizmo_input(void *userData, BObolInputAction action,
    const BObolInputEvent *event)
{
    NavgizmoState *state = static_cast<NavgizmoState *>(userData);
    if (!state || !state->gizmo || !event)
	return BOBOL_INPUT_RESULT_UNHANDLED;

    if (action == NAVGIZMO_ACTION_PRESS) {
	const SoBRLNavigationGizmo::Part part = navgizmo_hit(state, event);
	if (part == SoBRLNavigationGizmo::PART_NONE)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	state->pressedPart = part;
	state->pressX = event->x;
	state->pressY = event->y;
	state->dragging = 1;
	state->moved = 0;
	state->gizmo->setHoverPart(part);
	state->gizmo->setActivePart(part);
	navgizmo_request_frame(state, "navigation-gizmo-press");
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (action == NAVGIZMO_ACTION_MOTION) {
	const SoBRLNavigationGizmo::Part part = navgizmo_hit(state, event);
	if (!state->dragging) {
	    const int changed = state->gizmo->hoverPart.getValue() !=
		static_cast<int>(part);
	    state->gizmo->setHoverPart(part);
	    if (changed)
		navgizmo_request_frame(state, "navigation-gizmo-hover");
	    return part == SoBRLNavigationGizmo::PART_NONE ?
		BOBOL_INPUT_RESULT_UNHANDLED : BOBOL_INPUT_RESULT_HANDLED;
	}

	if (std::abs(event->x - state->pressX) > 2 ||
	    std::abs(event->y - state->pressY) > 2)
	    state->moved = 1;
	(void)navgizmo_rotate(state, event);
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (action == NAVGIZMO_ACTION_RELEASE) {
	if (!state->dragging)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	const SoBRLNavigationGizmo::Part pressed = state->pressedPart;
	const SoBRLNavigationGizmo::Part hover = navgizmo_hit(state, event);
	const int click = !state->moved && hover == pressed;
	state->dragging = 0;
	state->pressedPart = SoBRLNavigationGizmo::PART_NONE;
	state->gizmo->setActivePart(SoBRLNavigationGizmo::PART_NONE);
	state->gizmo->setHoverPart(hover);
	if (!click || !navgizmo_orient(state, pressed))
	    navgizmo_request_frame(state, "navigation-gizmo-release");
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    return BOBOL_INPUT_RESULT_UNHANDLED;
}

static int
navgizmo_input_install(NavgizmoState *state)
{
    if (!state || !state->endpoint)
	return 0;
    static const unsigned int allModifiers = BOBOL_INPUT_MOD_SHIFT |
	BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META;
    static const BObolInputBinding bindings[] = {
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0, 0,
	 allModifiers, 1000, NAVGIZMO_ACTION_PRESS},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0, 0,
	 0, 1000, NAVGIZMO_ACTION_RELEASE},
	{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY, 0,
	 0, 1000, NAVGIZMO_ACTION_MOTION}
    };
    static const BObolInputActionLayer layer = {
	"libged-navigation-gizmo", bindings,
	sizeof(bindings) / sizeof(bindings[0]), navgizmo_input
    };
    return bobol_display_endpoint_input_action_layer_set(state->endpoint,
	&layer, state, state);
}

static void
navgizmo_feature_result(const BObolCommandResult &result, void *userData)
{
    if (result.status != BObolCommandResultStatus::Removed)
	return;
    NavgizmoState *state = static_cast<NavgizmoState *>(userData);
    if (!state)
	return;
    if (state->endpoint)
	(void)bobol_display_endpoint_input_action_layer_clear_if(
	    state->endpoint, state);
    state->gizmo = NULL;
    state->controller = NULL;
    state->endpoint = NULL;
    state->viewContext = NULL;
    delete state;
}

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

static int
navgizmo_enable(struct ged *gedp, struct ged_view_context *viewContext,
    bobol_display_endpoint_t *endpoint, BObolViewController *controller)
{
    NavgizmoState *existing = navgizmo_state(controller);
    if (existing) {
	existing->gizmo->visible = TRUE;
	existing->gizmo->setHoverPart(SoBRLNavigationGizmo::PART_NONE);
	existing->gizmo->setActivePart(SoBRLNavigationGizmo::PART_NONE);
	(void)existing->gizmo->rebuildGeometry();
	if (!navgizmo_input_install(existing)) {
	    existing->gizmo->visible = FALSE;
	    (void)existing->gizmo->rebuildGeometry();
	    bu_vls_printf(gedp->ged_result_str,
		"unable to install the navigation gizmo input layer");
	    return BRLCAD_ERROR;
	}
	navgizmo_request_frame(existing, "navigation-gizmo-enable");
	return BRLCAD_OK;
    }
    if (controller->features().exists(navgizmo_feature_name,
	BOBOL_FEATURE_SCOPE_LOCAL)) {
	bu_vls_printf(gedp->ged_result_str,
	    "a different owner already published %s", navgizmo_feature_name);
	return BRLCAD_ERROR;
    }

    NavgizmoState *state = new (std::nothrow) NavgizmoState;
    if (!state)
	return BRLCAD_ERROR;
    state->viewContext = viewContext;
    state->endpoint = endpoint;
    state->controller = controller;

    SoBRLNavigationGizmo *gizmo = new SoBRLNavigationGizmo;
    gizmo->ref();
    state->gizmo = gizmo;
    gizmo->setCamera(controller->getCamera());
    (void)gizmo->rebuildGeometry();

    BObolFeatureStyle style;
    style.hasVisible = TRUE;
    style.visible = TRUE;
    style.hasSelectable = TRUE;
    style.selectable = FALSE;
    style.hud = TRUE;

    BObolFeatureOwner owner;
    owner.ownerToken = state;
    owner.ownerId = "libged::navgizmo";
    owner.ownerRole = navgizmo_owner_role;
    owner.resultCallback = navgizmo_feature_result;
    owner.callbackUserData = state;

    const BObolFeatureHandle handle = controller->features().publishCustomNode(
	navgizmo_feature_name, BObolFeatureScope::Local, gizmo, &style, &owner);
    gizmo->unref();
    if (!handle.isValid()) {
	delete state;
	return BRLCAD_ERROR;
    }

    if (!navgizmo_input_install(state)) {
	(void)controller->features().remove(handle);
	bu_vls_printf(gedp->ged_result_str,
	    "unable to install the navigation gizmo input layer");
	return BRLCAD_ERROR;
    }

    BObolOverlayInfo overlay;
    overlay.isOverlay = TRUE;
    overlay.ownerToken = state;
    overlay.role = BObolOverlayRole::Screen;
    overlay.overlayClass = BObolOverlayClass::Faceplate;
    overlay.lifecycle = BObolOverlayLifecycle::PerView;
    overlay.order = BObolOverlayOrder::Screen;
    overlay.sortOrder = 500;
    (void)controller->features().setOverlayInfo(handle, overlay);
    navgizmo_request_frame(state, "navigation-gizmo-enable");
    return BRLCAD_OK;
}

static int
navgizmo_disable(BObolViewController *controller)
{
    NavgizmoState *state = navgizmo_state(controller);
    if (!state)
	return BRLCAD_OK;
    if (state->endpoint)
	(void)bobol_display_endpoint_input_action_layer_clear_if(
	    state->endpoint, state);
    state->dragging = 0;
    state->moved = 0;
    state->pressedPart = SoBRLNavigationGizmo::PART_NONE;
    state->gizmo->hoverPart = SoBRLNavigationGizmo::PART_NONE;
    state->gizmo->activePart = SoBRLNavigationGizmo::PART_NONE;
    state->gizmo->visible = FALSE;
    (void)state->gizmo->rebuildGeometry();
    navgizmo_request_frame(state, "navigation-gizmo-disable");
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
    const int enabled = state && state->gizmo->visible.getValue() ? 1 : 0;
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
	    const int current = state ? state->gizmo->style.getValue() :
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
	const int wasEnabled = enabled;
	if (!state) {
	    if (navgizmo_enable(gedp, viewContext, endpoint, controller) !=
		BRLCAD_OK)
		return BRLCAD_ERROR;
	    state = navgizmo_state(controller);
	}
	if (!state)
	    return BRLCAD_ERROR;
	state->gizmo->style = requested;
	(void)state->gizmo->rebuildGeometry();
	navgizmo_request_frame(state, "navigation-gizmo-style");
	if (!wasEnabled)
	    return navgizmo_disable(controller);
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
