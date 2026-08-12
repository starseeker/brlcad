/*          T E S T _ N A V I G A T I O N _ G I Z M O _ P L U G I N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BNavigationGizmo.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "bu/app.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "bv.h"
#include "ged.h"
#include "ged/plugin/obol.h"

#include <Inventor/SbMatrix.h>
#include <Inventor/nodes/SoNode.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

static int failures = 0;

#define CHECK(_condition, _message) do { \
    if (!(_condition)) { \
	bu_log("FAIL [%s:%d] %s\n", __FILE__, __LINE__, (_message)); \
	failures++; \
    } \
} while (0)

static int fallbackPresses = 0;

static int
fallback_input(void *, BObolInputAction, const BObolInputEvent *)
{
    fallbackPresses++;
    return BOBOL_INPUT_RESULT_HANDLED;
}

static int
run_command(struct ged *gedp, const char *argument,
    const char *value = NULL)
{
    const char *argv[6] = {"view", "faceplate", "navgizmo", argument,
	value, NULL};
    bu_vls_trunc(gedp->ged_result_str, 0);
    return ged_exec(gedp, value ? 5 : (argument ? 4 : 3), argv);
}

int
main(int argc, const char **argv)
{
    bu_setprogname(argv[0]);
    if (argc != 2) {
	bu_log("Usage: %s <directory-containing-moss.g>\n", argv[0]);
	return EXIT_FAILURE;
    }

    struct bu_vls path = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&path, "%s/moss.g", argv[1]);
    struct ged *gedp = ged_open("db", bu_vls_cstr(&path), 1);
    CHECK(gedp != NULL, "test database must open");
    if (!gedp) {
	bu_vls_free(&path);
	return EXIT_FAILURE;
    }

    struct ged_view_context *viewContext = ged_view_active_ctx(gedp);
    if (!viewContext) {
	viewContext = ged_view_context_create();
	CHECK(viewContext != NULL, "test view must be created");
	if (viewContext) {
	    CHECK(ged_view_set_context_add(ged_view_set_ctx(gedp), viewContext),
		"test view must join the GED view set");
	    ged_view_context_owned_add(gedp, viewContext);
	    ged_view_active_ctx_set(gedp, viewContext);
	}
    }
    CHECK(viewContext != NULL, "an active GED view is required");
    CHECK(ged_view_context_display_endpoint_ensure(viewContext),
	"active view must acquire an Obol endpoint");
    bobol_display_endpoint_t *endpoint =
	ged_plugin_obol_endpoint_get(viewContext);
    BObolViewController *controller =
	ged_plugin_obol_view_controller(viewContext);
    CHECK(endpoint != NULL && controller != NULL,
	"installed plugin extension must expose the active endpoint/controller");
    if (!viewContext || !endpoint || !controller) {
	ged_close(gedp);
	bu_vls_free(&path);
	return EXIT_FAILURE;
    }

    CHECK(bv_context_dimensions_set(
	reinterpret_cast<struct bv_context *>(viewContext), 800, 600),
	"test view dimensions must be set");
    CHECK(bobol_display_endpoint_resize(endpoint, 800, 600, 1.0),
	"test endpoint dimensions must be set");
    CHECK(bobol_display_endpoint_view_sync(endpoint, viewContext),
	"test camera must synchronize from BRL-CAD view state");

    CHECK(!ged_cmd_exists("navgizmo") &&
	ged_cmd_exists("view.faceplate.navgizmo"),
	"plugin command must live only in the qualified view namespace");
    const char * const *topCommands = NULL;
    const size_t topCommandCount = ged_cmd_list(&topCommands);
    bool qualifiedLeaked = false;
    for (size_t i = 0; i < topCommandCount; i++)
	qualifiedLeaked = qualifiedLeaked ||
	    (topCommands[i] && strchr(topCommands[i], '.') != NULL);
    CHECK(!qualifiedLeaked,
	"qualified plugin entries must not pollute GED's top-level command list");
    const char *helpArgv[] = {"view", "faceplate", "-h", NULL};
    CHECK(ged_exec(gedp, 3, helpArgv) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "navgizmo") != NULL,
	"faceplate help must discover dynamically registered plugin subcommands");

    static const BObolInputBinding fallbackBindings[] = {
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0, 0, 0, 500,
	 0x46414c4cu}
    };
    static const BObolInputActionLayer fallbackLayer = {
	"navigation-gizmo-test-fallback", fallbackBindings,
	sizeof(fallbackBindings) / sizeof(fallbackBindings[0]), fallback_input
    };
    int fallbackOwner = 0;
    CHECK(bobol_display_endpoint_input_action_layer_set(endpoint,
	&fallbackLayer, &fallbackOwner, NULL),
	"pre-existing plugin input layer must install");

    CHECK(run_command(gedp, NULL) == BRLCAD_OK &&
	BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "0"),
	"dynamic navgizmo command should initially report disabled");
    CHECK(run_command(gedp, "style") == BRLCAD_OK &&
	BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "cube"),
	"navigation gizmo should initially report cube style");
    CHECK(run_command(gedp, "style", "circles") == BRLCAD_OK &&
	run_command(gedp, NULL) == BRLCAD_OK &&
	BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "0"),
	"setting style while hidden must not implicitly show the gizmo");
    CHECK(run_command(gedp, "style") == BRLCAD_OK &&
	BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "circles"),
	"circle style selection must persist while the gizmo is hidden");
    CHECK(run_command(gedp, "on") == BRLCAD_OK,
	"dynamic navgizmo command should enable the feature");

    BObolFeatureHandle handle = controller->features().find(
	"faceplate::navigation_gizmo", BOBOL_FEATURE_SCOPE_LOCAL);
    SoNode *node = controller->features().node(handle);
    CHECK(handle.isValid() && node &&
	node->isOfType(SoBRLNavigationGizmo::getClassTypeId()),
	"plugin must publish the concrete retained gizmo through FeatureStore");
    CHECK(controller->getFramebufferOverlayRoot()->findChild(node) >= 0,
	"screen faceplate must attach to the post-CAD framebuffer overlay root");
    SoBRLNavigationGizmo *gizmo = node &&
	node->isOfType(SoBRLNavigationGizmo::getClassTypeId()) ?
	static_cast<SoBRLNavigationGizmo *>(node) : NULL;
    CHECK(gizmo && gizmo->style.getValue() == SoBRLNavigationGizmo::CIRCLES,
	"published gizmo must realize the selected circles style");

    if (gizmo) {
	SoBRLNavigationGizmo::Part part = SoBRLNavigationGizmo::PART_NONE;
	int targetX = 740;
	int targetY = 60;
	float expectedAzimuth = 0.0f;
	float expectedElevation = 0.0f;
	for (int y = 12; y <= 108 &&
	    part == SoBRLNavigationGizmo::PART_NONE; y += 2) {
	    for (int x = 692; x <= 788; x += 2) {
		const SoBRLNavigationGizmo::Part candidate =
		    gizmo->hitTest(x, y, 800, 600);
		if (!SoBRLNavigationGizmo::partAet(candidate,
		    expectedAzimuth, expectedElevation))
		    continue;
		part = candidate;
		targetX = x;
		targetY = y;
		break;
	    }
	}
	CHECK(part != SoBRLNavigationGizmo::PART_NONE,
	    "the circles style must expose a canonical BRL-CAD axis target");

	BObolInputEvent press;
	press.type = BOBOL_INPUT_POINTER_PRESS;
	press.button = 0;
	press.buttons = 1u;
	press.x = targetX;
	press.y = targetY;
	CHECK(bobol_display_endpoint_input_dispatch(endpoint, &press) ==
	    BOBOL_INPUT_RESULT_HANDLED && fallbackPresses == 0,
	    "gizmo must consume only a press inside its footprint");
	BObolInputEvent release;
	release.type = BOBOL_INPUT_POINTER_RELEASE;
	release.button = 0;
	release.x = press.x;
	release.y = press.y;
	CHECK(bobol_display_endpoint_input_dispatch(endpoint, &release) ==
	    BOBOL_INPUT_RESULT_HANDLED,
	    "gizmo click release must be consumed");
	vect_t aet = VINIT_ZERO;
	CHECK(bv_aet_get(aet, bv_context_view(
	    reinterpret_cast<struct bv_context *>(viewContext))) &&
	    std::fabs(aet[0] - expectedAzimuth) < 1.0e-4 &&
	    std::fabs(aet[1] - expectedElevation) < 1.0e-4,
	    "gizmo click must use BRL-CAD AET rather than foreign Y-up orbit state");

	mat_t before;
	mat_t beforeRotation;
	mat_t beforeCenter;
	CHECK(bv_model2view_get(before, bv_context_view(
	    reinterpret_cast<struct bv_context *>(viewContext))),
	    "pre-drag BRL-CAD view matrix must be readable");
	CHECK(bv_rotation_get(beforeRotation, bv_context_view(
	    reinterpret_cast<struct bv_context *>(viewContext))) &&
	    bv_center_mat_get(beforeCenter, bv_context_view(
	    reinterpret_cast<struct bv_context *>(viewContext))),
	    "pre-drag BRL-CAD rotation and center must be readable");
	CHECK(bobol_display_endpoint_input_dispatch(endpoint, &press) ==
	    BOBOL_INPUT_RESULT_HANDLED,
	    "gizmo drag press must start inside the control");
	BObolInputEvent motion;
	motion.type = BOBOL_INPUT_POINTER_MOTION;
	motion.button = 0;
	motion.buttons = 1u;
	motion.x = targetX + 12;
	motion.y = targetY + 8;
	motion.dx = 12;
	motion.dy = 8;
	CHECK(bobol_display_endpoint_input_dispatch(endpoint, &motion) ==
	    BOBOL_INPUT_RESULT_HANDLED,
	    "gizmo drag must be consumed through the endpoint layer");
	release.x = motion.x;
	release.y = motion.y;
	CHECK(bobol_display_endpoint_input_dispatch(endpoint, &release) ==
	    BOBOL_INPUT_RESULT_HANDLED,
	    "gizmo drag release must complete the gesture");
	mat_t after;
	CHECK(bv_model2view_get(after, bv_context_view(
	    reinterpret_cast<struct bv_context *>(viewContext))),
	    "post-drag BRL-CAD view matrix must be readable");
	bool matrixChanged = false;
	for (int i = 0; i < 16; i++)
	    matrixChanged = matrixChanged ||
		std::fabs(before[i] - after[i]) > 1.0e-8;
	CHECK(matrixChanged,
	    "gizmo drag must use BRL-CAD's incremental view rotation path");

	struct bv *view = bv_context_view(
	    reinterpret_cast<struct bv_context *>(viewContext));
	CHECK(bv_rotation_set(view, beforeRotation) &&
	    bv_center_mat_set(view, beforeCenter) && bv_update(view),
	    "test must restore the pre-drag BRL-CAD view state");
	point_t expectedCenter = VINIT_ZERO;
	mat_t expected;
	CHECK(bv_center_get(expectedCenter, view) &&
	    bv_mouse_delta_adjust(view, motion.dx, motion.dy, expectedCenter,
		BV_ADJUST_ROT) && bv_model2view_get(expected, view),
	    "ordinary BRL-CAD mouse rotation must produce an expected matrix");
	bool exactOrbit = true;
	for (int i = 0; i < 16; i++)
	    exactOrbit = exactOrbit && std::fabs(after[i] - expected[i]) < 1.0e-8;
	CHECK(exactOrbit,
	    "gizmo drag must exactly match the normal MGED-style OrbitCamera path");
    }

    CHECK(run_command(gedp, "off") == BRLCAD_OK &&
	controller->features().exists("faceplate::navigation_gizmo",
	    BOBOL_FEATURE_SCOPE_LOCAL) && gizmo &&
	!gizmo->visible.getValue() && gizmo->getNumChildren() == 0,
	"disabling must hide the retained styled feature without HUD geometry");
    BObolInputEvent outsidePress;
    outsidePress.type = BOBOL_INPUT_POINTER_PRESS;
    outsidePress.button = 0;
    outsidePress.buttons = 1u;
    outsidePress.x = 740;
    outsidePress.y = 60;
    CHECK(bobol_display_endpoint_input_dispatch(endpoint, &outsidePress) ==
	BOBOL_INPUT_RESULT_HANDLED && fallbackPresses == 1,
	"hiding must clear only the gizmo layer and preserve other owners");

    CHECK(run_command(gedp, "on") == BRLCAD_OK,
	"gizmo should support a second enable cycle");
    CHECK(controller->features().node(handle) == node && gizmo &&
	gizmo->style.getValue() == SoBRLNavigationGizmo::CIRCLES,
	"enable must reuse the retained feature and preserve its chosen style");
    controller->features().clear();
    CHECK(bobol_display_endpoint_input_dispatch(endpoint, &outsidePress) ==
	BOBOL_INPUT_RESULT_HANDLED && fallbackPresses == 2,
	"view-store teardown must synchronously remove the gizmo input layer");

    (void)bobol_display_endpoint_input_action_layer_clear_if(endpoint,
	&fallbackOwner);
    ged_close(gedp);
    bu_vls_free(&path);
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
