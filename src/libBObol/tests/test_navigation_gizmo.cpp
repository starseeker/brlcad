/*             T E S T _ N A V I G A T I O N _ G I Z M O . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BInit.h"
#include "BObol/BNavigationGizmo.h"

#include <Inventor/SbBox3f.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/nodes/SoOrthographicCamera.h>

#include <cmath>
#include <cstdio>

#define CHECK(_condition, _message) do { \
    if (!(_condition)) { \
	std::fprintf(stderr, "FAIL: %s\n", (_message)); \
	failures++; \
    } \
} while (0)

static bool
close_value(float a, float b)
{
    return std::fabs(a - b) < 1.0e-4f;
}

int
main(void)
{
    int failures = 0;
    bobol_init(NULL);

    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->ref();
    SoBRLNavigationGizmo *gizmo = new SoBRLNavigationGizmo;
    gizmo->ref();
    gizmo->setCamera(camera);
    CHECK(gizmo->getCamera() == camera,
	"gizmo should retain and track its assigned camera");
    CHECK(gizmo->rebuildGeometry() != NULL && gizmo->getNumChildren() == 1,
	"visible gizmo should build one retained HUD subtree");
    CHECK(gizmo->style.getValue() == SoBRLNavigationGizmo::CUBE,
	"navigation gizmo should default to the cube presentation");

    SbVec3f direction;
    CHECK(SoBRLNavigationGizmo::partDirection(
	SoBRLNavigationGizmo::PART_POS_X_POS_Y_NEG_Z, direction) &&
	close_value(direction[0], 1.0f) && close_value(direction[1], 1.0f) &&
	close_value(direction[2], -1.0f),
	"part direction encoding should preserve all three BRL-CAD axes");
    CHECK(!SoBRLNavigationGizmo::partDirection(
	SoBRLNavigationGizmo::PART_NONE, direction),
	"the empty part must not decode to an orientation");
    CHECK(!SoBRLNavigationGizmo::partDirection(
	SoBRLNavigationGizmo::PART_ORBIT, direction),
	"the free-orbit field must not decode to a snap orientation");

    float azimuth = 0.0f;
    float elevation = 0.0f;
    CHECK(SoBRLNavigationGizmo::partAet(
	SoBRLNavigationGizmo::PART_POS_X, azimuth, elevation) &&
	close_value(azimuth, 0.0f) && close_value(elevation, 0.0f),
	"+X must use BRL-CAD's canonical front AET");
    CHECK(SoBRLNavigationGizmo::partAet(
	SoBRLNavigationGizmo::PART_POS_Y, azimuth, elevation) &&
	close_value(azimuth, 90.0f) && close_value(elevation, 0.0f),
	"+Y must use BRL-CAD's canonical left AET");
    CHECK(SoBRLNavigationGizmo::partAet(
	SoBRLNavigationGizmo::PART_NEG_Y, azimuth, elevation) &&
	close_value(azimuth, 270.0f) && close_value(elevation, 0.0f),
	"-Y must use BRL-CAD's canonical right AET");
    CHECK(SoBRLNavigationGizmo::partAet(
	SoBRLNavigationGizmo::PART_POS_Z, azimuth, elevation) &&
	close_value(azimuth, 270.0f) && close_value(elevation, 90.0f),
	"+Z must use BRL-CAD's canonical top AET pole convention");

    /* Default upper-right placement: center is (740, 60) in normalized input
     * coordinates for an 800x600 viewport.  An identity camera exposes +Z. */
    CHECK(gizmo->hitTest(740, 60, 800, 600) ==
	SoBRLNavigationGizmo::PART_POS_Z,
	"center hit should resolve the visible face");
    CHECK(gizmo->hitTest(760, 60, 800, 600) ==
	SoBRLNavigationGizmo::PART_POS_X_POS_Z,
	"face boundary hit should resolve the visible edge");
    CHECK(gizmo->hitTest(760, 40, 800, 600) ==
	SoBRLNavigationGizmo::PART_POS_X_POS_Y_POS_Z,
	"two face-boundary coordinates should resolve the visible corner");
    CHECK(gizmo->hitTest(100, 100, 800, 600) ==
	SoBRLNavigationGizmo::PART_NONE,
	"pixels outside the control must fall through to other input owners");

    gizmo->style = SoBRLNavigationGizmo::CIRCLES;
    CHECK(gizmo->rebuildGeometry() != NULL && gizmo->getNumChildren() == 1,
	"circle style should build one retained HUD subtree");
    CHECK(gizmo->hitTest(781, 60, 800, 600) ==
	SoBRLNavigationGizmo::PART_POS_X,
	"positive circle-axis endpoint should expose the BRL-CAD +X target");
    CHECK(gizmo->hitTest(699, 60, 800, 600) ==
	SoBRLNavigationGizmo::PART_NEG_X,
	"negative circle-axis endpoint should expose the BRL-CAD -X target");
    CHECK(gizmo->hitTest(760, 50, 800, 600) ==
	SoBRLNavigationGizmo::PART_ORBIT,
	"the circle field should expose a draggable free-orbit target");
    CHECK(gizmo->hitTest(100, 100, 800, 600) ==
	SoBRLNavigationGizmo::PART_NONE,
	"circle style must preserve input fallthrough outside its footprint");

    gizmo->setHoverPart(SoBRLNavigationGizmo::PART_POS_Z);
    gizmo->setActivePart(SoBRLNavigationGizmo::PART_POS_X_POS_Z);
    CHECK(gizmo->hoverPart.getValue() == SoBRLNavigationGizmo::PART_POS_Z &&
	gizmo->activePart.getValue() ==
	SoBRLNavigationGizmo::PART_POS_X_POS_Z,
	"retained hover and active fields should follow interaction state");

    SoGetBoundingBoxAction bounds(SbViewportRegion(800, 600));
    bounds.apply(gizmo);
    CHECK(bounds.getBoundingBox().isEmpty(),
	"screen-locked gizmo must not alter CAD autoview extents");

    gizmo->visible = FALSE;
    CHECK(gizmo->rebuildGeometry() == NULL && gizmo->getNumChildren() == 0 &&
	gizmo->hitTest(740, 60, 800, 600) ==
	    SoBRLNavigationGizmo::PART_NONE,
	"hidden gizmo should retain no HUD geometry and accept no picks");

    gizmo->unref();
    camera->unref();
    return failures ? 1 : 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
