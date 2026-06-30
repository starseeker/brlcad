/*                  T E S T _ S C E N E _ R E F . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file test_scene_ref.cpp
 *
 * Tests Obol-tagged neutral scene refs exported by SoBRLSceneController.
 */

#include "common.h"

#include "brlobol/init.h"
#include "brlobol/scene_controller.h"
#include "bu/app.h"
#include "rt/view.h"

#include <Inventor/nodes/SoSeparator.h>

#include <stdio.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

int
main(int argc, char **argv)
{
    (void)argc;
    bu_setprogname(argv[0]);
    brlobol_init(NULL);

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLSceneController controller(root);
    const SoBRLSceneController &constController = controller;

    rt_view_scene_ref ref = controller.getSceneRef();
    rt_view_scene_ref constRef = constController.getSceneRef();
    if (rt_view_scene_ref_is_null(ref))
	FAIL("Obol scene-controller ref should be non-null");
    if (rt_view_scene_ref_backend(ref) != RT_VIEW_SCENE_BACKEND_OBOL)
	FAIL("Obol scene-controller ref should carry the Obol backend tag");
    if (rt_view_scene_ref_context(ref) != &controller)
	FAIL("Obol scene-controller ref should carry the controller pointer");
    if (!rt_view_scene_ref_equal(ref, constRef))
	FAIL("const and mutable Obol scene-controller refs should compare equal");
    if (!SoBRLSceneController::sceneRefIsObol(ref))
	FAIL("sceneRefIsObol should accept controller refs");
    if (SoBRLSceneController::fromSceneRef(ref) != &controller)
	FAIL("fromSceneRef should recover the controller from an Obol ref");
    if (SoBRLSceneController::fromConstSceneRef(ref) != &controller)
	FAIL("fromConstSceneRef should recover the controller from an Obol ref");

    rt_view_scene_ref nullRef = rt_view_scene_ref_null();
    if (SoBRLSceneController::sceneRefIsObol(nullRef))
	FAIL("sceneRefIsObol should reject null refs");
    if (SoBRLSceneController::fromSceneRef(nullRef))
	FAIL("fromSceneRef should reject null refs");

    rt_view_scene_ref bsgTagged =
	rt_view_scene_ref_make(&controller, RT_VIEW_SCENE_BACKEND_BSG);
    if (SoBRLSceneController::sceneRefIsObol(bsgTagged))
	FAIL("sceneRefIsObol should reject BSG-tagged refs");
    if (SoBRLSceneController::fromSceneRef(bsgTagged))
	FAIL("fromSceneRef should reject BSG-tagged refs");
    if (rt_view_scene_ref_equal(ref, bsgTagged))
	FAIL("neutral scene-ref equality should include the backend tag");

    rt_view_scene_ref legacyUntagged =
	rt_view_scene_ref_make(&controller, RT_VIEW_SCENE_BACKEND_NONE);
    if (SoBRLSceneController::sceneRefIsObol(legacyUntagged))
	FAIL("sceneRefIsObol should reject opaque legacy refs");
    if (SoBRLSceneController::fromSceneRef(legacyUntagged))
	FAIL("fromSceneRef should reject opaque legacy refs");
    if (!rt_view_scene_ref_equal(ref, legacyUntagged))
	FAIL("legacy opaque scene refs should still compare by pointer");

    root->unref();
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
