/*                  D R A W _ O B O L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file draw_obol.cpp
 *
 * libged bridge from neutral GED draw transactions to Obol scene state.
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BADC.h"
#include "BObol/BDrawCache.h"
#include "BObol/BGrid.h"
#include "BObol/BInit.h"
#include "BObol/BImageSource.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BLodRealization.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSceneController.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BSnapAction.h"
#include "BObol/BSourceRealization.h"
#include "BObol/BViewportImage.h"
#include "BObol/BVListShape.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "BObol/BViewQuery.h"
#include "BObol/BViewStore.h"
#include "bg/line_layer.h"
#include "bg/plane.h"
#include "bg/polygon.h"
#include "bv.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/path.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "bu/units.h"
#include "bu/vls.h"
#include "ged.h"
#include "ged/db_index.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "ged/selection.h"
#include "ged/view.h"
#include "icv.h"
#include "rt/db5.h"
#include "rt/db_fullpath.h"
#include "rt/search.h"
#include "rt/view.h"
#include "vmath.h"

#include "./ged_bobol_private.hpp"
#include "./draw_obol_bridge_private.hpp"
#include "./ged_draw_view_private.h"
#include "./ged_private.h"

#include <algorithm>
#include <Inventor/SbVec2f.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Obol/cad/SoCADAssembly.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <atomic>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>


extern "C" int
ged_draw_obol_progressive_available(struct ged *gedp, struct ged_view_context *view_ctx)
{
    (void)gedp;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    return (controller && controller->findProgressiveProviderData(
	ged_obol_progressive_advance_provider)) ? 1 : 0;
}

/* Progressive (coarse-first / deferred-leaf) drawing is driven by the part of
 * the view policy used by the requested representation.  Shaded and
 * hidden-line modes consume mesh LoD.  Wire may consume CSG LoD for primitive
 * plotting and mesh LoD for PoP-backed BoTs, so either policy can make its
 * production path progressive. */
extern "C" int
ged_draw_obol_view_lod_enabled(struct ged *gedp,
			      struct ged_view_context *view_ctx,
			      int draw_mode)
{
    (void)gedp;
    if (!ged_bobol_view_controller(view_ctx))
	return 0;

    ged_view_lod_policy policy;
    if (!ged_view_lod_policy_get(&policy, view_ctx) ||
	policy.policy == BV_LOD_OFF)
	return 0;

    switch (draw_mode) {
	case GED_DRAW_MODE_SHADED_BOTS:
	case GED_DRAW_MODE_SHADED:
	case GED_DRAW_MODE_HIDDEN_LINE:
	    return policy.mesh_enabled ? 1 : 0;
	case GED_DRAW_MODE_WIRE:
	    return (policy.csg_enabled || policy.mesh_enabled) ? 1 : 0;
	case GED_DRAW_MODE_EVAL_WIRE:
	case GED_DRAW_MODE_EVAL_POINTS:
	default:
	    return 0;
    }
}

int
ged_obol_bind_view_render_root(struct ged_view_context *view_ctx,
			       BObolSceneController *shared_scene,
			       BObolViewController *view_controller)
{
    if (!shared_scene || !view_controller)
	return 0;

    SoNode *shared_root = shared_scene->getSceneRoot();
    if (!shared_root)
	return 0;

    view_controller->getSceneController()->shareRealizationRepository(
	shared_scene);

    SoNode *local_root = view_controller->getSceneRoot();
    if (!local_root || local_root == shared_root) {
	SoBRLSceneGroup *view_group = new SoBRLSceneGroup;
	std::string group_path("_view/");
	group_path += ged_obol_view_scope_name(view_ctx);
	view_group->groupPath = group_path.c_str();
	view_controller->setSceneRoot(view_group, TRUE);
	local_root = view_group;
    }

    SoSeparator *render_root = new SoSeparator;
    if (!ged_obol_view_scope_is_independent(view_ctx))
	render_root->addChild(shared_root);
    /* The interlay belongs after model geometry but before the view-local
     * faceplate and interactive feature root. */
    SoGroup *interlay = view_controller->getFramebufferInterlayRoot();
    if (!interlay)
	return 0;
    render_root->addChild(interlay);
    if (local_root && local_root != shared_root)
	render_root->addChild(local_root);
    /* The shared/interlay/local wrapper changes presentation ownership, not
     * the modeled source.  Retain an in-flight cold PoP build across this
     * endpoint composition install; source replacement still uses the normal
     * destructive render-root API. */
    view_controller->setRenderSceneRoot(render_root, TRUE);
    return 1;
}

int
ged_obol_node_contains(SoNode *root, SoNode *target)
{
    if (!root || !target)
	return 0;
    if (root == target)
	return 1;
    if (!root->isOfType(SoGroup::getClassTypeId()))
	return 0;

    SoGroup *group = static_cast<SoGroup *>(root);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (ged_obol_node_contains(group->getChild(i), target))
	    return 1;
    }
    return 0;
}

static int
ged_draw_obol_attach_view_common(struct ged *gedp,
				 struct ged_view_context *view_ctx,
				 BObolSceneController *scene_controller,
				 BObolViewController *view_controller,
				 int sync_current_scene)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp || !scene_controller || !view_controller)
	return 0;

    ged_obol_configure_compact_realization(gedp, view_ctx,
	    scene_controller);

    if (!ged_obol_observer_ensure(gedp, gdp))
	return 0;

    if (!ged_draw_obol_scene_controller_ensure_owned(gedp,
	    sync_current_scene))
	return 0;
    BObolSceneController *shared_scene =
	ged_draw_obol_scene_controller(gedp);
    if (!shared_scene)
	return 0;

    (void)ged_view_context_host_attach(gedp, view_ctx);
    if (!ged_view_context_obol_attachment_bind(view_ctx,
	    view_controller->getViewAttachment()))
	return 0;

    /* Keep the renderer-neutral view record synchronized with the retained
     * root it now renders.  This is an identity attachment only: all
     * non-independent views continue to share the same source nodes and
     * geometry repository. */
    (void)ged_draw_source_root_attach_view_contexts(gedp, view_ctx, NULL);

    /* The shared controller may itself be a single-view endpoint.  Its scene
     * root already is the shared composition, so replacing it with a local
     * render composition would make the shared scene self-referential. */
    if (view_controller != ged_draw_obol_controller(gedp)) {
	void *existing = view_controller->findProgressiveProviderData(
	    ged_obol_progressive_advance_provider);
	SoNode *shared_root = shared_scene->getSceneRoot();
	const int shared_visible =
	    !ged_obol_view_scope_is_independent(view_ctx);
	const int shared_bound = ged_obol_node_contains(
	    view_controller->getRenderSceneRoot(), shared_root);
	if (!existing || shared_visible != shared_bound) {
	    if (!ged_obol_bind_view_render_root(view_ctx, shared_scene,
		    view_controller)) {
		(void)ged_view_context_obol_attachment_unbind(view_ctx,
		    view_controller->getViewAttachment());
		return 0;
	    }
	}

	if (sync_current_scene &&
		ged_obol_view_scope_is_independent(view_ctx))
	    (void)ged_draw_obol_scene_sync_full_scene(gedp, view_ctx,
		ged_obol_fold_revision(ged_draw_scene_revision(gedp)),
		view_controller->getSceneController());
    }

    if (!ged_obol_register_progressive_provider(gedp, view_ctx,
	    view_controller)) {
	(void)ged_view_context_obol_attachment_unbind(view_ctx,
	    view_controller->getViewAttachment());
	return 0;
    }

    /* Every graphical endpoint uses the same automatic display-LoD service.
     * Starting it only when cold structural bounds happened to be missing made
     * cold cache misses work while warm caches (and small models with directly
     * available bounds) never traversed their valid source-mesh requests.
     * Failure is non-fatal: the standing compact boxes remain a useful
     * fallback and an explicit service start may be retried later. */
    (void)ged_obol_lod_service_ensure(view_controller);

    return 1;
}

extern "C" int
ged_draw_obol_controller_attach_opaque_for_view(struct ged *gedp,
	struct ged_view_context *view_ctx,
	void *controller,
	int sync_current_scene)
{
    BObolViewController *view_controller =
	static_cast<BObolViewController *>(controller);
    if (!view_controller)
	return 0;
    return ged_draw_obol_attach_view_common(gedp, view_ctx,
	view_controller->getSceneController(), view_controller,
	sync_current_scene);
}

extern "C" int
ged_draw_obol_render_endpoint_ensure_for_view(struct ged *gedp,
	struct ged_view_context *view_ctx,
	int sync_current_scene)
{
    (void)sync_current_scene;
    if (!gedp || !view_ctx)
	return 0;
    if (ged_view_context_obol_endpoint_get(view_ctx))
	return 1;

    bobol_init(NULL);
    BObolViewController *controller =
	new (std::nothrow) BObolViewController(new SoBRLSceneGroup);
    if (!controller)
	return 0;

    bobol_display_endpoint_t *endpoint =
	bobol_display_endpoint_create(controller,
	    BOBOL_ENDPOINT_OWN_CONTROLLER);
    if (!endpoint) {
	delete controller;
	return 0;
    }
    if (!ged_view_context_obol_endpoint_set(view_ctx, endpoint, 1)) {
	bobol_display_endpoint_destroy(endpoint);
	return 0;
    }
    return 1;
}

void
ged_draw_obol_scene_controller_detach(struct ged *gedp)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp)
	return;

    ged_obol_state_destroy(gdp);
}

extern "C" void
ged_draw_obol_controller_detach_for_view(struct ged *gedp,
					 struct ged_view_context *view_ctx)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp || !view_ctx)
	return;

    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	return;
    const uint64_t token = controller->findProgressiveProviderToken(
	ged_obol_progressive_advance_provider);
    if (token)
	controller->unregisterProgressiveProvider(token);
    controller->stopManagedLodService();
    (void)ged_view_context_obol_attachment_unbind(view_ctx,
	controller->getViewAttachment());
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
