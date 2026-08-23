/*              G E D _ B O B O L _ P R I V A T E . H P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_bobol_private.hpp
 *
 * Private C++ integration between drawing-aware libged commands and libBObol.
 * This is deliberately not an installed libged facade: commands may include
 * libBObol headers and use the borrowed scene and controller objects directly.
 * The only state packaged here is the explicit, stack-owned context needed
 * while a database walk publishes several source occurrences as one batch.
 */

#ifndef LIBGED_GED_BOBOL_PRIVATE_HPP
#define LIBGED_GED_BOBOL_PRIVATE_HPP

#include <cstdio>
#include <string>

#include "BObol/BViewStore.h"
#include "BObol/BViewController.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BFramebuffer.h"
#include "bu/str.h"
#include "bv.h"
#include "./ged_draw_types_private.h"
#include "ged/view.h"

struct db_i;
struct ged;
struct ged_view_context;
class BObolViewController;
class BObolFramebufferStream;
class BObolSceneController;

struct ged_bobol_publication_context {
    struct ged *gedp = nullptr;
    struct ged_view_context *view_ctx = nullptr;
    BObolSceneController *scene = nullptr;
    int draw_mode = -1;
    std::string group_path;
    const struct ged_draw_appearance_settings *appearance = nullptr;
    bool active = false;
};

/* Borrowed accessors.  A view controller is present exactly when its display
 * endpoint carries one; there is no primary-controller fallback. */
BObolSceneController *ged_bobol_scene(struct ged *gedp);
BObolViewController *ged_bobol_view_controller(
    struct ged_view_context *view_ctx);
BObolViewController *ged_bobol_shared_view_controller(struct ged *gedp);

typedef int (*ged_bobol_view_controller_visit_t)(
    struct ged_view_context *view_ctx,
    BObolViewController *controller,
    void *userdata);
void ged_bobol_view_controllers_foreach(
    struct ged *gedp,
    ged_bobol_view_controller_visit_t callback,
    void *userdata);

/* A shared feature belongs to the GED-wide scene root, while each graphical
 * view owns its own retained presentation.  Wake every view that renders the
 * shared root after an already-published shared-content mutation. */
void ged_bobol_shared_view_presentation_request(struct ged *gedp,
    const char *reason);

/* Source inventory is GED-wide; instance identities say whether an entry is
 * shared or belongs to an independent view.  Keep that value interpretation
 * beside the borrowed accessors so command queries do not depend on a render
 * cache having been populated. */
inline bool
ged_bobol_source_in_view(struct ged_view_context *view_ctx,
    const BObolDatabaseSourceSummary &summary)
{
    if (!summary.valid || summary.path.getLength() == 0)
	return false;

    std::string key = summary.instanceKey.getString();
    const std::string mode_marker(":ged-draw-mode:");
    const size_t marker = key.rfind(mode_marker);
    if (marker != std::string::npos)
	key.erase(marker);

    if (!view_ctx || !ged_view_context_is_independent(view_ctx)) {
	if (key.empty())
	    return true;
	if (key.compare(0, 9, "ged-view:") == 0)
	    return false;
	if (key.compare(0, 14, "brlcad-direct:") == 0)
	    return true;
	const char *source_path = summary.path.getString();
	while (*source_path == '/')
	    source_path++;
	const char *instance_path = key.c_str();
	while (*instance_path == '/')
	    instance_path++;
	return BU_STR_EQUAL(instance_path, source_path);
    }

    const char *view_name = bv_name_get(
	bv_context_view_const(reinterpret_cast<const struct bv_context *>(view_ctx)));
    char fallback[64] = {0};
    if (!view_name || !view_name[0]) {
	snprintf(fallback, sizeof(fallback), "%p", static_cast<void *>(view_ctx));
	view_name = fallback;
    }
    std::string prefix("ged-view:");
    prefix += view_name;
    prefix += ":";
    return key.compare(0, prefix.size(), prefix) == 0;
}

/* Execute while the Obol framebuffer transport lock is held.  The stream is
 * borrowed only for the callback; command code performs the actual libBObol
 * operation directly. */
typedef int (*ged_bobol_framebuffer_operation_t)(
    BObolFramebufferStream &stream, void *userdata);
int ged_bobol_framebuffer_apply(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    ged_bobol_framebuffer_operation_t operation,
    void *userdata,
    bool publish);

/* Small value-only helpers for direct feature-store use.  These do not own,
 * attach, or preserve any controller state. */
struct ged_bobol_feature_binding {
    BObolViewController *controller = nullptr;
    BObolFeatureHandle handle;
};

struct ged_bobol_polygon_binding {
    BObolViewController *controller = nullptr;
    BObolPolygonHandle handle;
};

inline BObolFeatureOwner
ged_bobol_view_feature_owner(struct ged_view_context *view_ctx)
{
    BObolFeatureOwner owner;
    owner.ownerToken = view_ctx;
    char owner_id[64] = {0};
    snprintf(owner_id, sizeof(owner_id), "%p", static_cast<void *>(view_ctx));
    owner.ownerId = owner_id;
    owner.ownerRole = "view";
    return owner;
}

inline BObolViewController *
ged_bobol_feature_controller(struct ged *gedp,
    struct ged_view_context *view_ctx, int local)
{
    return local ? ged_bobol_view_controller(view_ctx) :
	ged_bobol_shared_view_controller(gedp);
}

inline struct ged_bobol_feature_binding
ged_bobol_feature_find(struct ged *gedp, struct ged_view_context *view_ctx,
    const char *name)
{
    struct ged_bobol_feature_binding binding;
    if (!name || !name[0])
	return binding;

    BObolViewController *local = ged_bobol_view_controller(view_ctx);
    if (local) {
	const BObolFeatureOwner owner =
	    ged_bobol_view_feature_owner(view_ctx);
	binding.handle = local->features().findOwned(name,
	    BOBOL_FEATURE_SCOPE_LOCAL, &owner);
	if (!binding.handle.isValid())
	    binding.handle = local->features().find(name,
		BOBOL_FEATURE_SCOPE_SHARED);
	if (binding.handle.isValid()) {
	    binding.controller = local;
	    return binding;
	}
    }

    BObolViewController *shared = ged_bobol_shared_view_controller(gedp);
    if (shared && shared != local) {
	binding.handle = shared->features().find(name,
	    BOBOL_FEATURE_SCOPE_SHARED);
	if (binding.handle.isValid())
	    binding.controller = shared;
    }
    return binding;
}

inline struct ged_bobol_polygon_binding
ged_bobol_polygon_find(struct ged *gedp, struct ged_view_context *view_ctx,
    const char *name, bool local_only = false)
{
    struct ged_bobol_polygon_binding binding;
    if (!gedp || !view_ctx || !name || !name[0])
	return binding;

    BObolViewController *local = ged_bobol_view_controller(view_ctx);
    if (local) {
	binding.handle = local->polygons().find(name,
	    BOBOL_FEATURE_SCOPE_LOCAL);
	if (binding.handle.isValid()) {
	    binding.controller = local;
	    return binding;
	}
    }
    if (local_only)
	return binding;

    BObolViewController *shared = ged_bobol_shared_view_controller(gedp);
    if (shared) {
	binding.handle = shared->polygons().find(name,
	    BOBOL_FEATURE_SCOPE_SHARED);
	if (binding.handle.isValid())
	    binding.controller = shared;
    }
    return binding;
}

bool ged_bobol_publication_begin(
    struct ged_bobol_publication_context *publication,
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    int draw_mode);
void ged_bobol_publication_appearance_set(
    struct ged_bobol_publication_context *publication,
    const struct ged_draw_appearance_settings *settings);
void ged_bobol_publication_group_set(
    struct ged_bobol_publication_context *publication,
    const char *group_path);
void ged_bobol_publication_end(
    struct ged_bobol_publication_context *publication);

int ged_bobol_database_source_ensure_for_path_with_placement(
    struct ged_bobol_publication_context *publication,
    const char *path,
    struct db_i *dbip,
    int draw_mode,
    uint64_t source_revision,
    int draw_mat_valid,
    const mat_t draw_mat,
    int draw_center_valid,
    const point_t draw_center,
    int draw_size_valid,
    fastf_t draw_size);

int ged_bobol_database_source_ensure_for_path(
    struct ged_bobol_publication_context *publication,
    const char *path,
    struct db_i *dbip,
    int draw_mode,
    uint64_t source_revision);

int ged_bobol_database_source_publish_indexed_face_set_for_path(
    struct ged_bobol_publication_context *publication,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count);

int ged_bobol_database_source_update_display_for_path(
    struct ged_bobol_publication_context *publication,
    const char *path,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision);

#endif /* LIBGED_GED_BOBOL_PRIVATE_HPP */
