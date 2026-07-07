/*                    D R A W _ O B O L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file ged/draw_obol.h
 *
 * C++ bridge between libged's neutral draw transaction API and an
 * Obol/libbrlobol view controller.
 */

#ifndef GED_DRAW_OBOL_H
#define GED_DRAW_OBOL_H

#include "common.h"

#include <stddef.h>
#include <stdint.h>

#include "ged/defines.h"
#include "vmath.h"

__BEGIN_DECLS

struct rt_db_internal;
struct bg_tess_tol;
struct bn_tol;

struct ged_draw_obol_lod_service_status {
    int attached;
    int running;
    int auto_submit;
    size_t worker_count;
    size_t in_flight;
    size_t pending_tasks;
    size_t queued_results;
    size_t queued_cache_writes;
    size_t delayed_tasks;
    unsigned int last_visited_mesh_count;
    unsigned int last_submitted_task_count;
    unsigned int last_skipped_mesh_count;
    size_t last_result_count;
    unsigned int last_matched_result_count;
    unsigned int last_applied_result_count;
    unsigned int last_rejected_result_count;
    unsigned int last_unmatched_result_count;
    size_t active_mesh_payloads;
    size_t active_aabb_proxy_payloads;
    size_t active_obb_proxy_payloads;
    char last_diagnostics[2048];
};

struct ged_draw_obol_source_expansion_status {
    size_t child_count;
    size_t considered;
    size_t expanded;
    size_t existing;
    size_t skipped_non_union;
    size_t skipped_duplicate_instance;
    size_t skipped_invalid;
    size_t remaining;
    size_t proxy_published;
    size_t metadata_applied;
    size_t comb_sources;
    size_t leaf_sources;
};

struct ged_draw_obol_source_prewarm_status {
    size_t child_count;
    size_t considered;
    size_t submitted;
    size_t already_cached;
    size_t skipped_non_union;
    size_t skipped_duplicate_instance;
    size_t skipped_invalid;
    size_t remaining;
    size_t comb_sources;
    size_t leaf_sources;
};

/**
 * C-compatible form of ged_draw_obol_controller_attach.
 *
 * @p controller is a borrowed BRLObolViewController pointer passed as opaque
 * storage so C callers such as the GED dm command can attach a libdm-owned
 * Obol backend without exposing C++ types.
 */
GED_EXPORT int
ged_draw_obol_controller_attach_opaque(struct ged *gedp,
				       void *controller,
				       int sync_current_scene);

/**
 * Attach a borrowed Obol view controller for one GED view.
 *
 * Unlike ged_draw_obol_controller_attach_opaque, this does not replace other
 * attached Obol controllers.  It is intended for libdm-owned Obol display
 * managers, where each view has its own renderer/controller.
 */
GED_EXPORT int
ged_draw_obol_controller_attach_opaque_for_view(struct ged *gedp,
	void *view_ctx,
	void *controller,
	int sync_current_scene);

/**
 * Detach a previously borrowed opaque Obol view controller.
 */
GED_EXPORT void
ged_draw_obol_controller_detach_opaque(struct ged *gedp,
				       void *controller);

/**
 * Present pending Obol framebuffer stream state into the attached view scene.
 *
 * Applications should call this from their view/update thread before render
 * traversal.  fbserv producer threads only update image-stream state; they do
 * not mutate Coin/Obol scene graph fields directly.
 */
GED_EXPORT int
ged_draw_obol_framebuffer_present(struct ged *gedp);

/**
 * Release the libged-owned Obol framebuffer bridge, if one is active.
 *
 * Call this before destroying an attached Obol display manager/controller.
 */
GED_EXPORT void
ged_draw_obol_framebuffer_release(struct ged *gedp);

GED_EXPORT int
ged_draw_obol_lod_service_start(struct ged *gedp,
				void *view_ctx,
				size_t worker_count);

GED_EXPORT int
ged_draw_obol_lod_service_stop(struct ged *gedp,
			       void *view_ctx);

GED_EXPORT int
ged_draw_obol_lod_service_poll(struct ged *gedp,
			       void *view_ctx,
			       size_t max_results,
			       struct ged_draw_obol_lod_service_status *status);

GED_EXPORT int
ged_draw_obol_lod_service_status(struct ged *gedp,
				 void *view_ctx,
				 struct ged_draw_obol_lod_service_status *status);

GED_EXPORT int
ged_draw_obol_view_lod_policy_changed(struct ged *gedp,
				      void *view_ctx);

GED_EXPORT size_t
ged_draw_obol_lod_service_prewarm(struct ged *gedp,
				  void *view_ctx,
				  int argc,
				  const char * const *argv,
				  struct ged_draw_obol_lod_service_status *status);

GED_EXPORT int
ged_draw_obol_database_source_expand_children(
    struct ged *gedp,
    void *view_ctx,
    const char *path,
    int ged_draw_mode,
    size_t max_children,
    struct ged_draw_obol_source_expansion_status *status);

GED_EXPORT size_t
ged_draw_obol_database_source_prewarm_child_aabb_proxies(
    struct ged *gedp,
    void *view_ctx,
    const char *path,
    int ged_draw_mode,
    size_t max_children,
    struct ged_draw_obol_source_prewarm_status *status);

GED_EXPORT size_t
ged_draw_obol_database_source_prewarm_visible_child_aabb_proxies(
    struct ged *gedp,
    void *view_ctx,
    const char *root_path,
    int ged_draw_mode,
    size_t max_sources,
    size_t max_children_per_source,
    struct ged_draw_obol_source_prewarm_status *status);

GED_EXPORT int
ged_draw_obol_database_source_expand_visible_children(
    struct ged *gedp,
    void *view_ctx,
    const char *root_path,
    int ged_draw_mode,
    size_t max_sources,
    size_t max_children_per_source,
    struct ged_draw_obol_source_expansion_status *status);

GED_EXPORT int
ged_draw_obol_database_source_count(struct ged *gedp,
				    int skip_overlay_groups,
				    size_t *out);

GED_EXPORT int
ged_draw_obol_database_source_remove_for_path(struct ged *gedp,
					      const char *path);

GED_EXPORT int
ged_draw_obol_scene_controller_ensure_owned(struct ged *gedp,
					    int sync_current_scene);

GED_EXPORT int
ged_draw_view_context_obol_scene_adapter_attach(struct ged *gedp,
						void *view_ctx);

GED_EXPORT int
ged_draw_view_context_obol_feature_adapter_attach(struct ged *gedp,
						  void *view_ctx);

GED_EXPORT int
ged_draw_view_context_obol_polygon_adapter_attach(struct ged *gedp,
						  void *view_ctx);

GED_EXPORT int
ged_draw_view_context_obol_selection_adapter_attach(struct ged *gedp,
						    void *view_ctx);

GED_EXPORT size_t
ged_draw_obol_view_context_clear(void *view_ctx,
				 int flags);

GED_EXPORT int
ged_draw_obol_database_source_publish_line_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    const int *commands,
    size_t point_count);

GED_EXPORT int
ged_draw_obol_database_source_publish_point_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count);

GED_EXPORT int
ged_draw_obol_database_source_publish_indexed_face_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count);

GED_EXPORT int
ged_draw_obol_database_source_publish_lod_indexed_face_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count);

GED_EXPORT int
ged_draw_obol_database_source_publish_primitive_wireframe_for_path(
    struct ged *gedp,
    const char *path,
    struct rt_db_internal *ip,
    const struct bg_tess_tol *ttol,
    const struct bn_tol *tol);

__END_DECLS

#ifdef __cplusplus

class BRLObolViewController;
class SoBRLSceneController;
struct ged_draw_transaction;
struct ged_draw_transaction_result;

/**
 * Attach a borrowed Obol scene controller to @p gedp's draw state.
 *
 * The controller remains owned by the caller.  While attached, successful GED
 * draw transactions are mirrored into the controller's database-source scene.
 * When @p sync_current_scene is non-zero, the current retained GED draw set is
 * also synchronized immediately.
 *
 * This is the preferred Capability 1 boundary: libged mutates the
 * libbrlobol/Obol scene model without depending on GUI view policy.
 *
 * Returns 1 on successful attachment, 0 on invalid input or observer
 * registration failure.
 */
GED_EXPORT int
ged_draw_obol_scene_controller_attach(struct ged *gedp,
				      SoBRLSceneController *controller,
				      int sync_current_scene = 1);

/**
 * Stop mirroring GED draw transactions to the currently attached Obol scene
 * controller.
 */
GED_EXPORT void
ged_draw_obol_scene_controller_detach(struct ged *gedp);

/**
 * Return the borrowed scene controller currently attached to @p gedp, or NULL.
 */
GED_EXPORT SoBRLSceneController *
ged_draw_obol_scene_controller(struct ged *gedp);

/**
 * Return the active Obol scene controller, creating a libged-owned controller
 * when none is attached.
 *
 * The owned controller is synchronized from the current GED draw set when
 * @p sync_current_scene is non-zero, receives subsequent draw transactions
 * through the same observer path as borrowed controllers, and is destroyed
 * automatically by ged_draw_obol_scene_controller_detach or GED teardown.
 */
GED_EXPORT SoBRLSceneController *
ged_draw_obol_scene_controller_ensure(struct ged *gedp,
				      int sync_current_scene = 1);

/**
 * Return non-zero when @p gedp's active Obol scene controller is owned by
 * libged rather than borrowed from qtcad/qged or another caller.
 */
GED_EXPORT int
ged_draw_obol_scene_controller_owned(struct ged *gedp);

/**
 * Mirror one completed GED draw transaction into @p controller.
 *
 * If @p controller is NULL, the currently attached scene controller is used.
 * Returns 1 if the Obol scene changed, 0 otherwise.
 */
GED_EXPORT int
ged_draw_obol_scene_sync_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    SoBRLSceneController *controller = NULL);

/**
 * Rebuild @p controller's database-source scene from the current GED draw set.
 *
 * If @p controller is NULL, the currently attached scene controller is used.
 * A zero @p source_revision requests a folded value derived from GED's current
 * draw scene revision.  Returns 1 if the Obol scene changed, 0 otherwise.
 */
GED_EXPORT int
ged_draw_obol_scene_sync_full_scene(struct ged *gedp,
				    void *view_ctx,
				    uint32_t source_revision = 0,
				    SoBRLSceneController *controller = NULL);

/**
 * Attach a borrowed Obol view controller to @p gedp's draw state.
 *
 * This compatibility wrapper attaches the view controller's scene controller.
 * New libged integrations should prefer ged_draw_obol_scene_controller_attach.
 */
GED_EXPORT int
ged_draw_obol_controller_attach(struct ged *gedp,
				BRLObolViewController *controller,
				int sync_current_scene = 1);

GED_EXPORT int
ged_draw_obol_controller_attach_for_view(struct ged *gedp,
	void *view_ctx,
	BRLObolViewController *controller,
	int sync_current_scene = 1);

/**
 * Stop mirroring GED draw transactions to the currently attached controller.
 */
GED_EXPORT void
ged_draw_obol_controller_detach(struct ged *gedp);

/**
 * Return the borrowed controller currently attached to @p gedp, or NULL.
 */
GED_EXPORT BRLObolViewController *
ged_draw_obol_controller(struct ged *gedp);

/**
 * Mirror one completed GED draw transaction into @p controller.
 *
 * This compatibility wrapper delegates to ged_draw_obol_scene_sync_transaction.
 * If @p controller is NULL, the currently attached scene controller is used.
 * Returns 1 if the Obol scene changed, 0 otherwise.
 */
GED_EXPORT int
ged_draw_obol_sync_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    BRLObolViewController *controller = NULL);

/**
 * Rebuild @p controller's database-source scene from the current GED draw set.
 *
 * This compatibility wrapper delegates to ged_draw_obol_scene_sync_full_scene.
 * If @p controller is NULL, the currently attached scene controller is used.  A
 * zero @p source_revision requests a folded value derived from GED's current
 * draw scene revision.  Returns 1 if the Obol scene changed, 0 otherwise.
 */
GED_EXPORT int
ged_draw_obol_sync_full_scene(struct ged *gedp,
			      void *view_ctx,
			      uint32_t source_revision = 0,
			      BRLObolViewController *controller = NULL);

#endif /* __cplusplus */

#endif /* GED_DRAW_OBOL_H */
