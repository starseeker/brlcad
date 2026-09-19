/*             G E D _ S C E N E _ B A C K E N D _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_scene_backend_private.h
 *
 * Private semantic-scene/backend boundary.  The reducer knows only this
 * contract; concrete Obol synchronization remains behind its default
 * implementation.  Tests may install a recording backend on one GED owner.
 */

#ifndef LIBGED_GED_SCENE_BACKEND_PRIVATE_H
#define LIBGED_GED_SCENE_BACKEND_PRIVATE_H

#include "common.h"

#include "./ged_scene_reducer_private.h"

__BEGIN_DECLS

struct ged_scene_backend_ops {
    int (*apply)(struct ged *gedp,
	const struct ged_scene_reducer_request *transaction,
	const struct ged_scene_reducer_result *result,
	void *client_data);
    int (*snapshot)(struct ged *gedp, void *client_data);
    int (*selection)(struct ged *gedp,
	const char *const *added_paths, size_t added_count,
	const char *const *removed_paths, size_t removed_count,
	const char *const *selected_paths, size_t selected_count,
	void *client_data);
    int (*view_policy)(struct ged *gedp,
	struct ged_view_context *view, void *client_data);
    void (*detach)(struct ged *gedp, void *client_data);
};

extern int ged_scene_backend_apply_private(
    struct ged *gedp,
    const struct ged_scene_reducer_request *transaction,
    const struct ged_scene_reducer_result *result);

extern int ged_scene_backend_snapshot_private(struct ged *gedp);

extern int ged_scene_backend_selection_private(
    struct ged *gedp,
    const char *const *added_paths, size_t added_count,
    const char *const *removed_paths, size_t removed_count,
    const char *const *selected_paths, size_t selected_count);

/* Notify the active adapter after a renderer-neutral view policy changes. */
extern int ged_scene_backend_view_policy_private(
    struct ged *gedp, struct ged_view_context *view);

/* Tear down an installed override without activating or snapshotting the
 * default adapter.  This is reserved for GED owner destruction. */
extern void ged_scene_backend_detach_private(struct ged *gedp);

/* This one entry point is the deliberately exported in-tree backend adapter
 * hook.  Libged copies the operations table, but client_data remains borrowed
 * until detach is called.  Replacing an adapter detaches the old adapter, then
 * gives the new adapter exactly one snapshot.  Passing NULL restores the
 * private default adapter and snapshots it if it has an attached endpoint.
 * Calls are owner-thread only; detach callbacks must not replace the adapter
 * reentrantly.  The apply/snapshot/selection calls remain hidden reducer
 * details. */
COMPILER_DLLEXPORT extern void ged_scene_backend_set_private(
    struct ged *gedp,
    const struct ged_scene_backend_ops *operations,
    void *client_data);

__END_DECLS

#endif /* LIBGED_GED_SCENE_BACKEND_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
