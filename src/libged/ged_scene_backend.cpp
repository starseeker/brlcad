/*                G E D _ S C E N E _ B A C K E N D . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged_scene_backend.cpp
 *
 * Owner-scoped dispatch and lifetime management for the private semantic
 * scene-backend contract.  Concrete default rendering is isolated behind the
 * narrow adapter declared in ged_scene_backend_obol_private.h.
 */

#include "common.h"

#include <new>

#include "./ged_private.h"
#include "./ged_scene_backend_private.h"
#include "./ged_scene_backend_obol_private.h"

extern "C" int
ged_scene_backend_apply_private(
    struct ged *gedp,
    const struct ged_draw_transaction *transaction,
    const struct ged_draw_transaction_result *result)
{
    if (!gedp || !transaction)
	return 0;
    Ged_Internal *state = gedp->i && gedp->i->i ? gedp->i->i : nullptr;
    if (!state || state->scene_backend_transition)
	return 0;

    const int changed = state->scene_backend_ops ?
	(state->scene_backend_ops->apply ?
	    state->scene_backend_ops->apply(gedp, transaction, result,
		state->scene_backend_data) : 0) :
	ged_scene_backend_obol_apply_private(gedp, transaction, result);
    return changed;
}


extern "C" int
ged_scene_backend_snapshot_private(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return 0;
    Ged_Internal *state = gedp->i->i;
    if (state->scene_backend_transition)
	return 0;

    const int changed = state->scene_backend_ops ?
	(state->scene_backend_ops->snapshot ?
	    state->scene_backend_ops->snapshot(gedp,
		state->scene_backend_data) : 0) :
	ged_scene_backend_obol_snapshot_private(gedp);
    return changed;
}


extern "C" int
ged_scene_backend_selection_private(
    struct ged *gedp,
    const char *const *added_paths, size_t added_count,
    const char *const *removed_paths, size_t removed_count,
    const char *const *selected_paths, size_t selected_count)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return 0;
    Ged_Internal *state = gedp->i->i;
    if (state->scene_backend_transition)
	return 0;

    const int changed = state->scene_backend_ops ?
	(state->scene_backend_ops->selection ?
	    state->scene_backend_ops->selection(gedp, added_paths,
		added_count, removed_paths, removed_count, selected_paths,
		selected_count, state->scene_backend_data) : 0) :
	ged_scene_backend_obol_selection_private(gedp, added_paths,
	    added_count, removed_paths, removed_count, selected_paths,
	    selected_count);
    return changed;
}


extern "C" int
ged_scene_backend_view_policy_private(
    struct ged *gedp, struct ged_view_context *view)
{
    if (!gedp || !view || !gedp->i || !gedp->i->i)
	return 0;
    Ged_Internal *state = gedp->i->i;
    if (state->scene_backend_transition)
	return 0;

    const int changed = state->scene_backend_ops ?
	(state->scene_backend_ops->view_policy ?
	    state->scene_backend_ops->view_policy(gedp, view,
		state->scene_backend_data) : 0) :
	ged_scene_backend_obol_view_policy_private(gedp, view);
    return changed;
}


extern "C" void
ged_scene_backend_detach_private(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;

    Ged_Internal *state = gedp->i->i;
    if (state->scene_backend_transition)
	return;

    struct ged_scene_backend_ops *old_operations = state->scene_backend_ops;
    void *old_data = state->scene_backend_data;

    /* Clear the binding before invoking client teardown.  Reducer activity
     * caused by teardown therefore cannot dispatch back into an adapter whose
     * lifetime is ending. */
    state->scene_backend_transition = true;
    state->scene_backend_ops = nullptr;
    state->scene_backend_data = nullptr;
    if (old_operations && old_operations->detach)
	old_operations->detach(gedp, old_data);
    delete old_operations;
    state->scene_backend_transition = false;
}


extern "C" void
ged_scene_backend_set_private(
    struct ged *gedp,
    const struct ged_scene_backend_ops *operations,
    void *client_data)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;

    Ged_Internal *state = gedp->i->i;
    if (state->scene_backend_transition)
	return;

    struct ged_scene_backend_ops *operations_copy = operations ?
	new (std::nothrow) ged_scene_backend_ops(*operations) : nullptr;
    if (operations && !operations_copy)
	return;

    ged_scene_backend_detach_private(gedp);

    if (!operations) {
	(void)ged_scene_backend_snapshot_private(gedp);
	return;
    }

    state->scene_backend_ops = operations_copy;
    state->scene_backend_data = client_data;
    (void)ged_scene_backend_snapshot_private(gedp);
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
