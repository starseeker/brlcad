/*          G E D _ S C E N E _ B A C K E N D _ O B O L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged_scene_backend_obol.cpp
 *
 * Default renderer adapter for the semantic scene-backend contract.
 */

#include "common.h"

#include "./ged_draw_private.h"
#include "./ged_scene_backend_obol_private.h"

extern "C" int
ged_scene_backend_obol_apply_private(
    struct ged *gedp,
    const struct ged_scene_reducer_request *transaction,
    const struct ged_scene_reducer_result *result)
{
    return ged_draw_obol_backend_apply_transaction(gedp, transaction, result);
}


extern "C" int
ged_scene_backend_obol_snapshot_private(struct ged *gedp)
{
    return ged_draw_obol_scene_sync_full_scene(gedp, nullptr, 0, nullptr);
}


extern "C" int
ged_scene_backend_obol_selection_private(
    struct ged *gedp,
    const char *const *added_paths, size_t added_count,
    const char *const *removed_paths, size_t removed_count,
    const char *const *selected_paths, size_t selected_count)
{
    return ged_draw_obol_database_sources_apply_selection_delta(gedp,
	added_paths, added_count, removed_paths, removed_count,
	selected_paths, selected_count);
}


extern "C" int
ged_scene_backend_obol_view_policy_private(
    struct ged *gedp, struct ged_view_context *view)
{
    return ged_draw_obol_view_lod_policy_changed(gedp, view);
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
