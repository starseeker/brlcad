/*         G E D _ S C E N E _ B A C K E N D _ O B O L _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_scene_backend_obol_private.h
 *
 * Narrow private entry points implemented by the default Obol adapter.  The
 * semantic reducer and public scene API do not include renderer headers.
 */

#ifndef LIBGED_GED_SCENE_BACKEND_OBOL_PRIVATE_H
#define LIBGED_GED_SCENE_BACKEND_OBOL_PRIVATE_H

#include "common.h"

#include "./ged_scene_reducer_private.h"

__BEGIN_DECLS

extern int ged_scene_backend_obol_apply_private(
    struct ged *gedp,
    const struct ged_draw_transaction *transaction,
    const struct ged_draw_transaction_result *result);

extern int ged_scene_backend_obol_snapshot_private(struct ged *gedp);

extern int ged_scene_backend_obol_selection_private(
    struct ged *gedp,
    const char *const *added_paths, size_t added_count,
    const char *const *removed_paths, size_t removed_count,
    const char *const *selected_paths, size_t selected_count);

extern int ged_scene_backend_obol_view_policy_private(
    struct ged *gedp, struct ged_view_context *view);

__END_DECLS

#endif /* LIBGED_GED_SCENE_BACKEND_OBOL_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
