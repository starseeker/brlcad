/*             G E D _ S C E N E _ R E D U C E R _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_scene_reducer_private.h
 *
 * Private migration representation used inside the semantic scene reducer
 * and Obol adapter.  Public clients use the typed ged_scene_* requests and
 * committed ged_scene_delta accessors.
 */

#ifndef LIBGED_GED_SCENE_REDUCER_PRIVATE_H
#define LIBGED_GED_SCENE_REDUCER_PRIVATE_H

#include "common.h"

#include <stdint.h>

#include "bu/vls.h"
#include "./ged_draw_types_private.h"
#include "./ged_scene_record_private.h"
#include "ged/scene_types.h"

__BEGIN_DECLS

typedef enum ged_draw_transaction_kind {
    GED_DRAW_TXN_NONE = 0,
    GED_DRAW_TXN_DRAW,
    GED_DRAW_TXN_ERASE,
    GED_DRAW_TXN_ERASE_PREFIX,
    GED_DRAW_TXN_REDRAW,
    GED_DRAW_TXN_VISIBILITY,
    GED_DRAW_TXN_HIGHLIGHT,
    GED_DRAW_TXN_HIGHLIGHTS_CLEAR,
    GED_DRAW_TXN_HIGHLIGHT_OCCURRENCE,
    GED_DRAW_TXN_MATERIAL_CHANGED,
    GED_DRAW_TXN_TRANSPARENCY,
    GED_DRAW_TXN_DEFAULT_DRAW_MODE,
    GED_DRAW_TXN_STALE_SOURCE,
    GED_DRAW_TXN_SOURCE_UPDATED,
    GED_DRAW_TXN_SOURCE_RENAMED,
    GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED,
    GED_DRAW_TXN_CLEAR,
    GED_DRAW_TXN_CLEAR_SCOPE,
    GED_DRAW_TXN_TEARDOWN
} ged_draw_transaction_kind;

struct ged_draw_transaction {
    ged_draw_transaction_kind kind;
    const char *path;
    const char *new_path;
    const char **paths;
    int path_count;
    const struct db_full_path *dfp;
    struct ged_view_context *view;
    int mode;
    enum ged_scene_path_match match;
    ged_draw_shape_ref shape_ref;
    const void *appearance;
    int autoview;
    double value;
    ged_draw_stale_reason stale_reason;
    int removed;
    int redraw;
};

struct ged_draw_transaction_result {
    int status;
    int affected_groups;
    int affected_shapes;
    int stale_count;
    int redrawn_count;
    int error_count;
    int presentation_only;
    int progressive_data_complete;
    uint64_t scene_revision_before;
    uint64_t scene_revision_after;
    struct bu_vls names;
    struct bu_vls errors;
};

extern struct ged_draw_transaction ged_draw_transaction_make(
    ged_draw_transaction_kind kind, const char *path);
extern struct ged_draw_transaction ged_draw_transaction_make_value(
    ged_draw_transaction_kind kind, const char *path, double value);
extern void ged_draw_transaction_result_init(
    struct ged_draw_transaction_result *result);
extern void ged_draw_transaction_result_free(
    struct ged_draw_transaction_result *result);
extern int ged_draw_apply_transaction(struct ged *gedp,
    const struct ged_draw_transaction *txn,
    struct ged_draw_transaction_result *result);

__END_DECLS

#endif /* LIBGED_GED_SCENE_REDUCER_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
