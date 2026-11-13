/*             G E D _ S C E N E _ R E D U C E R _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_scene_reducer_private.h
 *
 * Private request and result types for the semantic scene reducer and backend
 * adapters.  Public clients use typed ged_scene_* requests and committed
 * ged_scene_delta accessors.
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

typedef enum ged_scene_reducer_operation {
    GED_SCENE_REDUCER_NONE = 0,
    GED_SCENE_REDUCER_DRAW,
    GED_SCENE_REDUCER_ERASE,
    GED_SCENE_REDUCER_ERASE_PREFIX,
    GED_SCENE_REDUCER_REDRAW,
    GED_SCENE_REDUCER_VISIBILITY,
    GED_SCENE_REDUCER_HIGHLIGHT,
    GED_SCENE_REDUCER_HIGHLIGHTS_CLEAR,
    GED_SCENE_REDUCER_HIGHLIGHT_OCCURRENCE,
    GED_SCENE_REDUCER_MATERIAL_CHANGED,
    GED_SCENE_REDUCER_TRANSPARENCY,
    GED_SCENE_REDUCER_DEFAULT_DRAW_MODE,
    GED_SCENE_REDUCER_STALE_SOURCE,
    GED_SCENE_REDUCER_SOURCE_UPDATED,
    GED_SCENE_REDUCER_SOURCE_RENAMED,
    GED_SCENE_REDUCER_SOURCE_REFERENCES_REMOVED,
    GED_SCENE_REDUCER_CLEAR,
    GED_SCENE_REDUCER_CLEAR_SCOPE,
    GED_SCENE_REDUCER_TEARDOWN
} ged_scene_reducer_operation;

struct ged_scene_reducer_request {
    ged_scene_reducer_operation kind;
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

struct ged_scene_reducer_result {
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
    char **paths;
    size_t path_count;
    size_t path_capacity;
    struct bu_vls errors;
};

extern struct ged_scene_reducer_request ged_scene_reducer_request_make(
    ged_scene_reducer_operation kind, const char *path);
extern struct ged_scene_reducer_request ged_scene_reducer_request_make_value(
    ged_scene_reducer_operation kind, const char *path, double value);
extern void ged_scene_reducer_result_init(
    struct ged_scene_reducer_result *result);
extern void ged_scene_reducer_result_free(
    struct ged_scene_reducer_result *result);
extern int ged_scene_reduce(struct ged *gedp,
    const struct ged_scene_reducer_request *txn,
    struct ged_scene_reducer_result *result);

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
