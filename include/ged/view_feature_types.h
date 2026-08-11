/*             V I E W _ F E A T U R E _ T Y P E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/view_feature_types.h
 *
 * Renderer-neutral retained view-feature value types.  Database draw records,
 * polygon editing records, and renderer realization state deliberately do not
 * belong to this domain.
 */

#ifndef GED_VIEW_FEATURE_TYPES_H
#define GED_VIEW_FEATURE_TYPES_H

#include <stddef.h>
#include <stdint.h>

#include "vmath.h"
#include "ged/defines.h"

struct bg_line_layer_builder;
struct bu_vls;
struct bv_axes_state;
struct bv_view_info;
struct ged;
struct ged_view_context;

__BEGIN_DECLS

struct ged_view_feature_style {
    int visible;
    int selectable;
    int color_valid;
    unsigned char color[3];
    int line_width;
    int line_style;
    int arrow;
    fastf_t arrow_tip_length;
    fastf_t arrow_tip_width;
};

#define GED_VIEW_FEATURE_STYLE_INIT { -1, -1, 0, {0, 0, 0}, -1, -1, -1, -1.0, -1.0 }

struct ged_view_feature_label {
    const char *text;
    point_t point;
    int color_valid;
    unsigned char color[3];
    int line_flag;
    point_t target;
    int anchor;
    int arrow;
    fastf_t font_size;
};

/** One independently styled line layer in a retained feature batch. */
struct ged_view_feature_line_layer {
    const char *name;
    const point_t *points;
    const int *commands;
    size_t point_count;
    struct ged_view_feature_style style;
};

#define GED_VIEW_FEATURE_LINE_LAYER_INIT { NULL, NULL, NULL, 0, GED_VIEW_FEATURE_STYLE_INIT }

struct ged_view_feature_batch;

enum ged_view_feature_batch_status {
    GED_VIEW_FEATURE_BATCH_NONE = 0,
    GED_VIEW_FEATURE_BATCH_ACCEPTED,
    GED_VIEW_FEATURE_BATCH_UPDATED,
    GED_VIEW_FEATURE_BATCH_REMOVED,
    GED_VIEW_FEATURE_BATCH_FAILED
};

struct ged_view_feature_batch_event {
    enum ged_view_feature_batch_status status;
    const char *feature_name;
    const char *command;
    const char *diagnostic;
    uint64_t feature_id;
    uint64_t feature_revision;
};

#define GED_VIEW_FEATURE_BATCH_EVENT_INIT { GED_VIEW_FEATURE_BATCH_NONE, NULL, NULL, NULL, 0, 0 }

typedef void (*ged_view_feature_batch_callback)(
    const struct ged_view_feature_batch_event *result,
    void *data);

struct ged_view_feature_batch_desc {
    const char *owner_id;
    const char *owner_role;
    const char *run_id;
    uint64_t generation;
    int local;
    int overlay_class;
    int lifecycle;
    int overlay_order;
    ged_view_feature_batch_callback event_cb;
    void *event_cb_data;
};

#define GED_VIEW_FEATURE_BATCH_DESC_INIT { NULL, NULL, NULL, 0, 0, GED_VIEW_FEATURE_OVERLAY_CLASS_COMMAND_RESULT, GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND, GED_VIEW_FEATURE_OVERLAY_ORDER_POST_TRANSPARENT, NULL, NULL }

struct ged_view_feature_metadata {
    const char *key;
    const char *value;
};

#define GED_VIEW_FEATURE_METADATA_INIT { NULL, NULL }

/** Renderer-neutral commands associated with retained line-set points. */
enum ged_draw_view_line_command {
    GED_DRAW_VIEW_LINE_MOVE = 0,
    GED_DRAW_VIEW_LINE_DRAW = 1,
    GED_DRAW_VIEW_LINE_POINT_DRAW = 12
};

enum ged_view_feature_kind {
    GED_VIEW_FEATURE_KIND_UNKNOWN = 0,
    GED_VIEW_FEATURE_KIND_LINES,
    GED_VIEW_FEATURE_KIND_INDEXED_LINES,
    GED_VIEW_FEATURE_KIND_POINTS,
    GED_VIEW_FEATURE_KIND_LABELS,
    GED_VIEW_FEATURE_KIND_ARROW,
    GED_VIEW_FEATURE_KIND_AXES,
    GED_VIEW_FEATURE_KIND_LINE_LAYER,
    GED_VIEW_FEATURE_KIND_EDIT_PREVIEW,
    GED_VIEW_FEATURE_KIND_INDEXED_FACE_SET,
    GED_VIEW_FEATURE_KIND_POLYGON_OVERLAY,
    GED_VIEW_FEATURE_KIND_HUD_LABEL,
    GED_VIEW_FEATURE_KIND_CUSTOM_NODE
};

enum ged_view_feature_scope {
    GED_VIEW_FEATURE_SCOPE_UNKNOWN = 0,
    GED_VIEW_FEATURE_SCOPE_SHARED,
    GED_VIEW_FEATURE_SCOPE_LOCAL
};

enum ged_view_feature_overlay_class {
    GED_VIEW_FEATURE_OVERLAY_CLASS_UNKNOWN = 0,
    GED_VIEW_FEATURE_OVERLAY_CLASS_NONE,
    GED_VIEW_FEATURE_OVERLAY_CLASS_FACEPLATE,
    GED_VIEW_FEATURE_OVERLAY_CLASS_EDIT_HANDLE,
    GED_VIEW_FEATURE_OVERLAY_CLASS_MEASURE,
    GED_VIEW_FEATURE_OVERLAY_CLASS_SELECTION_RUBBER_BAND,
    GED_VIEW_FEATURE_OVERLAY_CLASS_SNAP_GUIDE,
    GED_VIEW_FEATURE_OVERLAY_CLASS_COMMAND_RESULT,
    GED_VIEW_FEATURE_OVERLAY_CLASS_DIAGNOSTIC,
    GED_VIEW_FEATURE_OVERLAY_CLASS_TCL_OVERLAY,
    GED_VIEW_FEATURE_OVERLAY_CLASS_POLYGON_EDIT,
    GED_VIEW_FEATURE_OVERLAY_CLASS_SKETCH_EDIT,
    GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION
};

enum ged_view_feature_lifecycle {
    GED_VIEW_FEATURE_LIFECYCLE_UNKNOWN = 0,
    GED_VIEW_FEATURE_LIFECYCLE_NONE,
    GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT,
    GED_VIEW_FEATURE_LIFECYCLE_PER_FRAME,
    GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND,
    GED_VIEW_FEATURE_LIFECYCLE_PER_TOOL,
    GED_VIEW_FEATURE_LIFECYCLE_PER_VIEW,
    GED_VIEW_FEATURE_LIFECYCLE_SHARED_VIEW_SET,
    GED_VIEW_FEATURE_LIFECYCLE_AUTO_REMOVE_ON_SOURCE
};

enum ged_view_feature_overlay_order {
    GED_VIEW_FEATURE_OVERLAY_ORDER_UNKNOWN = 0,
    GED_VIEW_FEATURE_OVERLAY_ORDER_MODEL,
    GED_VIEW_FEATURE_OVERLAY_ORDER_SCREEN,
    GED_VIEW_FEATURE_OVERLAY_ORDER_XRAY,
    GED_VIEW_FEATURE_OVERLAY_ORDER_POST_TRANSPARENT
};

#define GED_VIEW_FEATURE_OWNER_ID_MAX 64
#define GED_VIEW_FEATURE_OWNER_ROLE_MAX 64

struct ged_view_feature_summary {
    int exists;
    int is_overlay;
    int is_label;
    int is_transient_preview;
    int is_command_result;
    int visible;
    int kind;
    int scope;
    int overlay_class;
    int lifecycle;
    unsigned char color[3];
    size_t child_count;
    size_t geometry_command_count;
    size_t metadata_count;
    size_t primitive_metadata_count;
    size_t selected_primitive_count;
    size_t highlighted_primitive_count;
    uint64_t owner_generation;
    char owner_id[GED_VIEW_FEATURE_OWNER_ID_MAX];
    char owner_role[GED_VIEW_FEATURE_OWNER_ROLE_MAX];
};

#define GED_VIEW_FEATURE_SUMMARY_INIT { 0, 0, 0, 0, 0, 0, GED_VIEW_FEATURE_KIND_UNKNOWN, GED_VIEW_FEATURE_SCOPE_UNKNOWN, GED_VIEW_FEATURE_OVERLAY_CLASS_UNKNOWN, GED_VIEW_FEATURE_LIFECYCLE_UNKNOWN, {0, 0, 0}, 0, 0, 0, 0, 0, 0, 0, {'\0'}, {'\0'} }

__END_DECLS

#endif /* GED_VIEW_FEATURE_TYPES_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
