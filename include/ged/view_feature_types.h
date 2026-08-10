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
struct bg_tess_tol;
struct bn_tol;
struct bu_vls;
struct bv_axes_state;
struct bv_view_info;
struct db_i;
struct ged;
struct ged_view_context;
struct rt_db_internal;

__BEGIN_DECLS

struct ged_draw_view_line_style {
    int color[3];
    int line_width;
};

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

/**
 * Opaque value reference to a managed view feature.
 *
 * The values identify the issuing view owner, object, and backing-store
 * generation.  A reference does not retain its owner and contains no pointer.
 */
typedef struct ged_view_feature_ref {
    uint64_t owner;
    uint64_t id;
    uint64_t generation;
} ged_view_feature_ref;

#define GED_VIEW_FEATURE_REF_NULL_INIT {0, 0, 0}
#ifdef __cplusplus
#  define GED_VIEW_FEATURE_REF_NULL ged_view_feature_ref{0, 0, 0}
#else
#  define GED_VIEW_FEATURE_REF_NULL ((ged_view_feature_ref){0, 0, 0})
#endif

enum ged_view_feature_family {
    GED_VIEW_FEATURE_UNKNOWN = 0,
    GED_VIEW_FEATURE_TRANSIENT_PREVIEW = 1
};

enum ged_view_edit_preview_event {
    GED_VIEW_EDIT_PREVIEW_BEGIN = 0,
    GED_VIEW_EDIT_PREVIEW_UPDATE,
    GED_VIEW_EDIT_PREVIEW_COMMIT,
    GED_VIEW_EDIT_PREVIEW_CANCEL,
    GED_VIEW_EDIT_PREVIEW_REPLACE_SOURCE,
    GED_VIEW_EDIT_PREVIEW_DISCARD
};

struct ged_view_edit_transaction {
    enum ged_view_edit_preview_event event;
    ged_view_feature_ref feature;
    const char *feature_name;
    const void *owner;
    const char *source_path;
    const char *edit_intent_id;
    const char *edit_intent_role;
    const point_t *points;
    const int *commands;
    size_t point_count;
    struct db_i *dbip;
    struct rt_db_internal *internal;
    const fastf_t *matrix;
    const struct bg_tess_tol *ttol;
    const struct bn_tol *tol;
    uint32_t source_revision;
    uint32_t inputs_revision;
    int color_valid;
    unsigned char color[3];
};

#define GED_VIEW_EDIT_TRANSACTION_INIT { \
    GED_VIEW_EDIT_PREVIEW_UPDATE, GED_VIEW_FEATURE_REF_NULL, \
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, \
    NULL, 0, 0, 0, {255, 255, 255} }

struct ged_view_feature_label {
    const char *text;
    point_t point;
    int color_valid;
    unsigned char color[3];
    fastf_t font_size;
};

struct ged_annotation_label {
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

#define GED_ANNOTATION_LABEL_INIT { NULL, VINIT_ZERO, 0, {0, 0, 0}, 0, VINIT_ZERO, 0, 0, 0.0 }

struct ged_annotation_axes {
    point_t position;
    fastf_t size;
    int line_width;
    int color[3];
};

struct ged_annotation_line_layer {
    const char *name;
    const point_t *points;
    const int *commands;
    size_t point_count;
    struct ged_view_feature_style style;
};

#define GED_ANNOTATION_LINE_LAYER_INIT { NULL, NULL, NULL, 0, GED_VIEW_FEATURE_STYLE_INIT }

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
    ged_view_feature_batch_callback event_cb;
    void *event_cb_data;
};

#define GED_VIEW_FEATURE_BATCH_DESC_INIT { NULL, NULL, NULL, 0, 0, NULL, NULL }

struct ged_view_feature_metadata {
    const char *key;
    const char *value;
};

#define GED_VIEW_FEATURE_METADATA_INIT { NULL, NULL }

typedef int (*ged_view_feature_depth_cb)(fastf_t depth, void *data);

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
