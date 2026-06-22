/*                    Q G L E G A C Y V I E W . H
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
/** @file QgLegacyView.h
 *
 * Opaque transitional qtcad view handle.
 */

#ifndef QGLEGACYVIEW_H
#define QGLEGACYVIEW_H

#include <stdint.h>

#include "bg/polygon.h"
#include "qtcad/defines.h"
#include "rt/view.h"
#include "vmath.h"

struct bu_color;
struct bu_vls;
struct db_i;
struct directory;
struct ged;
struct qg_legacy_view;
struct qg_legacy_view_draw_appearance;
struct qg_legacy_view_bounds_update_state;
struct qg_legacy_view_pick_result;
typedef struct qg_legacy_view_draw_transaction
	qg_legacy_view_draw_transaction;
typedef struct qg_legacy_view_draw_transaction_result
	qg_legacy_view_draw_transaction_result;

typedef int (*qg_legacy_view_polygon_record_callback_t)(
	rt_view_polygon_ref ref,
	const struct rt_view_polygon_record *record,
	void *data);

typedef void (*qg_legacy_view_selection_path_callback_t)(
	const char *path,
	void *data);

typedef void (*qg_legacy_view_draw_observer_callback_t)(
	struct ged *gedp,
	const qg_legacy_view_draw_transaction *txn,
	const qg_legacy_view_draw_transaction_result *result,
	void *client_data);

typedef enum qg_legacy_view_draw_transaction_kind {
    QG_LEGACY_VIEW_DRAW_TXN_NONE = 0,
    QG_LEGACY_VIEW_DRAW_TXN_DRAW,
    QG_LEGACY_VIEW_DRAW_TXN_ERASE,
    QG_LEGACY_VIEW_DRAW_TXN_ERASE_PREFIX,
    QG_LEGACY_VIEW_DRAW_TXN_REDRAW,
    QG_LEGACY_VIEW_DRAW_TXN_VISIBILITY,
    QG_LEGACY_VIEW_DRAW_TXN_HIGHLIGHT,
    QG_LEGACY_VIEW_DRAW_TXN_MATERIAL_CHANGED,
    QG_LEGACY_VIEW_DRAW_TXN_REFRESH_MATERIAL_COLORS,
    QG_LEGACY_VIEW_DRAW_TXN_TRANSPARENCY,
    QG_LEGACY_VIEW_DRAW_TXN_DEFAULT_DRAW_MODE,
    QG_LEGACY_VIEW_DRAW_TXN_STALE_SOURCE,
    QG_LEGACY_VIEW_DRAW_TXN_SOURCE_UPDATED,
    QG_LEGACY_VIEW_DRAW_TXN_SOURCE_RENAMED,
    QG_LEGACY_VIEW_DRAW_TXN_SOURCE_REFERENCES_REMOVED,
    QG_LEGACY_VIEW_DRAW_TXN_CLEAR,
    QG_LEGACY_VIEW_DRAW_TXN_CLEAR_SCOPE,
    QG_LEGACY_VIEW_DRAW_TXN_TEARDOWN
} qg_legacy_view_draw_transaction_kind;

#define QG_LEGACY_VIEW_DRAW_TRANSACTION_STORAGE_WORDS 16

struct qg_legacy_view_draw_transaction {
    uintptr_t storage[QG_LEGACY_VIEW_DRAW_TRANSACTION_STORAGE_WORDS];
};

struct qg_legacy_view_draw_transaction_result {
    uintptr_t impl;
};

typedef struct qg_legacy_view_feature_ref {
    uintptr_t token;
    uint64_t revision;
} qg_legacy_view_feature_ref;

#define QG_LEGACY_VIEW_FEATURE_REF_NULL_INIT {0, 0}

#ifdef __cplusplus
#  define QG_LEGACY_VIEW_FEATURE_REF_NULL qg_legacy_view_feature_ref{0, 0}
#else
#  define QG_LEGACY_VIEW_FEATURE_REF_NULL ((qg_legacy_view_feature_ref){0, 0})
#endif

typedef enum qg_legacy_view_feature_family {
    QG_LEGACY_VIEW_FEATURE_UNKNOWN = 0,
    QG_LEGACY_VIEW_FEATURE_TRANSIENT_PREVIEW = 1
} qg_legacy_view_feature_family;

typedef enum qg_legacy_view_edit_preview_event {
    QG_LEGACY_VIEW_EDIT_PREVIEW_BEGIN = 0,
    QG_LEGACY_VIEW_EDIT_PREVIEW_UPDATE,
    QG_LEGACY_VIEW_EDIT_PREVIEW_COMMIT,
    QG_LEGACY_VIEW_EDIT_PREVIEW_CANCEL,
    QG_LEGACY_VIEW_EDIT_PREVIEW_REPLACE_SOURCE,
    QG_LEGACY_VIEW_EDIT_PREVIEW_DISCARD
} qg_legacy_view_edit_preview_event;

typedef struct qg_legacy_view_edit_preview_callbacks {
    uint64_t (*revision_cb)(void *);
    int (*update_cb)(void *);
    int (*pick_cb)(void *, int, int, void *);
} qg_legacy_view_edit_preview_callbacks;

#define QG_LEGACY_VIEW_EDIT_PREVIEW_CALLBACKS_INIT { 0, 0, 0 }

typedef struct qg_legacy_view_feature_label {
    const char *text;
    point_t point;
    int color_valid;
    unsigned char color[3];
} qg_legacy_view_feature_label;

#define QG_LEGACY_VIEW_SNAP_KIND_GRID 0x00000001u
#define QG_LEGACY_VIEW_SNAP_KIND_ENDPOINT 0x00000002u
#define QG_LEGACY_VIEW_SNAP_KIND_MIDPOINT 0x00000004u
#define QG_LEGACY_VIEW_SNAP_KIND_INTERSECTION 0x00000008u
#define QG_LEGACY_VIEW_SNAP_KIND_PERPENDICULAR 0x00000010u
#define QG_LEGACY_VIEW_SNAP_KIND_TANGENT 0x00000020u
#define QG_LEGACY_VIEW_SNAP_KIND_OVERLAY_HANDLE 0x00000040u

#define QG_LEGACY_VIEW_REFRESH_VIEW 0x00000001u
#define QG_LEGACY_VIEW_REFRESH_DRAW 0x00000002u
#define QG_LEGACY_VIEW_REFRESH_EDIT 0x00000004u
#define QG_LEGACY_VIEW_REFRESH_FRAMEBUFFER 0x00000008u
#define QG_LEGACY_VIEW_REFRESH_OVERLAY 0x00000010u
#define QG_LEGACY_VIEW_REFRESH_FORCE 0x80000000u
#define QG_LEGACY_VIEW_REFRESH_ALL 0xffffffffu

#define QG_LEGACY_VIEW_DRAW_MODE_ANY -1
#define QG_LEGACY_VIEW_DRAW_MODE_WIRE 0
#define QG_LEGACY_VIEW_DRAW_MODE_SHADED_BOTS 1
#define QG_LEGACY_VIEW_DRAW_MODE_SHADED 2
#define QG_LEGACY_VIEW_DRAW_MODE_EVAL_WIRE 3
#define QG_LEGACY_VIEW_DRAW_MODE_HIDDEN_LINE 4
#define QG_LEGACY_VIEW_DRAW_MODE_EVAL_POINTS 5

#define QG_LEGACY_VIEW_ADJUST_IDLE 0x000ULL
#define QG_LEGACY_VIEW_ADJUST_ROT 0x001ULL
#define QG_LEGACY_VIEW_ADJUST_TRANS 0x002ULL
#define QG_LEGACY_VIEW_ADJUST_SCALE 0x004ULL
#define QG_LEGACY_VIEW_ADJUST_CENTER 0x008ULL
#define QG_LEGACY_VIEW_ADJUST_CON_X 0x010ULL
#define QG_LEGACY_VIEW_ADJUST_CON_Y 0x020ULL
#define QG_LEGACY_VIEW_ADJUST_CON_Z 0x040ULL
#define QG_LEGACY_VIEW_ADJUST_CON_GRID 0x080ULL
#define QG_LEGACY_VIEW_ADJUST_CON_LINES 0x100ULL

QTCAD_EXPORT extern qg_legacy_view *qg_legacy_view_local_create(
	const char *name);

QTCAD_EXPORT extern void qg_legacy_view_local_free(qg_legacy_view *view);

QTCAD_EXPORT extern void qg_legacy_view_local_destroy(qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_dimensions_set(qg_legacy_view *view,
	int width,
	int height);

QTCAD_EXPORT extern int qg_legacy_view_width_get(const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_height_get(const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_unit_conversion_set(qg_legacy_view *view,
	double local2base,
	double base2local);

QTCAD_EXPORT extern int qg_legacy_view_name_set(qg_legacy_view *view,
	const char *name);

QTCAD_EXPORT extern qg_legacy_view *qg_legacy_view_ged_active_get(
	struct ged *gedp);

QTCAD_EXPORT extern void qg_legacy_view_ged_active_set(struct ged *gedp,
	qg_legacy_view *view);

QTCAD_EXPORT extern struct db_i *qg_legacy_view_ged_database(
	struct ged *gedp);

QTCAD_EXPORT extern int qg_legacy_view_ged_view_set_add(struct ged *gedp,
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_ged_view_set_remove(struct ged *gedp,
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_ged_view_set_attach(struct ged *gedp,
	qg_legacy_view *view);

QTCAD_EXPORT extern qg_legacy_view *qg_legacy_view_draw_transaction_view_get(
	const qg_legacy_view_draw_transaction *txn);

QTCAD_EXPORT extern void qg_legacy_view_draw_transaction_view_set(
	qg_legacy_view_draw_transaction *txn,
	qg_legacy_view *view);

QTCAD_EXPORT extern void qg_legacy_view_draw_transaction_init(
	qg_legacy_view_draw_transaction *txn,
	qg_legacy_view_draw_transaction_kind kind,
	const char *path);

QTCAD_EXPORT extern void qg_legacy_view_draw_transaction_paths_set(
	qg_legacy_view_draw_transaction *txn,
	const char **paths,
	int path_count);

QTCAD_EXPORT extern void qg_legacy_view_draw_transaction_autoview_set(
	qg_legacy_view_draw_transaction *txn,
	int autoview);

QTCAD_EXPORT extern qg_legacy_view_draw_appearance *
qg_legacy_view_draw_appearance_create(int draw_mode);

QTCAD_EXPORT extern void qg_legacy_view_draw_appearance_destroy(
	qg_legacy_view_draw_appearance *appearance);

QTCAD_EXPORT extern void qg_legacy_view_draw_transaction_appearance_set(
	qg_legacy_view_draw_transaction *txn,
	const qg_legacy_view_draw_appearance *appearance);

QTCAD_EXPORT extern int qg_legacy_view_draw_path_state(struct ged *gedp,
	qg_legacy_view *view,
	const char *path,
	int mode);

QTCAD_EXPORT extern int qg_legacy_view_draw_has_paths(struct ged *gedp,
	qg_legacy_view *view,
	int mode);

QTCAD_EXPORT extern size_t qg_legacy_view_draw_list_paths(struct ged *gedp,
	qg_legacy_view *view,
	int mode,
	int expanded,
	struct bu_vls *paths);

QTCAD_EXPORT extern uint64_t qg_legacy_view_draw_scene_revision(
	struct ged *gedp);

QTCAD_EXPORT extern int qg_legacy_view_scene_attached(
	const qg_legacy_view *view);

QTCAD_EXPORT extern uintptr_t qg_legacy_view_draw_observer_add(
	struct ged *gedp,
	qg_legacy_view_draw_observer_callback_t callback,
	void *client_data);

QTCAD_EXPORT extern int qg_legacy_view_draw_observer_remove(
	struct ged *gedp,
	uintptr_t token);

QTCAD_EXPORT extern int qg_legacy_view_draw_redraw(struct ged *gedp);

QTCAD_EXPORT extern int qg_legacy_view_draw_transaction_has_view(
	const qg_legacy_view_draw_transaction *txn);

QTCAD_EXPORT extern int qg_legacy_view_draw_transaction_is_view_only(
	const qg_legacy_view_draw_transaction *txn);

QTCAD_EXPORT extern qg_legacy_view_draw_transaction_kind
qg_legacy_view_draw_transaction_kind_get(
	const qg_legacy_view_draw_transaction *txn);

QTCAD_EXPORT extern const char *qg_legacy_view_draw_transaction_path(
	const qg_legacy_view_draw_transaction *txn);

QTCAD_EXPORT extern int qg_legacy_view_draw_transaction_path_count(
	const qg_legacy_view_draw_transaction *txn);

QTCAD_EXPORT extern const char *qg_legacy_view_draw_transaction_path_at(
	const qg_legacy_view_draw_transaction *txn,
	int index);

QTCAD_EXPORT extern int qg_legacy_view_draw_transaction_mode(
	struct ged *gedp,
	const qg_legacy_view_draw_transaction *txn);

QTCAD_EXPORT extern void qg_legacy_view_draw_result_init(
	qg_legacy_view_draw_transaction_result *result);

QTCAD_EXPORT extern void qg_legacy_view_draw_result_free(
	qg_legacy_view_draw_transaction_result *result);

QTCAD_EXPORT extern int qg_legacy_view_draw_transaction_apply(
	struct ged *gedp,
	const qg_legacy_view_draw_transaction *txn,
	qg_legacy_view_draw_transaction_result *result);

QTCAD_EXPORT extern int qg_legacy_view_draw_result_status(
	const qg_legacy_view_draw_transaction_result *result);

QTCAD_EXPORT extern uint64_t qg_legacy_view_draw_result_scene_revision_after(
	const qg_legacy_view_draw_transaction_result *result);

QTCAD_EXPORT extern const char *qg_legacy_view_draw_result_names(
	const qg_legacy_view_draw_transaction_result *result);

QTCAD_EXPORT extern const char *qg_legacy_view_draw_result_errors(
	const qg_legacy_view_draw_transaction_result *result);

QTCAD_EXPORT extern int qg_legacy_view_draw_paths_apply(struct ged *gedp,
	qg_legacy_view *view,
	const char **paths,
	int path_count,
	int autoview);

QTCAD_EXPORT extern int qg_legacy_view_draw_erase_path_apply(
	struct ged *gedp,
	qg_legacy_view *view,
	const char *path);

QTCAD_EXPORT extern int qg_legacy_view_draw_highlight_shape_by_name(
	struct ged *gedp,
	const char *name);

QTCAD_EXPORT extern int qg_legacy_view_draw_highlight_clear(
	struct ged *gedp);

QTCAD_EXPORT extern int qg_legacy_view_refresh_request(qg_legacy_view *view,
	unsigned flags);

QTCAD_EXPORT extern uint32_t qg_legacy_view_refresh_consume(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_refresh_complete(qg_legacy_view *view);

QTCAD_EXPORT extern struct qg_legacy_view_bounds_update_state *
qg_legacy_view_bounds_update_suspend(qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_bounds_update_restore(
	qg_legacy_view *view,
	struct qg_legacy_view_bounds_update_state *state,
	int refresh_bounds);

QTCAD_EXPORT extern int qg_legacy_view_adjust(qg_legacy_view *view,
	int dx,
	int dy,
	point_t keypoint,
	int mode,
	unsigned long long flags);

QTCAD_EXPORT extern double qg_legacy_view_local2base_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern double qg_legacy_view_base2local_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern void qg_legacy_view_info_get(const qg_legacy_view *view,
	struct rt_view_info *info);

QTCAD_EXPORT extern const char *qg_legacy_view_name_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_aet_get(const qg_legacy_view *view,
	vect_t aet);

QTCAD_EXPORT extern int qg_legacy_view_center_get(const qg_legacy_view *view,
	mat_t center);

QTCAD_EXPORT extern int qg_legacy_view_rotation_get(
	const qg_legacy_view *view,
	mat_t rotation);

QTCAD_EXPORT extern fastf_t qg_legacy_view_scale_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_scale_set(qg_legacy_view *view,
	fastf_t scale);

QTCAD_EXPORT extern fastf_t qg_legacy_view_perspective_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_feature_ref_is_null(
	qg_legacy_view_feature_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_edit_preview_publish_event(
	qg_legacy_view *view,
	qg_legacy_view_feature_ref feature,
	qg_legacy_view_edit_preview_event event,
	const char *source_path);

QTCAD_EXPORT extern qg_legacy_view_feature_ref
qg_legacy_view_feature_overlay_ensure(qg_legacy_view *view,
	const char *name,
	const void *owner,
	void *preview_ctx,
	const qg_legacy_view_edit_preview_callbacks *callbacks,
	const char *source_path);

QTCAD_EXPORT extern qg_legacy_view_feature_ref
qg_legacy_view_feature_label_ensure(qg_legacy_view *view,
	const char *name,
	const void *owner);

QTCAD_EXPORT extern int qg_legacy_view_feature_remove(
	qg_legacy_view *view,
	const char *name);

QTCAD_EXPORT extern void qg_legacy_view_feature_set_view(
	qg_legacy_view_feature_ref ref,
	qg_legacy_view *view);

QTCAD_EXPORT extern void qg_legacy_view_feature_set_visible(
	qg_legacy_view_feature_ref ref,
	int visible);

QTCAD_EXPORT extern void qg_legacy_view_feature_set_color(
	qg_legacy_view_feature_ref ref,
	int r,
	int g,
	int b);

QTCAD_EXPORT extern int qg_legacy_view_feature_touch(
	qg_legacy_view_feature_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_feature_labels_replace(
	qg_legacy_view_feature_ref ref,
	const qg_legacy_view_feature_label *labels,
	size_t label_count);

QTCAD_EXPORT extern int qg_legacy_view_feature_points_replace(
	qg_legacy_view_feature_ref ref,
	qg_legacy_view_feature_family family,
	const point_t *points,
	const int *cmds,
	size_t point_count);

QTCAD_EXPORT extern int qg_legacy_view_feature_clear_geometry(
	qg_legacy_view_feature_ref ref);

QTCAD_EXPORT extern unsigned long long qg_legacy_view_hash(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_unique_object_name(struct bu_vls *name,
	const char *seed,
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_framebuffer_mode_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_framebuffer_mode_set(
	qg_legacy_view *view,
	int mode);

QTCAD_EXPORT extern int qg_legacy_view_lod_policy_get(
	const qg_legacy_view *view,
	struct rt_view_lod_policy *policy);

QTCAD_EXPORT extern int qg_legacy_view_lod_policy_apply(qg_legacy_view *view,
	const struct rt_view_lod_policy *policy);

QTCAD_EXPORT extern int qg_legacy_view_lod_policy_copy(qg_legacy_view *dst,
	const qg_legacy_view *src);

QTCAD_EXPORT extern int qg_legacy_view_autoview_default(qg_legacy_view *view,
	int all_view_objs);

QTCAD_EXPORT extern int qg_legacy_view_lod_bounds_update(qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_adc_state_get(
	const qg_legacy_view *view,
	struct rt_view_adc_state *state);

QTCAD_EXPORT extern int qg_legacy_view_adc_state_set(qg_legacy_view *view,
	const struct rt_view_adc_state *state);

QTCAD_EXPORT extern int qg_legacy_view_center_dot_state_get(
	const qg_legacy_view *view,
	struct rt_view_other_state *state);

QTCAD_EXPORT extern int qg_legacy_view_center_dot_state_set(
	qg_legacy_view *view,
	const struct rt_view_other_state *state);

QTCAD_EXPORT extern int qg_legacy_view_grid_state_get(
	const qg_legacy_view *view,
	struct rt_view_grid_state *state);

QTCAD_EXPORT extern int qg_legacy_view_grid_state_set(qg_legacy_view *view,
	const struct rt_view_grid_state *state);

QTCAD_EXPORT extern int qg_legacy_view_model_axes_state_get(
	const qg_legacy_view *view,
	struct rt_view_axes_state *state);

QTCAD_EXPORT extern int qg_legacy_view_model_axes_state_set(
	qg_legacy_view *view,
	const struct rt_view_axes_state *state);

QTCAD_EXPORT extern int qg_legacy_view_scale_overlay_state_get(
	const qg_legacy_view *view,
	struct rt_view_other_state *state);

QTCAD_EXPORT extern int qg_legacy_view_scale_overlay_state_set(
	qg_legacy_view *view,
	const struct rt_view_other_state *state);

QTCAD_EXPORT extern int qg_legacy_view_view_axes_state_get(
	const qg_legacy_view *view,
	struct rt_view_axes_state *state);

QTCAD_EXPORT extern int qg_legacy_view_view_axes_state_set(
	qg_legacy_view *view,
	const struct rt_view_axes_state *state);

QTCAD_EXPORT extern int qg_legacy_view_params_state_get(
	const qg_legacy_view *view,
	struct rt_view_params_state *state);

QTCAD_EXPORT extern int qg_legacy_view_params_state_set(qg_legacy_view *view,
	const struct rt_view_params_state *state);

QTCAD_EXPORT extern int qg_legacy_view_snap_source_view_only_set(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_snap_source_db_set(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_snap_lines_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_snap_lines_set(qg_legacy_view *view,
	int enabled);

QTCAD_EXPORT extern int qg_legacy_view_db_snap_enabled(
	const qg_legacy_view *view);

QTCAD_EXPORT extern unsigned qg_legacy_view_snap_kind_mask_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern fastf_t qg_legacy_view_size_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_size_set(qg_legacy_view *view,
	fastf_t size);

QTCAD_EXPORT extern double qg_legacy_view_snap_tolerance_factor_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_snap_exclude_clear(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_center_vec_set(qg_legacy_view *view,
	const point_t center);

QTCAD_EXPORT extern int qg_legacy_view_aet_set(qg_legacy_view *view,
	const vect_t aet);

QTCAD_EXPORT extern int qg_legacy_view_update(qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_screen_to_view(
	qg_legacy_view *view,
	fastf_t *vx,
	fastf_t *vy,
	fastf_t sx,
	fastf_t sy);

QTCAD_EXPORT extern int qg_legacy_view_screen_point_get(
	qg_legacy_view *view,
	point_t point,
	fastf_t sx,
	fastf_t sy);

QTCAD_EXPORT extern int qg_legacy_view_current_point_get(
	const qg_legacy_view *view,
	point_t point);

QTCAD_EXPORT extern int qg_legacy_view_current_point_set(
	qg_legacy_view *view,
	const point_t point);

QTCAD_EXPORT extern int qg_legacy_view_mouse_state_set(qg_legacy_view *view,
	int x,
	int y);

QTCAD_EXPORT extern int qg_legacy_view_view2model_get(
	const qg_legacy_view *view,
	mat_t view2model);

QTCAD_EXPORT extern int qg_legacy_view_view2model_set(qg_legacy_view *view,
	const mat_t view2model);

QTCAD_EXPORT extern int qg_legacy_view_model2view_get(
	const qg_legacy_view *view,
	mat_t model2view);

QTCAD_EXPORT extern int qg_legacy_view_ray_from_screen(
	qg_legacy_view *view,
	int sx,
	int sy,
	point_t origin,
	vect_t direction);

QTCAD_EXPORT extern int qg_legacy_view_interactive_rect_state_get(
	qg_legacy_view *view,
	struct rt_view_interactive_rect_state *state);

QTCAD_EXPORT extern int qg_legacy_view_interactive_rect_state_set(
	qg_legacy_view *view,
	const struct rt_view_interactive_rect_state *state);

QTCAD_EXPORT extern qg_legacy_view_pick_result *
qg_legacy_view_pick_point(qg_legacy_view *view,
	int x,
	int y,
	int first_only);

QTCAD_EXPORT extern qg_legacy_view_pick_result *
qg_legacy_view_pick_nearest(qg_legacy_view *view,
	int x,
	int y);

QTCAD_EXPORT extern qg_legacy_view_pick_result *
qg_legacy_view_pick_rect(qg_legacy_view *view,
	int x0,
	int y0,
	int x1,
	int y1);

QTCAD_EXPORT extern qg_legacy_view_pick_result *
qg_legacy_view_pick_result_create(void);

QTCAD_EXPORT extern void qg_legacy_view_pick_result_free(
	qg_legacy_view_pick_result *result);

QTCAD_EXPORT extern size_t qg_legacy_view_pick_result_count(
	const qg_legacy_view_pick_result *result);

QTCAD_EXPORT extern int qg_legacy_view_pick_result_path(
	const qg_legacy_view_pick_result *result,
	size_t index,
	struct bu_vls *path_out);

QTCAD_EXPORT extern fastf_t qg_legacy_view_pick_result_hit_dist(
	const qg_legacy_view_pick_result *result,
	size_t index);

QTCAD_EXPORT extern int qg_legacy_view_pick_result_append_copy(
	qg_legacy_view_pick_result *dest,
	const qg_legacy_view_pick_result *src,
	size_t index,
	fastf_t hit_dist);

QTCAD_EXPORT extern qg_legacy_view_pick_result *
qg_legacy_view_pick_result_filter_first(
	const qg_legacy_view_pick_result *src);

QTCAD_EXPORT extern int qg_legacy_view_selection_set_pick_result_ref(
	qg_legacy_view *view,
	const qg_legacy_view_pick_result *result,
	qg_legacy_view_selection_path_callback_t callback,
	void *data);

QTCAD_EXPORT extern int qg_legacy_view_selection_clear(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_polygon_ref_is_null(
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_polygon_record_get(
	rt_view_polygon_ref ref,
	struct rt_view_polygon_record *record);

QTCAD_EXPORT extern rt_view_polygon_ref qg_legacy_view_polygon_create(
	qg_legacy_view *view,
	int type,
	point_t *first_point);

QTCAD_EXPORT extern rt_view_polygon_ref qg_legacy_view_polygon_select(
	qg_legacy_view *view,
	point_t *current_point);

QTCAD_EXPORT extern rt_view_polygon_ref qg_legacy_view_polygon_find(
	qg_legacy_view *view,
	const char *name);

QTCAD_EXPORT extern rt_view_polygon_ref qg_legacy_view_polygon_dup(
	qg_legacy_view *view,
	const char *name,
	const char *new_name);

QTCAD_EXPORT extern void qg_legacy_view_polygon_visit_records(
	qg_legacy_view *view,
	qg_legacy_view_polygon_record_callback_t callback,
	void *data);

QTCAD_EXPORT extern size_t qg_legacy_view_polygon_snap_count(
	qg_legacy_view *view,
	rt_view_polygon_ref exclude);

QTCAD_EXPORT extern int qg_legacy_view_polygon_clear_point_selection(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_polygon_update(
	rt_view_polygon_ref ref,
	qg_legacy_view *view,
	int utype);

QTCAD_EXPORT extern int qg_legacy_view_polygon_update_screen_point(
	rt_view_polygon_ref ref,
	qg_legacy_view *view,
	int x,
	int y,
	int utype);

QTCAD_EXPORT extern int qg_legacy_view_polygon_move(
	rt_view_polygon_ref ref,
	point_t *current_point,
	point_t *previous_point);

QTCAD_EXPORT extern int qg_legacy_view_polygon_set_name(
	rt_view_polygon_ref ref,
	const char *name);

QTCAD_EXPORT extern int qg_legacy_view_polygon_set_view(
	rt_view_polygon_ref ref,
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_polygon_set_visual(
	rt_view_polygon_ref ref,
	const struct bu_color *edge_color,
	const struct bu_color *fill_color,
	fastf_t fill_slope_x,
	fastf_t fill_slope_y,
	fastf_t fill_density,
	fastf_t vZ,
	int fill_flag);

QTCAD_EXPORT extern int qg_legacy_view_polygon_set_open(
	rt_view_polygon_ref ref,
	int open);

QTCAD_EXPORT extern int qg_legacy_view_polygon_close(
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_polygon_clear_selected_point(
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_polygon_remove(
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern void *qg_legacy_view_polygon_user_data(
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_polygon_user_data_set(
	rt_view_polygon_ref ref,
	void *user_data);

QTCAD_EXPORT extern int qg_legacy_view_polygon_csg(
	rt_view_polygon_ref target,
	rt_view_polygon_ref stencil,
	bg_clip_t op);

QTCAD_EXPORT extern rt_view_polygon_ref qg_legacy_view_polygon_import_sketch(
	const char *name,
	struct db_i *dbip,
	struct directory *dp,
	qg_legacy_view *view);

QTCAD_EXPORT extern struct directory *qg_legacy_view_polygon_export_sketch(
	struct db_i *dbip,
	const char *name,
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_polygon_snap_exclude_set(
	qg_legacy_view *view,
	rt_view_polygon_ref ref);

#endif /* QGLEGACYVIEW_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
