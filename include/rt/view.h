/*                          V I E W . H
 * BRL-CAD
 *
 * Copyright (c) 1993-2026 United States Government as represented by
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
/** @file rt/view.h
 *
 */

#ifndef RT_VIEW_H
#define RT_VIEW_H

#include "common.h"
#include "vmath.h"
#include "bg/polygon.h"
#include "bu/color.h"
#include "bu/list.h"
#include "bu/hash.h"
#include "bu/ptbl.h"
#include "bn/tol.h"
#include "rt/defines.h"

__BEGIN_DECLS

struct db_i;
struct directory;
struct bu_vls;

typedef void (*rt_view_context_update_callback_t)(void *view_ctx, void *data);
typedef void (*rt_view_selection_path_callback_t)(const char *path, void *data);

typedef enum rt_view_scene_backend {
    RT_VIEW_SCENE_BACKEND_NONE = 0,
    RT_VIEW_SCENE_BACKEND_BSG = 1,
    RT_VIEW_SCENE_BACKEND_OBOL = 2
} rt_view_scene_backend;

typedef struct rt_view_scene_ref {
    void *opaque;
    unsigned int backend;
} rt_view_scene_ref;

#define RT_VIEW_SCENE_REF_NULL_INIT { NULL, RT_VIEW_SCENE_BACKEND_NONE }

/* Historical normalized view-coordinate range. */
#define RT_VIEW_MAX 2047.0
#define RT_VIEW_MIN -2048.0
#define RT_VIEW_RANGE 4095.0
#define RT_INV_VIEW 0.00048828125
#define RT_INV_4096 0.000244140625
#define RT_VIEW_MIN_SIZE 0.0001
#define RT_VIEW_MIN_SCALE 0.00005
#define RT_VIEW_AUTOVIEW_SCALE_DEFAULT -1
#define RT_VIEW_ADJUST_IDLE   0x000ULL
#define RT_VIEW_ADJUST_ROT    0x001ULL
#define RT_VIEW_ADJUST_TRANS  0x002ULL
#define RT_VIEW_ADJUST_SCALE  0x004ULL
#define RT_VIEW_ADJUST_CENTER 0x008ULL
#define RT_VIEW_ADJUST_CON_X  0x010ULL
#define RT_VIEW_ADJUST_CON_Y  0x020ULL
#define RT_VIEW_ADJUST_CON_Z  0x040ULL
#define RT_VIEW_ADJUST_CON_GRID 0x080ULL
#define RT_VIEW_ADJUST_CON_LINES 0x100ULL
#define RT_VIEW_REFRESH_VIEW        0x00000001u
#define RT_VIEW_REFRESH_DRAW        0x00000002u
#define RT_VIEW_REFRESH_EDIT        0x00000004u
#define RT_VIEW_REFRESH_FRAMEBUFFER 0x00000008u
#define RT_VIEW_REFRESH_OVERLAY     0x00000010u
#define RT_VIEW_REFRESH_FORCE       0x80000000u
#define RT_VIEW_REFRESH_ALL         0xffffffffu
#define RT_VIEW_CLEAR_DB    0x01
#define RT_VIEW_CLEAR_VIEW  0x02
#define RT_VIEW_CLEAR_LOCAL 0x04

struct rt_view_lod_settings {
    fastf_t scale;
    fastf_t curve_scale;
    fastf_t point_scale;
    size_t bot_threshold;
};

enum rt_view_lod_policy_mode {
    RT_VIEW_LOD_OFF = 0,
    RT_VIEW_LOD_AUTO = 1,
    RT_VIEW_LOD_FORCE_LEVEL = 2
};

struct rt_view_lod_policy {
    int policy;
    int forced_level;
    int mesh_enabled;
    int csg_enabled;
    int zoom_refresh;
    size_t bot_threshold;
    fastf_t scale;
    fastf_t curve_scale;
    fastf_t point_scale;
};

struct rt_view_info {
    int width;
    int height;
    fastf_t size;
    struct rt_view_lod_settings lod;
};

struct rt_view_mouse_delta_settings {
    fastf_t min_delta;
    fastf_t max_delta;
    fastf_t rotate_scale;
    fastf_t scale_scale;
};

#define RT_VIEW_MOUSE_DELTA_SETTINGS_INIT { 0.0, 0.0, 0.0, 0.0 }

struct rt_view_knobs {
    vect_t rot_m;
    int rot_m_flag;
    char origin_m;
    void *rot_m_udata;

    vect_t rot_o;
    int rot_o_flag;
    char origin_o;
    void *rot_o_udata;

    vect_t rot_v;
    int rot_v_flag;
    char origin_v;
    void *rot_v_udata;

    fastf_t sca;
    int sca_flag;
    void *sca_udata;

    vect_t tra_m;
    int tra_m_flag;
    void *tra_m_udata;

    vect_t tra_v;
    int tra_v_flag;
    void *tra_v_udata;

    vect_t rot_m_abs;
    vect_t rot_m_abs_last;
    vect_t rot_o_abs;
    vect_t rot_o_abs_last;
    vect_t rot_v_abs;
    vect_t rot_v_abs_last;

    fastf_t sca_abs;

    vect_t tra_m_abs;
    vect_t tra_m_abs_last;
    vect_t tra_v_abs;
    vect_t tra_v_abs_last;
};

enum rt_view_knobs_category {
    RT_VIEW_KNOBS_ALL = 0,
    RT_VIEW_KNOBS_RATE = 1,
    RT_VIEW_KNOBS_ABS = 2
};

#define RT_VIEW_SNAP_SHARED 0x1
#define RT_VIEW_SNAP_LOCAL  0x2
#define RT_VIEW_SNAP_DB     0x4
#define RT_VIEW_SNAP_VIEW   0x8
#define RT_VIEW_SNAP_TCL    0x10

#define RT_VIEW_SNAP_KIND_GRID           0x01ULL
#define RT_VIEW_SNAP_KIND_ENDPOINT       0x02ULL
#define RT_VIEW_SNAP_KIND_MIDPOINT       0x04ULL
#define RT_VIEW_SNAP_KIND_INTERSECTION   0x08ULL
#define RT_VIEW_SNAP_KIND_PERPENDICULAR  0x10ULL
#define RT_VIEW_SNAP_KIND_TANGENT        0x20ULL
#define RT_VIEW_SNAP_KIND_OVERLAY_HANDLE 0x40ULL

struct rt_view_knob_values {
    vect_t rate_rotation;
    vect_t rate_translation;
    fastf_t rate_scale;
    vect_t absolute_rotation;
    vect_t absolute_translation;
    fastf_t absolute_scale;
};

RT_EXPORT extern int
rt_view_knobs_state_reset(struct rt_view_knobs *knobs, int category);

RT_EXPORT extern unsigned long long
rt_view_knobs_state_hash(const struct rt_view_knobs *knobs,
			 struct bu_data_hash_state *state);

struct rt_view_interactive_rect_state {
    int active;
    int draw;
    int line_width;
    int line_style;
    int pos[2];
    int dim[2];
    fastf_t x;
    fastf_t y;
    fastf_t width;
    fastf_t height;
    int bg[3];
    int color[3];
    int cdim[2];
    fastf_t aspect;
};

struct rt_view_measure_result {
    fastf_t distance;
    fastf_t projection;
    fastf_t normal_alignment;
    int valid;
};

struct rt_view_adc_state {
    int         draw;
    int         dv_x;
    int         dv_y;
    int         dv_a1;
    int         dv_a2;
    int         dv_dist;
    fastf_t     pos_model[3];
    fastf_t     pos_view[3];
    fastf_t     pos_grid[3];
    fastf_t     a1;
    fastf_t     a2;
    fastf_t     dst;
    int         anchor_pos;
    int         anchor_a1;
    int         anchor_a2;
    int         anchor_dst;
    fastf_t     anchor_pt_a1[3];
    fastf_t     anchor_pt_a2[3];
    fastf_t     anchor_pt_dst[3];
    int         line_color[3];
    int         tick_color[3];
    int         line_width;
};

struct rt_view_grid_state {
    int       rc;
    int       draw;
    int       adaptive;
    int       snap;
    fastf_t   anchor[3];
    fastf_t   res_h;
    fastf_t   res_v;
    int       res_major_h;
    int       res_major_v;
    int       color[3];
};

struct rt_view_axes_state {
    int       draw;
    point_t   axes_pos;
    fastf_t   axes_size;
    int       line_width;
    int       axes_color[3];
    int       pos_only;
    int       label_flag;
    int       label_color[3];
    int       triple_color;
    int       tick_enabled;
    int       tick_length;
    int       tick_major_length;
    fastf_t   tick_interval;
    int       ticks_per_major;
    int       tick_threshold;
    int       tick_color[3];
    int       tick_major_color[3];
};

struct rt_view_other_state {
    int gos_draw;
    int gos_line_color[3];
    int gos_text_color[3];
    int gos_font_size;
};

struct rt_view_params_state {
    int draw;
    int draw_size;
    int draw_center;
    int draw_az;
    int draw_el;
    int draw_tw;
    int draw_fps;
    int color[3];
    int font_size;
};

typedef struct rt_view_feature_ref {
    uintptr_t token;
    uint64_t revision;
} rt_view_feature_ref;

#define RT_VIEW_FEATURE_REF_NULL_INIT {0, 0}
#ifdef __cplusplus
#  define RT_VIEW_FEATURE_REF_NULL rt_view_feature_ref{0, 0}
#else
#  define RT_VIEW_FEATURE_REF_NULL ((rt_view_feature_ref){0, 0})
#endif

enum rt_view_feature_family {
    RT_VIEW_FEATURE_UNKNOWN = 0,
    RT_VIEW_FEATURE_TRANSIENT_PREVIEW = 1
};

enum rt_view_edit_preview_event {
    RT_VIEW_EDIT_PREVIEW_BEGIN = 0,
    RT_VIEW_EDIT_PREVIEW_UPDATE,
    RT_VIEW_EDIT_PREVIEW_COMMIT,
    RT_VIEW_EDIT_PREVIEW_CANCEL,
    RT_VIEW_EDIT_PREVIEW_REPLACE_SOURCE,
    RT_VIEW_EDIT_PREVIEW_DISCARD
};

struct rt_view_edit_preview_callbacks {
    uint64_t (*revision_cb)(void *);
    int (*update_cb)(void *);
    int (*pick_cb)(void *, int, int, void *);
};

#define RT_VIEW_EDIT_PREVIEW_CALLBACKS_INIT { 0, 0, 0 }

struct rt_view_feature_label {
    const char *text;
    point_t point;
    int color_valid;
    unsigned char color[3];
    fastf_t font_size;
};

typedef int (*rt_view_context_feature_owns_ref_callback_t)(
    rt_view_feature_ref ref,
    void *data);
typedef int (*rt_view_context_edit_preview_publish_event_callback_t)(
    void *ctx,
    rt_view_feature_ref feature,
    enum rt_view_edit_preview_event event,
    const char *source_path,
    void *data);
typedef rt_view_feature_ref(*rt_view_context_feature_overlay_ensure_callback_t)(
    void *ctx,
    const char *name,
    const void *owner,
    void *preview_ctx,
    const struct rt_view_edit_preview_callbacks *callbacks,
    const char *source_path,
    void *data);
typedef rt_view_feature_ref(*rt_view_context_feature_label_ensure_callback_t)(
    void *ctx,
    const char *name,
    const void *owner,
    void *data);
typedef int (*rt_view_context_feature_remove_callback_t)(
    void *ctx,
    const char *name,
    void *data);
typedef int (*rt_view_feature_set_context_callback_t)(
    rt_view_feature_ref ref,
    void *ctx,
    void *data);
typedef int (*rt_view_feature_set_visible_callback_t)(
    rt_view_feature_ref ref,
    int visible,
    void *data);
typedef int (*rt_view_feature_set_color_callback_t)(
    rt_view_feature_ref ref,
    int r,
    int g,
    int b,
    void *data);
typedef int (*rt_view_feature_touch_callback_t)(
    rt_view_feature_ref ref,
    void *data);
typedef int (*rt_view_feature_labels_replace_callback_t)(
    rt_view_feature_ref ref,
    const struct rt_view_feature_label *labels,
    size_t label_count,
    void *data);
typedef int (*rt_view_feature_points_replace_callback_t)(
    rt_view_feature_ref ref,
    enum rt_view_feature_family family,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    void *data);
typedef int (*rt_view_feature_clear_geometry_callback_t)(
    rt_view_feature_ref ref,
    void *data);

struct rt_view_context_feature_adapter {
    rt_view_context_feature_owns_ref_callback_t owns_ref;
    rt_view_context_edit_preview_publish_event_callback_t edit_preview_publish_event;
    rt_view_context_feature_overlay_ensure_callback_t overlay_ensure;
    rt_view_context_feature_label_ensure_callback_t label_ensure;
    rt_view_context_feature_remove_callback_t remove;
    rt_view_feature_set_context_callback_t set_context;
    rt_view_feature_set_visible_callback_t set_visible;
    rt_view_feature_set_color_callback_t set_color;
    rt_view_feature_touch_callback_t touch;
    rt_view_feature_labels_replace_callback_t labels_replace;
    rt_view_feature_points_replace_callback_t points_replace;
    rt_view_feature_clear_geometry_callback_t clear_geometry;
    void *data;
};

RT_EXPORT extern int
rt_view_context_feature_adapter_set(
    void *ctx,
    const struct rt_view_context_feature_adapter *adapter);

RT_EXPORT extern int
rt_view_context_feature_adapter_get(
    void *ctx,
    struct rt_view_context_feature_adapter *adapter);

RT_EXPORT extern int
rt_view_feature_ref_is_null(rt_view_feature_ref ref);

RT_EXPORT extern int
rt_view_context_edit_preview_publish_event(
    void *ctx,
    rt_view_feature_ref feature,
    enum rt_view_edit_preview_event event,
    const char *source_path);

RT_EXPORT extern rt_view_feature_ref
rt_view_context_feature_overlay_ensure(
    void *ctx,
    const char *name,
    const void *owner,
    void *preview_ctx,
    const struct rt_view_edit_preview_callbacks *callbacks,
    const char *source_path);

RT_EXPORT extern rt_view_feature_ref
rt_view_context_feature_label_ensure(void *ctx,
				     const char *name,
				     const void *owner);

RT_EXPORT extern int
rt_view_context_feature_remove(void *ctx, const char *name);

RT_EXPORT extern void
rt_view_feature_set_context(rt_view_feature_ref ref, void *ctx);

RT_EXPORT extern void
rt_view_feature_set_visible(rt_view_feature_ref ref, int visible);

RT_EXPORT extern void
rt_view_feature_set_color(rt_view_feature_ref ref, int r, int g, int b);

RT_EXPORT extern int
rt_view_feature_touch(rt_view_feature_ref ref);

RT_EXPORT extern int
rt_view_feature_labels_replace(rt_view_feature_ref ref,
			       const struct rt_view_feature_label *labels,
			       size_t label_count);

RT_EXPORT extern int
rt_view_feature_points_replace(rt_view_feature_ref ref,
			       enum rt_view_feature_family family,
			       const point_t *points,
			       const int *cmds,
			       size_t point_count);

RT_EXPORT extern int
rt_view_feature_clear_geometry(rt_view_feature_ref ref);

typedef struct rt_view_polygon_ref {
    uintptr_t token;
    uint64_t revision;
} rt_view_polygon_ref;

#define RT_VIEW_POLYGON_GENERAL 0

#define RT_VIEW_POLYGON_REF_NULL_INIT {0, 0}
#ifdef __cplusplus
#  define RT_VIEW_POLYGON_REF_NULL rt_view_polygon_ref{0, 0}
#else
#  define RT_VIEW_POLYGON_REF_NULL ((rt_view_polygon_ref){0, 0})
#endif

struct rt_view_polygon_record {
    rt_view_polygon_ref ref;
    const char *name;
    int type;
    int fill_flag;
    vect2d_t fill_dir;
    fastf_t fill_delta;
    struct bu_color fill_color;
    unsigned char edge_color[3];
    long curr_contour_i;
    long curr_point_i;
    int first_contour_open;
    size_t contour_count;
    size_t point_count;
    point_t origin_point;
    plane_t vp;
    fastf_t vZ;
    void *user_data;
};

typedef int (*rt_view_polygon_record_callback_t)(rt_view_polygon_ref ref,
	const struct rt_view_polygon_record *record,
	void *data);

typedef int (*rt_view_context_polygon_owns_ref_callback_t)(
    rt_view_polygon_ref ref,
    void *data);
typedef int (*rt_view_context_polygon_record_get_callback_t)(
    rt_view_polygon_ref ref,
    struct rt_view_polygon_record *record,
    void *data);
typedef rt_view_polygon_ref(*rt_view_context_polygon_create_callback_t)(
    void *ctx,
    int type,
    point_t *fp,
    void *data);
typedef rt_view_polygon_ref(*rt_view_context_polygon_select_callback_t)(
    void *ctx,
    point_t *cp,
    void *data);
typedef rt_view_polygon_ref(*rt_view_context_polygon_find_callback_t)(
    void *ctx,
    const char *name,
    void *data);
typedef rt_view_polygon_ref(*rt_view_context_polygon_dup_callback_t)(
    void *ctx,
    const char *name,
    const char *new_name,
    void *data);
typedef void (*rt_view_context_polygon_visit_records_callback_t)(
    void *ctx,
    rt_view_polygon_record_callback_t callback,
    void *callback_data,
    void *data);
typedef size_t (*rt_view_context_polygon_snap_count_callback_t)(
    void *ctx,
    rt_view_polygon_ref exclude,
    void *data);
typedef int (*rt_view_context_polygon_clear_point_selection_callback_t)(
    void *ctx,
    void *data);
typedef int (*rt_view_context_polygon_update_callback_t)(
    rt_view_polygon_ref ref,
    void *ctx,
    int utype,
    void *data);
typedef int (*rt_view_context_polygon_update_screen_pt_callback_t)(
    rt_view_polygon_ref ref,
    void *ctx,
    int x,
    int y,
    int utype,
    void *data);
typedef int (*rt_view_polygon_move_callback_t)(
    rt_view_polygon_ref ref,
    point_t *current_point,
    point_t *previous_point,
    void *data);
typedef int (*rt_view_polygon_set_name_callback_t)(
    rt_view_polygon_ref ref,
    const char *name,
    void *data);
typedef int (*rt_view_polygon_set_context_callback_t)(
    rt_view_polygon_ref ref,
    void *ctx,
    void *data);
typedef int (*rt_view_polygon_set_visual_callback_t)(
    rt_view_polygon_ref ref,
    const struct bu_color *edge_color,
    const struct bu_color *fill_color,
    fastf_t fill_slope_x,
    fastf_t fill_slope_y,
    fastf_t fill_density,
    fastf_t vZ,
    int fill_flag,
    void *data);
typedef int (*rt_view_polygon_set_open_callback_t)(
    rt_view_polygon_ref ref,
    int open,
    void *data);
typedef int (*rt_view_polygon_close_callback_t)(
    rt_view_polygon_ref ref,
    void *data);
typedef int (*rt_view_polygon_clear_selected_point_callback_t)(
    rt_view_polygon_ref ref,
    void *data);
typedef int (*rt_view_polygon_remove_callback_t)(
    rt_view_polygon_ref ref,
    void *data);
typedef void *(*rt_view_polygon_user_data_callback_t)(
    rt_view_polygon_ref ref,
    void *data);
typedef int (*rt_view_polygon_user_data_set_callback_t)(
    rt_view_polygon_ref ref,
    void *user_data,
    void *data);
typedef int (*rt_view_polygon_csg_callback_t)(
    rt_view_polygon_ref target,
    rt_view_polygon_ref stencil,
    bg_clip_t op,
    void *data);
typedef rt_view_polygon_ref(*rt_view_polygon_import_sketch_context_callback_t)(
    const char *name,
    struct db_i *dbip,
    struct directory *dp,
    void *ctx,
    void *data);
typedef struct directory *(*rt_view_polygon_export_sketch_callback_t)(
    struct db_i *dbip,
    const char *name,
    rt_view_polygon_ref ref,
    void *data);
typedef int (*rt_view_context_polygon_snap_exclude_set_callback_t)(
    void *ctx,
    rt_view_polygon_ref ref,
    void *data);

struct rt_view_context_polygon_adapter {
    rt_view_context_polygon_owns_ref_callback_t owns_ref;
    rt_view_context_polygon_record_get_callback_t record_get;
    rt_view_context_polygon_create_callback_t create;
    rt_view_context_polygon_select_callback_t select;
    rt_view_context_polygon_find_callback_t find;
    rt_view_context_polygon_dup_callback_t dup;
    rt_view_context_polygon_visit_records_callback_t visit_records;
    rt_view_context_polygon_snap_count_callback_t snap_count;
    rt_view_context_polygon_clear_point_selection_callback_t clear_point_selection;
    rt_view_context_polygon_update_callback_t update;
    rt_view_context_polygon_update_screen_pt_callback_t update_screen_pt;
    rt_view_polygon_move_callback_t move;
    rt_view_polygon_set_name_callback_t set_name;
    rt_view_polygon_set_context_callback_t set_context;
    rt_view_polygon_set_visual_callback_t set_visual;
    rt_view_polygon_set_open_callback_t set_open;
    rt_view_polygon_close_callback_t close;
    rt_view_polygon_clear_selected_point_callback_t clear_selected_point;
    rt_view_polygon_remove_callback_t remove;
    rt_view_polygon_user_data_callback_t user_data;
    rt_view_polygon_user_data_set_callback_t user_data_set;
    rt_view_polygon_csg_callback_t csg;
    rt_view_polygon_import_sketch_context_callback_t import_sketch_context;
    rt_view_polygon_export_sketch_callback_t export_sketch;
    rt_view_context_polygon_snap_exclude_set_callback_t snap_exclude_set;
    void *data;
};

RT_EXPORT extern int
rt_view_context_polygon_adapter_set(
    void *ctx,
    const struct rt_view_context_polygon_adapter *adapter);

RT_EXPORT extern int
rt_view_context_polygon_adapter_get(
    void *ctx,
    struct rt_view_context_polygon_adapter *adapter);

RT_EXPORT extern int
rt_view_polygon_ref_is_null(rt_view_polygon_ref ref);

RT_EXPORT extern int
rt_view_polygon_record_get(rt_view_polygon_ref ref,
			   struct rt_view_polygon_record *record);

RT_EXPORT extern rt_view_polygon_ref
rt_view_context_polygon_create(void *ctx, int type, point_t *fp);

RT_EXPORT extern rt_view_polygon_ref
rt_view_context_polygon_select(void *ctx, point_t *cp);

RT_EXPORT extern rt_view_polygon_ref
rt_view_context_polygon_find(void *ctx, const char *name);

RT_EXPORT extern rt_view_polygon_ref
rt_view_context_polygon_dup(void *ctx,
			    const char *name,
			    const char *new_name);

RT_EXPORT extern void
rt_view_context_polygon_visit_records(
    void *ctx,
    rt_view_polygon_record_callback_t callback,
    void *data);

RT_EXPORT extern size_t
rt_view_context_polygon_snap_count(void *ctx, rt_view_polygon_ref exclude);

RT_EXPORT extern int
rt_view_context_polygon_clear_point_selection(void *ctx);

RT_EXPORT extern int
rt_view_polygon_update_context(rt_view_polygon_ref ref, void *ctx, int utype);

RT_EXPORT extern int
rt_view_polygon_update_screen_pt_context(rt_view_polygon_ref ref,
	void *ctx,
	int x,
	int y,
	int utype);

RT_EXPORT extern int
rt_view_polygon_move(rt_view_polygon_ref ref,
		     point_t *current_point,
		     point_t *previous_point);

RT_EXPORT extern int
rt_view_polygon_set_name(rt_view_polygon_ref ref, const char *name);

RT_EXPORT extern int
rt_view_polygon_set_context(rt_view_polygon_ref ref, void *ctx);

RT_EXPORT extern int
rt_view_polygon_set_visual(rt_view_polygon_ref ref,
			   const struct bu_color *edge_color,
			   const struct bu_color *fill_color,
			   fastf_t fill_slope_x,
			   fastf_t fill_slope_y,
			   fastf_t fill_density,
			   fastf_t vZ,
			   int fill_flag);

RT_EXPORT extern int
rt_view_polygon_set_open(rt_view_polygon_ref ref, int open);

RT_EXPORT extern int
rt_view_polygon_close(rt_view_polygon_ref ref);

RT_EXPORT extern int
rt_view_polygon_clear_selected_point(rt_view_polygon_ref ref);

RT_EXPORT extern int
rt_view_polygon_remove(rt_view_polygon_ref ref);

RT_EXPORT extern void *
rt_view_polygon_user_data(rt_view_polygon_ref ref);

RT_EXPORT extern int
rt_view_polygon_user_data_set(rt_view_polygon_ref ref, void *user_data);

RT_EXPORT extern int
rt_view_polygon_csg(rt_view_polygon_ref target,
		    rt_view_polygon_ref stencil,
		    bg_clip_t op);

RT_EXPORT extern rt_view_polygon_ref
rt_view_polygon_import_sketch_context(const char *name,
				      struct db_i *dbip,
				      struct directory *dp,
				      void *ctx);

RT_EXPORT extern struct directory *
rt_view_polygon_export_sketch(struct db_i *dbip,
			      const char *name,
			      rt_view_polygon_ref ref);

RT_EXPORT extern int
rt_view_context_polygon_snap_exclude_set(void *ctx, rt_view_polygon_ref ref);

struct rt_view_render_summary {
    int item_count;
    int highlighted_count;
};

struct rt_view_render_export_consistency {
    int export_record_found;
    int render_item_found;
    int backend_node_found;
    int export_render_consistent;
    int export_backend_consistent;
};

typedef void *(*rt_view_context_pick_semantic_path_callback_t)(
    void *ctx,
    const char *path_pattern,
    void *data);

typedef int (*rt_view_context_render_export_consistency_callback_t)(
    void *ctx,
    const char *drawn_prefix,
    struct rt_view_render_export_consistency *summary,
    void *data);

struct rt_view_context_scene_adapter {
    rt_view_context_pick_semantic_path_callback_t pick_semantic_path;
    rt_view_context_render_export_consistency_callback_t render_export_consistency;
    void *data;
};

typedef int (*rt_view_context_selection_available_callback_t)(
    void *ctx,
    void *data);
typedef size_t (*rt_view_context_selection_count_callback_t)(
    void *ctx,
    void *data);
typedef int (*rt_view_context_selection_set_pick_result_context_callback_t)(
    void *ctx,
    const void *result_ctx,
    rt_view_selection_path_callback_t callback,
    void *callback_data,
    void *data);
typedef int (*rt_view_context_selection_clear_callback_t)(
    void *ctx,
    void *data);

struct rt_view_context_selection_adapter {
    rt_view_context_selection_available_callback_t available;
    rt_view_context_selection_count_callback_t count;
    rt_view_context_selection_set_pick_result_context_callback_t set_pick_result_context;
    rt_view_context_selection_clear_callback_t clear;
    void *data;
};

struct rt_view_feature_geometry_summary {
    int exists;
    size_t point_count;
    size_t command_count;
};

#define RT_VIEW_LOD_SETTINGS_INIT { 1.0, 1.0, 1.0, 0 }
#define RT_VIEW_LOD_POLICY_INIT { RT_VIEW_LOD_AUTO, 0, 0, 0, 0, 0, 1.0, 1.0, 1.0 }
#define RT_VIEW_INFO_INIT { 1, 1, 1.0, RT_VIEW_LOD_SETTINGS_INIT }
#define RT_VIEW_INTERACTIVE_RECT_STATE_INIT { 0, 0, 0, 0, {0, 0}, {0, 0}, 0.0, 0.0, 0.0, 0.0, {0, 0, 0}, {0, 0, 0}, {0, 0}, 0.0 }
#define RT_VIEW_MEASURE_RESULT_INIT { 0.0, 0.0, 0.0, 0 }
#define RT_VIEW_ADC_STATE_INIT { 0, 0, 0, 0, 0, 0, VINIT_ZERO, VINIT_ZERO, VINIT_ZERO, 0.0, 0.0, 0.0, 0, 0, 0, 0, VINIT_ZERO, VINIT_ZERO, VINIT_ZERO, {0, 0, 0}, {0, 0, 0}, 0 }
#define RT_VIEW_GRID_STATE_INIT { 0, 0, 0, 0, VINIT_ZERO, 0.0, 0.0, 0, 0, {0, 0, 0} }
#define RT_VIEW_AXES_STATE_INIT { 0, VINIT_ZERO, 0.0, 0, {0, 0, 0}, 0, 0, {0, 0, 0}, 0, 0, 0, 0, 0.0, 0, 0, {0, 0, 0}, {0, 0, 0} }
#define RT_VIEW_OTHER_STATE_INIT { 0, {0, 0, 0}, {0, 0, 0}, 0 }
#define RT_VIEW_PARAMS_STATE_INIT { 0, 0, 0, 0, 0, 0, 0, {0, 0, 0}, 0 }
#define RT_VIEW_RENDER_SUMMARY_INIT { 0, 0 }
#define RT_VIEW_RENDER_EXPORT_CONSISTENCY_INIT { 0, 0, 0, 0, 0 }
#define RT_VIEW_FEATURE_GEOMETRY_SUMMARY_INIT { 0, 0, 0 }

RT_EXPORT extern void rt_view_info_init(struct rt_view_info *info);
RT_EXPORT extern void rt_view_info_sanitize(struct rt_view_info *info);
RT_EXPORT extern void rt_view_lod_policy_init(struct rt_view_lod_policy *policy);
RT_EXPORT extern void rt_view_lod_policy_sanitize(struct rt_view_lod_policy *policy);
RT_EXPORT extern fastf_t rt_view_lod_curve_scale(const struct rt_view_info *info);
RT_EXPORT extern size_t rt_view_lod_bot_threshold(const struct rt_view_info *info);
RT_EXPORT extern fastf_t rt_view_avg_sample_spacing(const struct rt_view_info *info);
RT_EXPORT extern fastf_t rt_view_solid_point_spacing(const struct rt_view_info *info, fastf_t solid_width);
RT_EXPORT extern int rt_view_context_is_valid(const void *ctx);
RT_EXPORT extern int rt_view_context_is_retained(const void *ctx);
RT_EXPORT extern const char *rt_view_context_name_get(const void *ctx);
RT_EXPORT extern int rt_view_context_name_set(void *ctx, const char *name);
RT_EXPORT extern int rt_view_context_width_get(const void *ctx);
RT_EXPORT extern int rt_view_context_height_get(const void *ctx);
RT_EXPORT extern int rt_view_context_dimensions_set(void *ctx, int width, int height);
RT_EXPORT extern int rt_view_context_unit_conversion_set(void *ctx, fastf_t local2base, fastf_t base2local);
RT_EXPORT extern int rt_view_context_frametime_set(void *ctx, uint64_t frametime);
RT_EXPORT extern uint64_t rt_view_context_frametime_get(const void *ctx);
RT_EXPORT extern fastf_t rt_view_context_local2base_get(const void *ctx);
RT_EXPORT extern fastf_t rt_view_context_base2local_get(const void *ctx);
RT_EXPORT extern fastf_t rt_view_context_scale_get(const void *ctx);
RT_EXPORT extern int rt_view_context_scale_set(void *ctx, fastf_t scale);
RT_EXPORT extern fastf_t *rt_view_context_scale_storage_get(void *ctx);
RT_EXPORT extern int rt_view_context_scale_state_set(void *ctx, fastf_t scale, fastf_t initial_scale, fastf_t absolute_scale, fastf_t size, fastf_t inverse_size);
RT_EXPORT extern fastf_t rt_view_context_initial_scale_get(const void *ctx);
RT_EXPORT extern int rt_view_context_initial_scale_set(void *ctx, fastf_t scale);
RT_EXPORT extern fastf_t rt_view_context_absolute_scale_get(const void *ctx);
RT_EXPORT extern int rt_view_context_absolute_scale_set(void *ctx, fastf_t scale);
RT_EXPORT extern fastf_t rt_view_context_size_get(const void *ctx);
RT_EXPORT extern int rt_view_context_size_set(void *ctx, fastf_t size);
RT_EXPORT extern fastf_t rt_view_context_inverse_size_get(const void *ctx);
RT_EXPORT extern fastf_t rt_view_context_perspective_get(const void *ctx);
RT_EXPORT extern int rt_view_context_perspective_set(void *ctx, fastf_t perspective);
RT_EXPORT extern fastf_t rt_view_context_radius_get(const void *ctx);
RT_EXPORT extern int rt_view_context_eye_pos_get(point_t eye_pos, const void *ctx);
RT_EXPORT extern int rt_view_context_eye_pos_set(void *ctx, const point_t eye_pos);
RT_EXPORT extern int rt_view_context_keypoint_get(point_t keypoint, const void *ctx);
RT_EXPORT extern int rt_view_context_keypoint_set(void *ctx, const point_t keypoint);
RT_EXPORT extern char rt_view_context_rotate_about_get(const void *ctx);
RT_EXPORT extern int rt_view_context_rotate_about_set(void *ctx, char rotate_about);
RT_EXPORT extern char rt_view_context_coord_get(const void *ctx);
RT_EXPORT extern int rt_view_context_coord_set(void *ctx, char coord);
RT_EXPORT extern int rt_view_context_zclip_get(const void *ctx);
RT_EXPORT extern int rt_view_context_zclip_set(void *ctx, int zclip);
RT_EXPORT extern int rt_view_context_model2view_get(mat_t model2view, const void *ctx);
RT_EXPORT extern int rt_view_context_model2view_set(void *ctx, const mat_t model2view);
RT_EXPORT extern int rt_view_context_view2model_get(mat_t view2model, const void *ctx);
RT_EXPORT extern int rt_view_context_view2model_set(void *ctx, const mat_t view2model);
RT_EXPORT extern int rt_view_context_pmodel2view_get(mat_t pmodel2view, const void *ctx);
RT_EXPORT extern int rt_view_context_pmodel2view_set(void *ctx, const mat_t pmodel2view);
RT_EXPORT extern int rt_view_context_pmat_get(mat_t pmat, const void *ctx);
RT_EXPORT extern int rt_view_context_pmat_set(void *ctx, const mat_t pmat);
RT_EXPORT extern int rt_view_context_center_get(mat_t center, const void *ctx);
RT_EXPORT extern int rt_view_context_center_set(void *ctx, const point_t center);
RT_EXPORT extern int rt_view_context_rotation_get(mat_t rotation, const void *ctx);
RT_EXPORT extern int rt_view_context_rotation_set(void *ctx, const mat_t rotation);
RT_EXPORT extern void rt_view_context_info_get(struct rt_view_info *info, const void *ctx);
RT_EXPORT extern int rt_view_context_obb_get(const void *ctx, point_t center, vect_t extent1, vect_t extent2, vect_t extent3);
RT_EXPORT extern int rt_view_context_orientation_quat_get(quat_t orientation, const void *ctx);
RT_EXPORT extern int rt_view_context_aet_get(vect_t aet, const void *ctx);
RT_EXPORT extern int rt_view_context_aet_set(void *ctx, const vect_t aet);
RT_EXPORT extern int rt_view_context_aet_state_set(void *ctx, const vect_t aet);
RT_EXPORT extern int rt_view_context_plane_get(plane_t *plane, const void *ctx);
RT_EXPORT extern int rt_view_context_refresh_request(void *ctx, uint32_t flags);
RT_EXPORT extern int rt_view_context_refresh_dirty_get(const void *ctx);
RT_EXPORT extern uint32_t rt_view_context_refresh_consume(void *ctx);
RT_EXPORT extern int rt_view_context_refresh_complete(void *ctx);
RT_EXPORT extern int rt_view_context_refresh_enabled_get(const void *ctx);
RT_EXPORT extern int rt_view_context_refresh_enabled_set(void *ctx, int enabled);
RT_EXPORT extern int rt_view_context_refresh_suppressed_get(const void *ctx);
RT_EXPORT extern int rt_view_context_refresh_suppress_begin(void *ctx);
RT_EXPORT extern int rt_view_context_refresh_suppress_end(void *ctx);
RT_EXPORT extern int rt_view_context_refresh_drawn_count_get(const void *ctx);
RT_EXPORT extern int rt_view_context_refresh_drawn_count_set(void *ctx, int count);
RT_EXPORT extern size_t rt_view_context_clear(void *ctx, int flags);
RT_EXPORT extern int rt_view_context_cleared_get(const void *ctx);
RT_EXPORT extern int rt_view_context_cleared_set(void *ctx, int cleared);
RT_EXPORT extern void *rt_view_context_user_data_get(const void *ctx);
RT_EXPORT extern int rt_view_context_user_data_set(void *ctx, void *user_data);
RT_EXPORT extern int rt_view_context_tclcad_data_set(void *ctx, void *tcl_data);
RT_EXPORT extern int rt_view_context_callbacks_set(void *ctx, struct bu_ptbl *callbacks);
RT_EXPORT extern int rt_view_context_update_callback_set(void *ctx, rt_view_context_update_callback_t callback, void *data);
RT_EXPORT extern void *rt_view_context_create(void);
RT_EXPORT extern void *rt_view_context_create_with_set(void *view_set_ctx);
RT_EXPORT extern void *rt_view_context_create_copy_with_set(const void *src_ctx, void *view_set_ctx);
RT_EXPORT extern void rt_view_context_free(void *ctx);
RT_EXPORT extern int rt_view_context_release_storage(void *ctx);
RT_EXPORT extern int rt_view_context_view_set_attach(void *ctx, void *view_set_ctx);
RT_EXPORT extern void *rt_view_context_line_set_create(void *ctx, const char *name, unsigned char r, unsigned char g, unsigned char b);
RT_EXPORT extern int rt_view_line_set_context_is_null(const void *line_set_ctx);
RT_EXPORT extern int rt_view_line_set_context_set_points(void *line_set_ctx, const point_t *points, const int *commands, size_t point_count);
RT_EXPORT extern void rt_view_line_set_context_destroy(void *line_set_ctx);
RT_EXPORT extern void *rt_view_context_display_manager_get(const void *ctx);
RT_EXPORT extern int rt_view_context_display_manager_set(void *ctx, void *dmp);
RT_EXPORT extern int rt_view_context_update(void *ctx);
RT_EXPORT extern int rt_view_context_autoview(void *ctx, fastf_t scale, int all_view_objs);
RT_EXPORT extern int rt_view_context_autoview_bounds(void *ctx, fastf_t scale, const point_t min, const point_t max);
RT_EXPORT extern int rt_view_context_adjust(void *ctx, int dx, int dy, point_t keypoint, int mode, unsigned long long flags);
RT_EXPORT extern unsigned long long rt_view_context_hash(const void *ctx);
RT_EXPORT extern int rt_view_context_unique_object_name(struct bu_vls *oname, const char *seed, void *ctx);
RT_EXPORT extern int rt_view_context_edit_matrix_set(void *ctx, matp_t edit_mat);
RT_EXPORT extern int rt_view_context_edit_matrix_clear(void *ctx);
RT_EXPORT extern uint64_t rt_view_context_frame_revision_get(const void *ctx);
RT_EXPORT extern uint64_t rt_view_context_frame_revision_bump(void *ctx);
RT_EXPORT extern int rt_view_context_knobs_state_get(struct rt_view_knobs *knobs, const void *ctx);
RT_EXPORT extern int rt_view_context_knobs_state_set(void *ctx, const struct rt_view_knobs *knobs);
RT_EXPORT extern unsigned long long rt_view_context_knobs_hash(void *ctx, struct bu_data_hash_state *state);
RT_EXPORT extern int rt_view_context_knob_values_get(struct rt_view_knob_values *values, const void *ctx);
RT_EXPORT extern int rt_view_context_knobs_reset(void *ctx, int category);
RT_EXPORT extern int rt_view_context_knobs_calibrate(void *ctx);
RT_EXPORT extern int rt_view_context_knobs_cmd_process(vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, void *ctx, const char *cmd, fastf_t factor, char origin, int model_flag, int incr_flag);
RT_EXPORT extern int rt_view_context_knobs_translate(void *ctx, const vect_t tvec, int model_flag);
RT_EXPORT extern int rt_view_context_knobs_rotate(void *ctx, const vect_t rvec, char origin, char coords, const matp_t obj_rot, const pointp_t pvt_pt);
RT_EXPORT extern int rt_view_context_knobs_update_rate_flags(void *ctx);
RT_EXPORT extern int rt_view_context_screen_to_view(fastf_t *fx, fastf_t *fy, void *ctx, fastf_t x, fastf_t y);
RT_EXPORT extern int rt_view_context_screen_point_get(point_t point, void *ctx, fastf_t x, fastf_t y);
RT_EXPORT extern int rt_view_context_current_point_get(point_t point, const void *ctx);
RT_EXPORT extern int rt_view_context_current_point_set(void *ctx, const point_t point);
RT_EXPORT extern int rt_view_context_previous_mouse_get(fastf_t *x, fastf_t *y, const void *ctx);
RT_EXPORT extern int rt_view_context_previous_mouse_set(void *ctx, fastf_t x, fastf_t y);
RT_EXPORT extern int rt_view_context_mouse_delta_settings_get(struct rt_view_mouse_delta_settings *settings, const void *ctx);
RT_EXPORT extern int rt_view_context_mouse_state_set(void *ctx, int x, int y);
RT_EXPORT extern int rt_view_context_interactive_rect_state_get(struct rt_view_interactive_rect_state *record, const void *ctx);
RT_EXPORT extern int rt_view_context_interactive_rect_state_set(void *ctx, const struct rt_view_interactive_rect_state *record);
RT_EXPORT extern void *rt_view_context_pick_point(void *ctx, int x, int y, int first_only);
RT_EXPORT extern void *rt_view_context_pick_nearest(void *ctx, int x, int y);
RT_EXPORT extern void *rt_view_context_pick_rect(void *ctx, int x0, int y0, int x1, int y1);
RT_EXPORT extern void *rt_view_context_pick_semantic_path(void *ctx, const char *path_pattern);
RT_EXPORT extern void *rt_view_pick_result_context_create(void);
RT_EXPORT extern void rt_view_pick_result_context_free(void *result_ctx);
RT_EXPORT extern size_t rt_view_pick_result_context_count(const void *result_ctx);
RT_EXPORT extern int rt_view_pick_result_context_path(const void *result_ctx, size_t index, struct bu_vls *path_out);
RT_EXPORT extern fastf_t rt_view_pick_result_context_hit_dist(const void *result_ctx, size_t index);
RT_EXPORT extern int rt_view_pick_result_context_append_path(void *result_ctx, void *view_ctx, int screen_x, int screen_y, const char *source_path, fastf_t hit_dist);
RT_EXPORT extern int rt_view_pick_result_context_append_copy(void *dest_ctx, const void *src_ctx, size_t index, fastf_t hit_dist);
RT_EXPORT extern void *rt_view_pick_result_context_filter_first(const void *src_ctx);
RT_EXPORT extern int rt_view_context_scene_adapter_set(void *ctx, const struct rt_view_context_scene_adapter *adapter);
RT_EXPORT extern int rt_view_context_scene_adapter_get(void *ctx, struct rt_view_context_scene_adapter *adapter);
RT_EXPORT extern int rt_view_context_selection_adapter_set(void *ctx, const struct rt_view_context_selection_adapter *adapter);
RT_EXPORT extern int rt_view_context_selection_adapter_get(void *ctx, struct rt_view_context_selection_adapter *adapter);
RT_EXPORT extern int rt_view_context_selection_available(void *ctx);
RT_EXPORT extern size_t rt_view_context_selection_count(void *ctx);
RT_EXPORT extern int rt_view_context_selection_set_pick_result_context(void *ctx, const void *result_ctx, rt_view_selection_path_callback_t callback, void *data);
RT_EXPORT extern int rt_view_context_selection_clear(void *ctx);
RT_EXPORT extern int rt_view_context_adc_state_get(struct rt_view_adc_state *record, const void *ctx);
RT_EXPORT extern int rt_view_context_adc_state_set(void *ctx, const struct rt_view_adc_state *record);
RT_EXPORT extern int rt_view_context_grid_state_get(struct rt_view_grid_state *record, const void *ctx);
RT_EXPORT extern int rt_view_context_grid_state_set(void *ctx, const struct rt_view_grid_state *record);
RT_EXPORT extern int rt_view_context_model_axes_state_get(struct rt_view_axes_state *record, const void *ctx);
RT_EXPORT extern int rt_view_context_model_axes_state_set(void *ctx, const struct rt_view_axes_state *record);
RT_EXPORT extern int rt_view_context_view_axes_state_get(struct rt_view_axes_state *record, const void *ctx);
RT_EXPORT extern int rt_view_context_view_axes_state_set(void *ctx, const struct rt_view_axes_state *record);
RT_EXPORT extern int rt_view_context_center_dot_state_get(struct rt_view_other_state *record, const void *ctx);
RT_EXPORT extern int rt_view_context_center_dot_state_set(void *ctx, const struct rt_view_other_state *record);
RT_EXPORT extern int rt_view_context_scale_overlay_state_get(struct rt_view_other_state *record, const void *ctx);
RT_EXPORT extern int rt_view_context_scale_overlay_state_set(void *ctx, const struct rt_view_other_state *record);
RT_EXPORT extern int rt_view_context_params_state_get(struct rt_view_params_state *record, const void *ctx);
RT_EXPORT extern int rt_view_context_params_state_set(void *ctx, const struct rt_view_params_state *record);
RT_EXPORT extern int rt_view_context_visible_render_summary(void *ctx, struct rt_view_render_summary *summary);
RT_EXPORT extern int rt_view_context_named_line_render_count(void *ctx, const char *name);
RT_EXPORT extern int rt_view_context_render_export_consistency(void *ctx, const char *drawn_prefix, struct rt_view_render_export_consistency *summary);
RT_EXPORT extern int rt_view_context_feature_geometry_summary(void *ctx, const char *name, struct rt_view_feature_geometry_summary *summary);
RT_EXPORT extern void *rt_view_set_context_create(void);
RT_EXPORT extern void rt_view_set_context_destroy(void *view_set_ctx);
RT_EXPORT extern void rt_view_set_context_init(void *view_set_ctx);
RT_EXPORT extern void rt_view_set_context_free(void *view_set_ctx);
RT_EXPORT extern struct bu_ptbl *rt_view_set_context_views(void *view_set_ctx);
RT_EXPORT extern void *rt_view_set_context_recycle_pool(void *view_set_ctx);
RT_EXPORT extern void *rt_view_set_context_find_view(void *view_set_ctx, const char *name);
RT_EXPORT extern int rt_view_set_context_add(void *view_set_ctx, void *view_ctx);
RT_EXPORT extern int rt_view_set_context_remove(void *view_set_ctx, void *view_ctx);
RT_EXPORT extern int rt_view_context_snap_grid_2d(void *ctx, fastf_t *vx, fastf_t *vy);
RT_EXPORT extern unsigned long long rt_view_context_prepare_tcl_snap(void *ctx);
RT_EXPORT extern int rt_view_context_measure_candidates(void *ctx, point_t a, point_t b, struct rt_view_measure_result *out);
RT_EXPORT extern void *rt_view_snap_result_context_create(void);
RT_EXPORT extern void rt_view_snap_result_context_free(void *result_ctx);
RT_EXPORT extern size_t rt_view_snap_result_context_count(const void *result_ctx);
RT_EXPORT extern int rt_view_snap_result_context_point(const void *result_ctx, size_t index, point_t point_out);
RT_EXPORT extern fastf_t rt_view_snap_result_context_distance(const void *result_ctx, size_t index);
RT_EXPORT extern unsigned long long rt_view_snap_result_context_kind(const void *result_ctx, size_t index);
RT_EXPORT extern int rt_view_snap_result_context_source_path(const void *result_ctx, size_t index, struct bu_vls *path_out);
RT_EXPORT extern int rt_view_context_snap_candidates_result(void *ctx, point_t sample, double tol, unsigned long long kinds, void *out_ctx);
RT_EXPORT extern int rt_view_context_snap_point_2d(void *ctx, fastf_t *vx, fastf_t *vy, unsigned long long kinds);
RT_EXPORT extern int rt_view_context_snap_first_candidate(void *ctx, const point_t sample, unsigned long long kinds, point_t candidate);
RT_EXPORT extern int rt_view_context_snap_lines_get(const void *ctx);
RT_EXPORT extern int rt_view_context_snap_lines_set(void *ctx, int enabled);
RT_EXPORT extern int rt_view_context_snap_source_flags_get(const void *ctx);
RT_EXPORT extern int rt_view_context_snap_source_flags_set(void *ctx, int flags);
RT_EXPORT extern unsigned long long rt_view_context_snap_kind_mask_get(const void *ctx);
RT_EXPORT extern int rt_view_context_snap_exclude_feature_clear(void *ctx);
RT_EXPORT extern int rt_view_context_center_linesnap(void *ctx);
RT_EXPORT extern double rt_view_context_snap_tolerance_factor_get(const void *ctx);
RT_EXPORT extern int rt_view_context_snap_tolerance_factor_set(void *ctx, double factor);
RT_EXPORT extern int rt_view_context_lod_policy_get(struct rt_view_lod_policy *policy, const void *ctx);
RT_EXPORT extern int rt_view_context_lod_policy_apply(void *ctx, const struct rt_view_lod_policy *policy);
RT_EXPORT extern int rt_view_context_lod_policy_copy(void *dst_ctx, const void *src_ctx);
RT_EXPORT extern int rt_view_context_lod_bounds_update(void *ctx);
RT_EXPORT extern int rt_view_context_lod_bounds_callback_set(void *ctx);
RT_EXPORT extern int rt_view_context_lod_bounds_callback_is(const void *ctx);
RT_EXPORT extern void *rt_view_context_bounds_update_suspend(void *ctx);
RT_EXPORT extern int rt_view_context_bounds_update_restore(void *ctx, void *state_ctx, int refresh_bounds);
RT_EXPORT extern rt_view_scene_ref rt_view_scene_ref_null(void);
RT_EXPORT extern rt_view_scene_ref rt_view_scene_ref_make(void *opaque, unsigned int backend);
RT_EXPORT extern int rt_view_scene_ref_is_null(rt_view_scene_ref ref);
RT_EXPORT extern int rt_view_scene_ref_equal(rt_view_scene_ref a, rt_view_scene_ref b);
RT_EXPORT extern void *rt_view_scene_ref_context(rt_view_scene_ref ref);
RT_EXPORT extern unsigned int rt_view_scene_ref_backend(rt_view_scene_ref ref);
RT_EXPORT extern rt_view_scene_ref rt_view_context_independent_scope_ref(void *ctx, int create);
RT_EXPORT extern int rt_view_context_is_independent(const void *ctx);
RT_EXPORT extern int rt_view_context_independent_scope_is_null(void *ctx, int create);
RT_EXPORT extern void rt_view_context_independent_scope_destroy(void *ctx);
RT_EXPORT extern rt_view_scene_ref rt_view_context_scene_root_ref(const void *ctx);
RT_EXPORT extern int rt_view_context_scene_root_ref_attach(void *ctx, rt_view_scene_ref root_ref);
RT_EXPORT extern int rt_view_context_scene_attached(const void *ctx);
RT_EXPORT extern int rt_view_context_scene_anchor_ensure(void *ctx);
RT_EXPORT extern int rt_view_context_scene_shared(const void *a_ctx, const void *b_ctx);
RT_EXPORT extern int rt_view_context_framebuffer_mode_get(const void *ctx);
RT_EXPORT extern int rt_view_context_framebuffer_mode_set(void *ctx, int mode);
RT_EXPORT extern int rt_view_context_settings_shared(const void *a, const void *b);
RT_EXPORT extern void rt_view_adc_model_to_view(struct rt_view_adc_state *adcs, mat_t model2view, fastf_t amax);
RT_EXPORT extern void rt_view_adc_grid_to_view(struct rt_view_adc_state *adcs, mat_t model2view, fastf_t amax);
RT_EXPORT extern void rt_view_adc_view_to_grid(struct rt_view_adc_state *adcs, mat_t model2view);
RT_EXPORT extern void rt_view_adc_reset(struct rt_view_adc_state *adcs, mat_t view2model, mat_t model2view);

/**
 * NOTE: Normally, librt doesn't have a concept of a "display" of the geometry.
 * However for at least the plotting routines, view information is sometimes
 * needed to produce more intelligent output.  In those situations, the
 * application may pass in an rt_view_info snapshot.
 */

/**
 * Specifies a subset of a primitive's geometry as the target for an
 * operation.
 *
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_selection {
    void *obj; /**< @brief primitive-specific selection object */
};

/**
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_selection_set {
    struct bu_ptbl selections; /**< @brief holds struct rt_selection */

    /** selection-object-specific routine that will free all memory
     *  associated with any of the stored selections
     */
    void (*free_selection)(struct rt_selection *);
};

/**
 * Stores selections associated with an object. There is an entry in
 * the selections table for each kind of selection (e.g. "active",
 * "option"). The table entries are sets to allow more than one
 * selection of the same type (e.g. multiple "option" selections).
 *
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_object_selections {
    /** selection type -> struct rt_selection_set */
    struct bu_hash_tbl *sets;
};

/**
 * Analogous to a database query. Specifies how to filter and sort the
 * selectable components of a primitive in order to find the most
 * relevant selections for a particular application.
 *
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_selection_query {
    point_t start;     /**< @brief start point of query ray */
    vect_t dir;        /**< @brief direction of query ray */

#define RT_SORT_UNSORTED         0
#define RT_SORT_CLOSEST_TO_START 1
    int sorting;
};

/**
 * Parameters of a translation applied to a selection.
 *
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_selection_translation {
    fastf_t dx;
    fastf_t dy;
    fastf_t dz;
};

/**
 * Describes an operation that can be applied to a selection.
 *
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_selection_operation {
#define RT_SELECTION_NOP         0
#define RT_SELECTION_TRANSLATION 1
    int type;
    union {
	struct rt_selection_translation tran;
    } parameters;
};

__END_DECLS

#endif /* RT_VIEW_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
