/*                    B V / V I E W . H
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
/** @file bv/view.h
 *
 * Scene-graph-free view state and manipulation logic.  This is a new libbv
 * boundary, not compatibility for the historical bview data model or full
 * old libbv API surface.
 */

#ifndef BV_VIEW_H
#define BV_VIEW_H

#include "common.h"

#include <stdint.h>

#include "vmath.h"
#include "bu/hash.h"
#include "bu/ptbl.h"
#include "bu/vls.h"

#ifndef BV_EXPORT
#  if defined(BV_DLL_EXPORTS) && defined(BV_DLL_IMPORTS)
#    error "Only BV_DLL_EXPORTS or BV_DLL_IMPORTS can be defined, not both."
#  elif defined(BV_DLL_EXPORTS)
#    define BV_EXPORT COMPILER_DLLEXPORT
#  elif defined(BV_DLL_IMPORTS)
#    define BV_EXPORT COMPILER_DLLIMPORT
#  else
#    define BV_EXPORT
#  endif
#endif

__BEGIN_DECLS

#define BV_STATE_MAGIC 0x42565657u
#define BV_SET_MAGIC 0x42565653u
#define BV_CONTEXT_MAGIC 0x42565643u
#define BV_CONTEXT_SET_MAGIC 0x42565654u

#define BV_DEFAULT_SCALE 500.0
#define BV_MIN_SIZE 0.0001
#define BV_MIN_SCALE 0.00005
#define BV_AUTOVIEW_SCALE_DEFAULT -1.0

/* Historical normalized view-coordinate range used by ADC/grid style tools. */
#define BV_VIEW_MAX 2047.0
#define BV_VIEW_MIN -2048.0
#define BV_VIEW_RANGE 4095.0
#define BV_INV_VIEW 0.00048828125
#define BV_INV_4096 0.000244140625

#define BV_ADJUST_IDLE      0x000ULL
#define BV_ADJUST_ROT       0x001ULL
#define BV_ADJUST_TRANS     0x002ULL
#define BV_ADJUST_SCALE     0x004ULL
#define BV_ADJUST_CENTER    0x008ULL
#define BV_ADJUST_CON_X     0x010ULL
#define BV_ADJUST_CON_Y     0x020ULL
#define BV_ADJUST_CON_Z     0x040ULL
#define BV_ADJUST_CON_GRID  0x080ULL
#define BV_ADJUST_CON_LINES 0x100ULL

#define BV_SNAP_SHARED 0x1
#define BV_SNAP_LOCAL  0x2
#define BV_SNAP_DB     0x4
#define BV_SNAP_VIEW   0x8
#define BV_SNAP_TCL    0x10

#define BV_SNAP_KIND_GRID           0x01ULL
#define BV_SNAP_KIND_ENDPOINT       0x02ULL
#define BV_SNAP_KIND_MIDPOINT       0x04ULL
#define BV_SNAP_KIND_INTERSECTION   0x08ULL
#define BV_SNAP_KIND_PERPENDICULAR  0x10ULL
#define BV_SNAP_KIND_TANGENT        0x20ULL
#define BV_SNAP_KIND_OVERLAY_HANDLE 0x40ULL
#define BV_SNAP_KIND_DEFAULT_MASK (BV_SNAP_KIND_GRID | BV_SNAP_KIND_ENDPOINT | \
	BV_SNAP_KIND_MIDPOINT | BV_SNAP_KIND_INTERSECTION | \
	BV_SNAP_KIND_PERPENDICULAR | BV_SNAP_KIND_TANGENT | \
	BV_SNAP_KIND_OVERLAY_HANDLE)

#define BV_CONTEXT_CHANGED_VIEW    0x00000001ULL
#define BV_CONTEXT_CHANGED_REFRESH 0x00000002ULL
#define BV_CONTEXT_CHANGED_ALL     0xffffffffffffffffULL

#define BV_REFRESH_VIEW        0x00000001u
#define BV_REFRESH_DRAW        0x00000002u
#define BV_REFRESH_EDIT        0x00000004u
#define BV_REFRESH_FRAMEBUFFER 0x00000008u
#define BV_REFRESH_OVERLAY     0x00000010u
#define BV_REFRESH_FORCE       0x80000000u
#define BV_REFRESH_ALL         0xffffffffu

enum bv_framebuffer_mode {
    BV_FRAMEBUFFER_MODE_OFF = 0,
    BV_FRAMEBUFFER_MODE_OVERLAY = 1,
    BV_FRAMEBUFFER_MODE_UNDERLAY = 2,
    /* Above model geometry but below view-local screen features. */
    BV_FRAMEBUFFER_MODE_INTERLAY = 3
};

struct bv_lod_settings {
    fastf_t scale;
    fastf_t curve_scale;
    fastf_t point_scale;
    size_t bot_threshold;
};

enum bv_lod_policy_mode {
    BV_LOD_OFF = 0,
    BV_LOD_AUTO = 1,
    BV_LOD_FORCE_LEVEL = 2
};

struct bv_lod_policy {
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

struct bv_view_info {
    int width;
    int height;
    fastf_t size;
    struct bv_lod_settings lod;
};

#define BV_LOD_SETTINGS_INIT { 1.0, 1.0, 1.0, 0 }
#define BV_LOD_POLICY_INIT { BV_LOD_AUTO, 0, 0, 0, 0, 0, 1.0, 1.0, 1.0 }
#define BV_VIEW_INFO_INIT { 1, 1, 1.0, BV_LOD_SETTINGS_INIT }

#define BV_INTERACTIVE_RECT_STATE_INIT { 0, 0, 0, 0, {0, 0}, {0, 0}, 0.0, 0.0, 0.0, 0.0, {0, 0, 0}, {0, 0, 0}, {0, 0}, 0.0 }
#define BV_MEASURE_RESULT_INIT { 0.0, 0.0, 0.0, 0 }
#define BV_ADC_STATE_INIT { 0, 0, 0, 0, 0, 0, VINIT_ZERO, VINIT_ZERO, VINIT_ZERO, 0.0, 0.0, 0.0, 0, 0, 0, 0, VINIT_ZERO, VINIT_ZERO, VINIT_ZERO, {0, 0, 0}, {0, 0, 0}, 0 }
#define BV_GRID_STATE_INIT { 0, 0, 0, 0, VINIT_ZERO, 0.0, 0.0, 0, 0, {0, 0, 0} }
#define BV_AXES_STATE_INIT { 0, VINIT_ZERO, 0.0, 0, {0, 0, 0}, 0, 0, {0, 0, 0}, 0, 0, 0, 0, 0.0, 0, 0, {0, 0, 0}, {0, 0, 0} }
#define BV_OTHER_STATE_INIT { 0, {0, 0, 0}, {0, 0, 0}, 0 }
#define BV_PARAMS_STATE_INIT { 0, 0, 0, 0, 0, 0, 0, {0, 0, 0}, 0 }
#define BV_SNAP_STATE_INIT { 0, 0, BV_SNAP_KIND_DEFAULT_MASK, 10.0 }
#define BV_MOUSE_DELTA_SETTINGS_INIT { 0.0, 0.0, 0.0, 0.0 }

enum bv_knobs_category {
    BV_KNOBS_ALL = 0,
    BV_KNOBS_RATE = 1,
    BV_KNOBS_ABS = 2
};

struct bv_knobs {
    vect_t rot_model;
    int rot_model_active;
    char rot_model_origin;
    void *rot_model_data;

    vect_t rot_object;
    int rot_object_active;
    char rot_object_origin;
    void *rot_object_data;

    vect_t rot_view;
    int rot_view_active;
    char rot_view_origin;
    void *rot_view_data;

    fastf_t scale_rate;
    int scale_rate_active;
    void *scale_rate_data;

    vect_t trans_model;
    int trans_model_active;
    void *trans_model_data;

    vect_t trans_view;
    int trans_view_active;
    void *trans_view_data;

    vect_t abs_rot_model;
    vect_t abs_rot_model_last;
    vect_t abs_rot_object;
    vect_t abs_rot_object_last;
    vect_t abs_rot_view;
    vect_t abs_rot_view_last;

    fastf_t abs_scale;

    vect_t abs_trans_model;
    vect_t abs_trans_model_last;
    vect_t abs_trans_view;
    vect_t abs_trans_view_last;
};

struct bv_knob_values {
    vect_t rate_rotation;
    vect_t rate_translation;
    fastf_t rate_scale;
    vect_t absolute_rotation;
    vect_t absolute_translation;
    fastf_t absolute_scale;
};

struct bv_interactive_rect_state {
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

struct bv_measure_result {
    fastf_t distance;
    fastf_t projection;
    fastf_t normal_alignment;
    int valid;
};

struct bv_adc_state {
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

struct bv_grid_state {
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

struct bv_axes_state {
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

struct bv_other_state {
    int gos_draw;
    int gos_line_color[3];
    int gos_text_color[3];
    int gos_font_size;
};

struct bv_params_state {
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

#define BV_BACKGROUND_STATE_INIT {{0, 0, 0}, {0, 0, 0}}
struct bv_background_state {
    int bottom[3];
    int top[3];
};

struct bv_snap_state {
    int lines;
    int source_flags;
    unsigned long long kind_mask;
    double tolerance_factor;
};

struct bv_mouse_delta_settings {
    fastf_t min_delta;
    fastf_t max_delta;
    fastf_t rotate_scale;
    fastf_t scale_scale;
};

struct bv_set;
struct bv_context;
struct bv_context_set;

typedef void (*bv_context_callback_t)(struct bv_context *ctx,
	uint64_t changed_flags,
	void *client_data);

struct bv {
    uint32_t magic;
    struct bu_vls name;
    void *user_data;
    struct bv_set *view_set;

    int width;
    int height;
    fastf_t scale;
    fastf_t initial_scale;
    fastf_t absolute_scale;
    fastf_t size;
    fastf_t inverse_size;
    fastf_t perspective;
    fastf_t local2base;
    fastf_t base2local;
    uint64_t frame_revision;
    uint32_t refresh_dirty;
    int refresh_enabled;
    int refresh_suppressed;
    int refresh_drawn_count;
    uint64_t frametime;
    int zclip;
    int framebuffer_mode;
    int cleared;

    vect_t aet;
    point_t eye_pos;
    point_t keypoint;
    point_t current_point;
    int mouse_x;
    int mouse_y;
    fastf_t previous_mouse_x;
    fastf_t previous_mouse_y;
    char coord_mode;
    char rotate_about;

    fastf_t min_mouse_delta;
    fastf_t max_mouse_delta;
    fastf_t rotate_scale;
    fastf_t scale_scale;

    mat_t rotation;
    mat_t center;
    mat_t model2view;
    mat_t pmodel2view;
    mat_t view2model;
    mat_t pmat;

    point_t obb_center;
    vect_t obb_extent1;
    vect_t obb_extent2;
    vect_t obb_extent3;
    fastf_t radius;

    struct bv_knobs knobs;
    struct bv_interactive_rect_state interactive_rect;
    struct bv_adc_state adc;
    struct bv_grid_state grid;
    struct bv_axes_state model_axes;
    struct bv_axes_state view_axes;
    struct bv_other_state center_dot;
    struct bv_other_state scale_overlay;
    struct bv_params_state params;
    struct bv_background_state background;
    struct bv_snap_state snap;
};

struct bv_set {
    uint32_t magic;
    struct bu_ptbl views;
};

struct bv_context {
    uint32_t magic;
    struct bv view;
    struct bv *view_ptr;
    int owns_view;
    struct bv_context_set *view_set;
    struct bu_ptbl callbacks;
};

struct bv_context_set {
    uint32_t magic;
    struct bv_set bvset;
    struct bu_ptbl views;
    struct bu_ptbl recycle_pool;
};

BV_EXPORT extern struct bv *bv_create(void);
BV_EXPORT extern void bv_destroy(struct bv *v);
BV_EXPORT extern void bv_init(struct bv *v);
BV_EXPORT extern void bv_free(struct bv *v);
BV_EXPORT extern int bv_is_valid(const struct bv *v);
BV_EXPORT extern int bv_copy(struct bv *dst, const struct bv *src);

BV_EXPORT extern int bv_name_set(struct bv *v, const char *name);
BV_EXPORT extern const char *bv_name_get(const struct bv *v);
BV_EXPORT extern int bv_user_data_set(struct bv *v, void *user_data);
BV_EXPORT extern void *bv_user_data_get(const struct bv *v);
BV_EXPORT extern int bv_dimensions_set(struct bv *v, int width, int height);
BV_EXPORT extern int bv_width_get(const struct bv *v);
BV_EXPORT extern int bv_height_get(const struct bv *v);
BV_EXPORT extern fastf_t bv_scale_get(const struct bv *v);
BV_EXPORT extern fastf_t *bv_scale_storage_get(struct bv *v);
BV_EXPORT extern fastf_t bv_initial_scale_get(const struct bv *v);
BV_EXPORT extern int bv_initial_scale_set(struct bv *v, fastf_t scale);
BV_EXPORT extern fastf_t bv_absolute_scale_get(const struct bv *v);
BV_EXPORT extern int bv_absolute_scale_set(struct bv *v, fastf_t scale);
BV_EXPORT extern fastf_t bv_size_get(const struct bv *v);
BV_EXPORT extern fastf_t bv_inverse_size_get(const struct bv *v);
BV_EXPORT extern fastf_t bv_radius_get(const struct bv *v);
BV_EXPORT extern fastf_t bv_perspective_get(const struct bv *v);
BV_EXPORT extern int bv_perspective_set(struct bv *v, fastf_t perspective);
BV_EXPORT extern fastf_t bv_local2base_get(const struct bv *v);
BV_EXPORT extern fastf_t bv_base2local_get(const struct bv *v);
BV_EXPORT extern int bv_unit_conversion_set(struct bv *v, fastf_t local2base, fastf_t base2local);
BV_EXPORT extern char bv_coord_get(const struct bv *v);
BV_EXPORT extern int bv_coord_set(struct bv *v, char coord);
BV_EXPORT extern char bv_rotate_about_get(const struct bv *v);
BV_EXPORT extern int bv_rotate_about_set(struct bv *v, char rotate_about);
BV_EXPORT extern int bv_refresh_request(struct bv *v, uint32_t flags);
BV_EXPORT extern int bv_refresh_dirty_get(const struct bv *v);
BV_EXPORT extern uint32_t bv_refresh_consume(struct bv *v);
BV_EXPORT extern int bv_refresh_complete(struct bv *v);
BV_EXPORT extern int bv_refresh_enabled_get(const struct bv *v);
BV_EXPORT extern int bv_refresh_enabled_set(struct bv *v, int enabled);
BV_EXPORT extern int bv_refresh_suppressed_get(const struct bv *v);
BV_EXPORT extern int bv_refresh_suppress_begin(struct bv *v);
BV_EXPORT extern int bv_refresh_suppress_end(struct bv *v);
BV_EXPORT extern int bv_refresh_drawn_count_get(const struct bv *v);
BV_EXPORT extern int bv_refresh_drawn_count_set(struct bv *v, int count);
BV_EXPORT extern uint64_t bv_frame_revision_get(const struct bv *v);
BV_EXPORT extern uint64_t bv_frame_revision_bump(struct bv *v);
BV_EXPORT extern unsigned long long bv_hash(const struct bv *v);
BV_EXPORT extern int bv_frametime_set(struct bv *v, uint64_t frametime);
BV_EXPORT extern uint64_t bv_frametime_get(const struct bv *v);
BV_EXPORT extern int bv_zclip_get(const struct bv *v);
BV_EXPORT extern int bv_zclip_set(struct bv *v, int zclip);
BV_EXPORT extern int bv_framebuffer_mode_get(const struct bv *v);
BV_EXPORT extern int bv_framebuffer_mode_set(struct bv *v, int mode);
BV_EXPORT extern int bv_cleared_get(const struct bv *v);
BV_EXPORT extern int bv_cleared_set(struct bv *v, int cleared);

BV_EXPORT extern void bv_faceplate_defaults(struct bv *v);
BV_EXPORT extern void bv_snap_defaults(struct bv *v);
BV_EXPORT extern int bv_interactive_rect_state_get(struct bv_interactive_rect_state *record, const struct bv *v);
BV_EXPORT extern int bv_interactive_rect_state_set(struct bv *v, const struct bv_interactive_rect_state *record);
BV_EXPORT extern int bv_adc_state_get(struct bv_adc_state *record, const struct bv *v);
BV_EXPORT extern int bv_adc_state_set(struct bv *v, const struct bv_adc_state *record);
BV_EXPORT extern void bv_adc_model_to_view(struct bv_adc_state *adcs, mat_t model2view, fastf_t amax);
BV_EXPORT extern void bv_adc_grid_to_view(struct bv_adc_state *adcs, mat_t model2view, fastf_t amax);
BV_EXPORT extern void bv_adc_view_to_grid(struct bv_adc_state *adcs, mat_t model2view);
BV_EXPORT extern void bv_adc_reset(struct bv_adc_state *adcs, mat_t view2model, mat_t model2view);
BV_EXPORT extern int bv_grid_state_get(struct bv_grid_state *record, const struct bv *v);
BV_EXPORT extern int bv_grid_state_set(struct bv *v, const struct bv_grid_state *record);
BV_EXPORT extern int bv_model_axes_state_get(struct bv_axes_state *record, const struct bv *v);
BV_EXPORT extern int bv_model_axes_state_set(struct bv *v, const struct bv_axes_state *record);
BV_EXPORT extern int bv_view_axes_state_get(struct bv_axes_state *record, const struct bv *v);
BV_EXPORT extern int bv_view_axes_state_set(struct bv *v, const struct bv_axes_state *record);
BV_EXPORT extern int bv_center_dot_state_get(struct bv_other_state *record, const struct bv *v);
BV_EXPORT extern int bv_center_dot_state_set(struct bv *v, const struct bv_other_state *record);
BV_EXPORT extern int bv_scale_overlay_state_get(struct bv_other_state *record, const struct bv *v);
BV_EXPORT extern int bv_scale_overlay_state_set(struct bv *v, const struct bv_other_state *record);
BV_EXPORT extern int bv_params_state_get(struct bv_params_state *record, const struct bv *v);
BV_EXPORT extern int bv_params_state_set(struct bv *v, const struct bv_params_state *record);
BV_EXPORT extern int bv_background_state_get(struct bv_background_state *record, const struct bv *v);
BV_EXPORT extern int bv_background_state_set(struct bv *v, const struct bv_background_state *record);
BV_EXPORT extern int bv_snap_state_get(struct bv_snap_state *record, const struct bv *v);
BV_EXPORT extern int bv_snap_state_set(struct bv *v, const struct bv_snap_state *record);
BV_EXPORT extern int bv_snap_lines_get(const struct bv *v);
BV_EXPORT extern int bv_snap_lines_set(struct bv *v, int enabled);
BV_EXPORT extern int bv_snap_source_flags_get(const struct bv *v);
BV_EXPORT extern int bv_snap_source_flags_set(struct bv *v, int flags);
BV_EXPORT extern unsigned long long bv_snap_kind_mask_get(const struct bv *v);
BV_EXPORT extern int bv_snap_kind_mask_set(struct bv *v, unsigned long long mask);
BV_EXPORT extern double bv_snap_tolerance_factor_get(const struct bv *v);
BV_EXPORT extern int bv_snap_tolerance_factor_set(struct bv *v, double factor);

BV_EXPORT extern void bv_view_info_init(struct bv_view_info *info);
BV_EXPORT extern void bv_view_info_sanitize(struct bv_view_info *info);
BV_EXPORT extern void bv_lod_policy_init(struct bv_lod_policy *policy);
BV_EXPORT extern void bv_lod_policy_sanitize(struct bv_lod_policy *policy);
BV_EXPORT extern fastf_t bv_view_lod_curve_scale(const struct bv_view_info *info);
BV_EXPORT extern size_t bv_view_lod_bot_threshold(const struct bv_view_info *info);
BV_EXPORT extern fastf_t bv_view_avg_sample_spacing(const struct bv_view_info *info);
BV_EXPORT extern fastf_t bv_view_solid_point_spacing(const struct bv_view_info *info, fastf_t solid_width);

BV_EXPORT extern int bv_update(struct bv *v);
BV_EXPORT extern int bv_model2view_get(mat_t model2view, const struct bv *v);
BV_EXPORT extern int bv_model2view_set(struct bv *v, const mat_t model2view);
BV_EXPORT extern int bv_view2model_get(mat_t view2model, const struct bv *v);
BV_EXPORT extern int bv_view2model_set(struct bv *v, const mat_t view2model);
BV_EXPORT extern int bv_pmodel2view_get(mat_t pmodel2view, const struct bv *v);
BV_EXPORT extern int bv_pmodel2view_set(struct bv *v, const mat_t pmodel2view);
BV_EXPORT extern int bv_pmat_get(mat_t pmat, const struct bv *v);
BV_EXPORT extern int bv_pmat_set(struct bv *v, const mat_t pmat);
BV_EXPORT extern int bv_rotation_get(mat_t rotation, const struct bv *v);
BV_EXPORT extern int bv_rotation_set(struct bv *v, const mat_t rotation);
BV_EXPORT extern int bv_orientation_quat_get(quat_t orientation, const struct bv *v);
BV_EXPORT extern int bv_center_mat_get(mat_t center, const struct bv *v);
BV_EXPORT extern int bv_center_mat_set(struct bv *v, const mat_t center);
BV_EXPORT extern int bv_center_get(point_t center, const struct bv *v);
BV_EXPORT extern int bv_center_set(struct bv *v, const point_t center);
BV_EXPORT extern int bv_scale_set(struct bv *v, fastf_t scale);
BV_EXPORT extern int bv_scale_state_set(struct bv *v, fastf_t scale, fastf_t initial_scale, fastf_t absolute_scale, fastf_t size, fastf_t inverse_size);
BV_EXPORT extern int bv_size_set(struct bv *v, fastf_t size);
BV_EXPORT extern int bv_aet_get(vect_t aet, const struct bv *v);
BV_EXPORT extern int bv_aet_set(struct bv *v, const vect_t aet);
BV_EXPORT extern int bv_aet_state_set(struct bv *v, const vect_t aet);
BV_EXPORT extern int bv_eye_pos_get(point_t eye_pos, const struct bv *v);
BV_EXPORT extern int bv_eye_pos_set(struct bv *v, const point_t eye_pos);
BV_EXPORT extern int bv_keypoint_get(point_t keypoint, const struct bv *v);
BV_EXPORT extern int bv_keypoint_set(struct bv *v, const point_t keypoint);
BV_EXPORT extern int bv_current_point_get(point_t current_point, const struct bv *v);
BV_EXPORT extern int bv_current_point_set(struct bv *v, const point_t current_point);

BV_EXPORT extern int bv_autoview_bounds(struct bv *v, fastf_t scale, const point_t min, const point_t max);
BV_EXPORT extern int bv_screen_to_view(fastf_t *vx, fastf_t *vy, const struct bv *v, fastf_t sx, fastf_t sy);
BV_EXPORT extern int bv_snap_grid_2d(const struct bv *v, fastf_t *vx, fastf_t *vy);
BV_EXPORT extern int bv_screen_to_model(point_t model_point, const struct bv *v, fastf_t sx, fastf_t sy);
BV_EXPORT extern int bv_plane_get(plane_t *plane, const struct bv *v);
BV_EXPORT extern int bv_previous_mouse_get(fastf_t *x, fastf_t *y, const struct bv *v);
BV_EXPORT extern int bv_previous_mouse_set(struct bv *v, fastf_t x, fastf_t y);
BV_EXPORT extern int bv_mouse_delta_settings_get(struct bv_mouse_delta_settings *settings, const struct bv *v);
/** Set validated, passive pointer-delta navigation policy for one view. */
BV_EXPORT extern int bv_mouse_delta_settings_set(struct bv *v, const struct bv_mouse_delta_settings *settings);
/** Clamp a signed pointer delta according to the view's navigation policy. */
BV_EXPORT extern int bv_mouse_delta_clamp(fastf_t *dx, fastf_t *dy, const struct bv *v);
/** Apply a normalized pointer drag using the view's navigation policy. */
BV_EXPORT extern int bv_mouse_delta_adjust(struct bv *v, int dx, int dy, const point_t keypoint, unsigned long long flags);
BV_EXPORT extern int bv_mouse_state_set(struct bv *v, int x, int y);
BV_EXPORT extern int bv_adjust(struct bv *v, int dx, int dy, const point_t keypoint, int mode, unsigned long long flags);

BV_EXPORT extern int bv_knobs_reset(struct bv_knobs *knobs, int category);
BV_EXPORT extern unsigned long long bv_knobs_hash(const struct bv_knobs *knobs, struct bu_data_hash_state *state);
BV_EXPORT extern int bv_knobs_cmd_process(vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, struct bv *v, const char *cmd, fastf_t factor, char origin, int model_flag, int incr_flag);
BV_EXPORT extern int bv_knobs_translate(struct bv *v, const vect_t tvec, int model_flag);
BV_EXPORT extern int bv_knobs_rotate(struct bv *v, const vect_t rvec, char origin, char coords, const matp_t obj_rot, const pointp_t pvt_pt);
BV_EXPORT extern int bv_knobs_update_rate_flags(struct bv *v);
BV_EXPORT extern int bv_knob_values_get(struct bv_knob_values *values, const struct bv *v);
BV_EXPORT extern int bv_knobs_calibrate(struct bv *v);

BV_EXPORT extern struct bv_context *bv_context_create(void);
BV_EXPORT extern void bv_context_destroy(struct bv_context *ctx);
BV_EXPORT extern void bv_context_init(struct bv_context *ctx, struct bv *view);
BV_EXPORT extern void bv_context_free(struct bv_context *ctx);
BV_EXPORT extern int bv_context_is_valid(const struct bv_context *ctx);
BV_EXPORT extern struct bv *bv_context_view(struct bv_context *ctx);
BV_EXPORT extern const struct bv *bv_context_view_const(const struct bv_context *ctx);
BV_EXPORT extern int bv_context_name_set(struct bv_context *ctx, const char *name);
BV_EXPORT extern const char *bv_context_name_get(const struct bv_context *ctx);
BV_EXPORT extern int bv_context_user_data_set(struct bv_context *ctx, void *user_data);
BV_EXPORT extern void *bv_context_user_data_get(const struct bv_context *ctx);
BV_EXPORT extern int bv_context_dimensions_set(struct bv_context *ctx, int width, int height);
BV_EXPORT extern int bv_context_width_get(const struct bv_context *ctx);
BV_EXPORT extern int bv_context_height_get(const struct bv_context *ctx);
BV_EXPORT extern int bv_context_refresh_request(struct bv_context *ctx, uint32_t flags);
BV_EXPORT extern int bv_context_refresh_complete(struct bv_context *ctx);
BV_EXPORT extern int bv_context_callback_add(struct bv_context *ctx, bv_context_callback_t callback, void *client_data);
BV_EXPORT extern int bv_context_callback_remove(struct bv_context *ctx, bv_context_callback_t callback, void *client_data);
BV_EXPORT extern void bv_context_callbacks_clear(struct bv_context *ctx);
BV_EXPORT extern int bv_context_notify(struct bv_context *ctx, uint64_t changed_flags);
BV_EXPORT extern int bv_context_update(struct bv_context *ctx, uint64_t changed_flags);
BV_EXPORT extern int bv_context_settings_shared(const struct bv_context *a, const struct bv_context *b);

BV_EXPORT extern struct bv_context_set *bv_context_set_create(void);
BV_EXPORT extern void bv_context_set_destroy(struct bv_context_set *s);
BV_EXPORT extern void bv_context_set_init(struct bv_context_set *s);
BV_EXPORT extern void bv_context_set_free(struct bv_context_set *s);
BV_EXPORT extern struct bu_ptbl *bv_context_set_views(struct bv_context_set *s);
BV_EXPORT extern void *bv_context_set_recycle_pool(struct bv_context_set *s);
BV_EXPORT extern struct bv_context *bv_context_set_find(struct bv_context_set *s, const char *name);
BV_EXPORT extern int bv_context_set_attach(struct bv_context_set *s, struct bv_context *ctx);
BV_EXPORT extern int bv_context_set_add(struct bv_context_set *s, struct bv_context *ctx);
BV_EXPORT extern int bv_context_set_remove(struct bv_context_set *s, struct bv_context *ctx);

BV_EXPORT extern struct bv_set *bv_set_create(void);
BV_EXPORT extern void bv_set_destroy(struct bv_set *s);
BV_EXPORT extern void bv_set_init(struct bv_set *s);
BV_EXPORT extern void bv_set_free(struct bv_set *s);
BV_EXPORT extern int bv_set_add(struct bv_set *s, struct bv *v);
BV_EXPORT extern int bv_set_remove(struct bv_set *s, struct bv *v);
BV_EXPORT extern struct bu_ptbl *bv_set_views(struct bv_set *s);
BV_EXPORT extern struct bv *bv_set_find(struct bv_set *s, const char *name);

__END_DECLS

#endif /* BV_VIEW_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
