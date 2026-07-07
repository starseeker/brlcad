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

#define BV_DEFAULT_SCALE 500.0
#define BV_MIN_SIZE 0.0001
#define BV_MIN_SCALE 0.00005
#define BV_AUTOVIEW_SCALE_DEFAULT -1.0

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

struct bv_set;

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
};

struct bv_set {
    uint32_t magic;
    struct bu_ptbl views;
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
BV_EXPORT extern int bv_frametime_set(struct bv *v, uint64_t frametime);
BV_EXPORT extern uint64_t bv_frametime_get(const struct bv *v);
BV_EXPORT extern int bv_zclip_get(const struct bv *v);
BV_EXPORT extern int bv_zclip_set(struct bv *v, int zclip);
BV_EXPORT extern int bv_framebuffer_mode_get(const struct bv *v);
BV_EXPORT extern int bv_framebuffer_mode_set(struct bv *v, int mode);
BV_EXPORT extern int bv_cleared_get(const struct bv *v);
BV_EXPORT extern int bv_cleared_set(struct bv *v, int cleared);

BV_EXPORT extern int bv_update(struct bv *v);
BV_EXPORT extern int bv_model2view_get(mat_t model2view, const struct bv *v);
BV_EXPORT extern int bv_view2model_get(mat_t view2model, const struct bv *v);
BV_EXPORT extern int bv_pmodel2view_get(mat_t pmodel2view, const struct bv *v);
BV_EXPORT extern int bv_pmat_get(mat_t pmat, const struct bv *v);
BV_EXPORT extern int bv_pmat_set(struct bv *v, const mat_t pmat);
BV_EXPORT extern int bv_rotation_get(mat_t rotation, const struct bv *v);
BV_EXPORT extern int bv_rotation_set(struct bv *v, const mat_t rotation);
BV_EXPORT extern int bv_center_get(point_t center, const struct bv *v);
BV_EXPORT extern int bv_center_set(struct bv *v, const point_t center);
BV_EXPORT extern int bv_scale_set(struct bv *v, fastf_t scale);
BV_EXPORT extern int bv_size_set(struct bv *v, fastf_t size);
BV_EXPORT extern int bv_aet_get(vect_t aet, const struct bv *v);
BV_EXPORT extern int bv_aet_set(struct bv *v, const vect_t aet);
BV_EXPORT extern int bv_aet_state_set(struct bv *v, const vect_t aet);
BV_EXPORT extern int bv_eye_pos_get(point_t eye_pos, const struct bv *v);
BV_EXPORT extern int bv_eye_pos_set(struct bv *v, const point_t eye_pos);
BV_EXPORT extern int bv_keypoint_get(point_t keypoint, const struct bv *v);
BV_EXPORT extern int bv_keypoint_set(struct bv *v, const point_t keypoint);

BV_EXPORT extern int bv_autoview_bounds(struct bv *v, fastf_t scale, const point_t min, const point_t max);
BV_EXPORT extern int bv_screen_to_view(fastf_t *vx, fastf_t *vy, const struct bv *v, fastf_t sx, fastf_t sy);
BV_EXPORT extern int bv_screen_to_model(point_t model_point, const struct bv *v, fastf_t sx, fastf_t sy);
BV_EXPORT extern int bv_plane_get(plane_t *plane, const struct bv *v);
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
