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
#include "bu/color.h"
#include "bu/list.h"
#include "bu/hash.h"
#include "bu/ptbl.h"
#include "bn/tol.h"
#include "rt/defines.h"

__BEGIN_DECLS

struct db_i;
struct rt_mesh_lod;

/* Historical normalized view-coordinate range. */
#define RT_VIEW_MAX 2047.0
#define RT_VIEW_MIN -2048.0
#define RT_VIEW_RANGE 4095.0
#define RT_INV_VIEW 0.00048828125
#define RT_INV_4096 0.000244140625
#define RT_VIEW_MIN_SIZE 0.0001
#define RT_VIEW_MIN_SCALE 0.00005
#define RT_VIEW_AUTOVIEW_SCALE_DEFAULT -1

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
};

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

/* Borrowed active LoD arrays; valid until the LoD is reloaded or destroyed. */
struct rt_mesh_lod_data {
    const int *faces;
    size_t face_count;
    const point_t *points;
    size_t point_count;
    const point_t *points_orig;
    size_t point_orig_count;
    const vect_t *normals;
    size_t normal_count;
    point_t bmin;
    point_t bmax;
};

/* Full-detail mesh arrays supplied by producer callbacks.  The callback owns
 * the array lifetimes; RT borrows them until the matching clear/free callback.
 * Normal payloads are optional, but when present must contain one normal per
 * triangle corner: normal_count == face_count * 3. */
struct rt_mesh_lod_detail {
    const int *faces;
    size_t face_count;
    const point_t *points;
    size_t point_count;
    const point_t *points_orig;
    size_t point_orig_count;
    const vect_t *normals;
    size_t normal_count;
};

typedef int (*rt_mesh_lod_detail_setup_callback)(struct rt_mesh_lod_detail *detail, void *cb_data);
typedef int (*rt_mesh_lod_detail_clear_callback)(void *cb_data);
typedef int (*rt_mesh_lod_detail_free_callback)(void *cb_data);

/* Summary of the active mesh LoD state.  This is stable to copy and does not
 * borrow array storage from the LoD provider. */
struct rt_mesh_lod_info {
    int active_level;
    size_t face_count;
    size_t point_count;
    size_t point_orig_count;
    size_t normal_count;
    int has_faces;
    int has_points;
    int has_original_points;
    int has_snapped_points;
    int has_normals;
    point_t bmin;
    point_t bmax;
};

/* Stable status for a database object's mesh LoD cache entry.  The key is an
 * opaque provider cache key; callers should treat it as diagnostic metadata,
 * not as an addressable storage handle. */
struct rt_mesh_lod_cache_status {
    int directory_found;
    int is_bot;
    int has_cache_key;
    int has_cached_payload;
    int stale_cache_entry;
    int cleared_cache_entry;
    int generated_cache_entry;
    unsigned long long cache_key;
    unsigned long long cleared_cache_key;
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
#define RT_MESH_LOD_INFO_INIT { -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, VINIT_ZERO, VINIT_ZERO }
#define RT_MESH_LOD_CACHE_STATUS_INIT { 0, 0, 0, 0, 0, 0, 0, 0, 0 }

RT_EXPORT extern void rt_view_info_init(struct rt_view_info *info);
RT_EXPORT extern void rt_view_info_sanitize(struct rt_view_info *info);
RT_EXPORT extern void rt_view_lod_policy_init(struct rt_view_lod_policy *policy);
RT_EXPORT extern void rt_view_lod_policy_sanitize(struct rt_view_lod_policy *policy);
RT_EXPORT extern fastf_t rt_view_lod_curve_scale(const struct rt_view_info *info);
RT_EXPORT extern size_t rt_view_lod_bot_threshold(const struct rt_view_info *info);
RT_EXPORT extern fastf_t rt_view_avg_sample_spacing(const struct rt_view_info *info);
RT_EXPORT extern fastf_t rt_view_solid_point_spacing(const struct rt_view_info *info, fastf_t solid_width);
RT_EXPORT extern void rt_view_adc_model_to_view(struct rt_view_adc_state *adcs, mat_t model2view, fastf_t amax);
RT_EXPORT extern void rt_view_adc_grid_to_view(struct rt_view_adc_state *adcs, mat_t model2view, fastf_t amax);
RT_EXPORT extern void rt_view_adc_view_to_grid(struct rt_view_adc_state *adcs, mat_t model2view);
RT_EXPORT extern void rt_view_adc_reset(struct rt_view_adc_state *adcs, mat_t view2model, mat_t model2view);

/* Routines for managing the mesh LoD cache */
RT_EXPORT extern void db_mesh_lod_init(struct db_i *dbip, int verbose);
RT_EXPORT extern void db_mesh_lod_clear(struct db_i *dbip);
RT_EXPORT extern void rt_mesh_lod_cache_clear_all(void);
RT_EXPORT extern int db_mesh_lod_update(struct db_i *dbip, const char *name);
RT_EXPORT extern void rt_mesh_lod_cache_status_init(struct rt_mesh_lod_cache_status *status);
RT_EXPORT extern int db_mesh_lod_status(struct db_i *dbip, const char *name, struct rt_mesh_lod_cache_status *status);
RT_EXPORT extern int db_mesh_lod_refresh(struct db_i *dbip, const char *name, struct rt_mesh_lod_cache_status *status);
RT_EXPORT extern int db_mesh_lod_invalidate(struct db_i *dbip, const char *name, struct rt_mesh_lod_cache_status *status);
/* Store caller-owned mesh arrays in the database LoD cache.  normals is
 * optional; when supplied it must contain one normal per triangle corner in
 * faces order, i.e. face_count * 3 normals. */
RT_EXPORT extern int db_mesh_lod_store_mesh(struct db_i *dbip, const char *name, const point_t *vertices, size_t vertex_count, const vect_t *normals, const int *faces, size_t face_count, unsigned long long user_key, fastf_t fidelity_ratio, struct rt_mesh_lod_cache_status *status);
RT_EXPORT extern struct rt_mesh_lod *db_mesh_lod_get(struct db_i *dbip, const char *name);
RT_EXPORT extern int rt_mesh_lod_load_level(struct rt_mesh_lod *lod, int level, int reset);
RT_EXPORT extern int rt_mesh_lod_load_view(struct rt_mesh_lod *lod, const struct rt_view_info *info, int reset);
RT_EXPORT extern int rt_mesh_lod_current_level(const struct rt_mesh_lod *lod);
RT_EXPORT extern int rt_mesh_lod_has_active_data(const struct rt_mesh_lod *lod);
RT_EXPORT extern int rt_mesh_lod_data_get(const struct rt_mesh_lod *lod, struct rt_mesh_lod_data *data);
RT_EXPORT extern void rt_mesh_lod_info_init(struct rt_mesh_lod_info *info);
RT_EXPORT extern int rt_mesh_lod_info_get(const struct rt_mesh_lod *lod, struct rt_mesh_lod_info *info);
RT_EXPORT extern void rt_mesh_lod_detail_init(struct rt_mesh_lod_detail *detail);
RT_EXPORT extern int rt_mesh_lod_detail_callbacks_set(struct rt_mesh_lod *lod, rt_mesh_lod_detail_setup_callback setup_clbk, rt_mesh_lod_detail_clear_callback clear_clbk, rt_mesh_lod_detail_free_callback free_clbk, void *cb_data);
RT_EXPORT extern void rt_mesh_lod_detail_callbacks_clear(struct rt_mesh_lod *lod);
RT_EXPORT extern void rt_mesh_lod_memshrink(struct rt_mesh_lod *lod);
RT_EXPORT extern void rt_mesh_lod_destroy(struct rt_mesh_lod *lod);

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
