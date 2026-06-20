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
 * the array lifetimes; RT borrows them until the matching clear/free callback. */
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

/* Routines for managing the mesh LoD cache */
RT_EXPORT extern void db_mesh_lod_init(struct db_i *dbip, int verbose);
RT_EXPORT extern void db_mesh_lod_clear(struct db_i *dbip);
RT_EXPORT extern void rt_mesh_lod_cache_clear_all(void);
RT_EXPORT extern int db_mesh_lod_update(struct db_i *dbip, const char *name);
RT_EXPORT extern void rt_mesh_lod_cache_status_init(struct rt_mesh_lod_cache_status *status);
RT_EXPORT extern int db_mesh_lod_status(struct db_i *dbip, const char *name, struct rt_mesh_lod_cache_status *status);
RT_EXPORT extern int db_mesh_lod_refresh(struct db_i *dbip, const char *name, struct rt_mesh_lod_cache_status *status);
RT_EXPORT extern int db_mesh_lod_invalidate(struct db_i *dbip, const char *name, struct rt_mesh_lod_cache_status *status);
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
