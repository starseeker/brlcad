/*              M E S H _ L O D _ C A C H E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/mesh_lod_cache.h
 *
 * Persistent mesh LoD cache support for Obol drawing.  This API is owned by
 * libbrlobol because the cached mesh payloads are display realizations, not
 * core librt raytrace state.
 */

#ifndef BRLOBOL_MESH_LOD_CACHE_H
#define BRLOBOL_MESH_LOD_CACHE_H

#include "common.h"

#include "brlobol/defines.h"

#include "bv/view.h"
#include "vmath.h"

#include <stddef.h>

__BEGIN_DECLS

struct db_i;
struct BRLObolMeshLod;

/* Borrowed active LoD arrays; valid until the LoD is reloaded or destroyed. */
struct BRLObolMeshLodData {
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
 * the array lifetimes; libbrlobol borrows them until the matching clear/free
 * callback. */
struct BRLObolMeshLodDetail {
    const int *faces;
    size_t face_count;
    const point_t *points;
    size_t point_count;
    const point_t *points_orig;
    size_t point_orig_count;
    const vect_t *normals;
    size_t normal_count;
};

typedef int (*BRLObolMeshLodDetailSetupCallback)(
	struct BRLObolMeshLodDetail *detail, void *cb_data);
typedef int (*BRLObolMeshLodDetailClearCallback)(void *cb_data);
typedef int (*BRLObolMeshLodDetailFreeCallback)(void *cb_data);

/* Stable summary of the active mesh LoD state. */
struct BRLObolMeshLodInfo {
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

/* Stable status for a database object's mesh LoD cache entry. */
struct BRLObolMeshLodCacheStatus {
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

#define BRLOBOL_MESH_LOD_INFO_INIT { -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, VINIT_ZERO, VINIT_ZERO }
#define BRLOBOL_MESH_LOD_CACHE_STATUS_INIT { 0, 0, 0, 0, 0, 0, 0, 0, 0 }

BRLOBOL_EXPORT void
brlobol_mesh_lod_cache_init(struct db_i *dbip, int verbose);

BRLOBOL_EXPORT void
brlobol_mesh_lod_cache_clear_database(struct db_i *dbip);

BRLOBOL_EXPORT void
brlobol_mesh_lod_cache_clear_all(void);

BRLOBOL_EXPORT int
brlobol_mesh_lod_cache_update(struct db_i *dbip, const char *name);

BRLOBOL_EXPORT void
brlobol_mesh_lod_cache_status_init(
	struct BRLObolMeshLodCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_mesh_lod_cache_status(struct db_i *dbip,
			      const char *name,
			      struct BRLObolMeshLodCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_mesh_lod_cache_refresh(struct db_i *dbip,
			       const char *name,
			       struct BRLObolMeshLodCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_mesh_lod_cache_invalidate(struct db_i *dbip,
				  const char *name,
				  struct BRLObolMeshLodCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_mesh_lod_cache_store_mesh(
	struct db_i *dbip,
	const char *name,
	const point_t *vertices,
	size_t vertex_count,
	const vect_t *normals,
	const int *faces,
	size_t face_count,
	unsigned long long user_key,
	fastf_t fidelity_ratio,
	struct BRLObolMeshLodCacheStatus *status);

BRLOBOL_EXPORT struct BRLObolMeshLod *
brlobol_mesh_lod_get(struct db_i *dbip, const char *name);

BRLOBOL_EXPORT int
brlobol_mesh_lod_load_level(struct BRLObolMeshLod *lod,
			    int level,
			    int reset);

BRLOBOL_EXPORT int
brlobol_mesh_lod_load_view(struct BRLObolMeshLod *lod,
			   const struct bv_view_info *info,
			   int reset);

BRLOBOL_EXPORT int
brlobol_mesh_lod_current_level(const struct BRLObolMeshLod *lod);

BRLOBOL_EXPORT int
brlobol_mesh_lod_has_active_data(const struct BRLObolMeshLod *lod);

BRLOBOL_EXPORT int
brlobol_mesh_lod_data_get(const struct BRLObolMeshLod *lod,
			  struct BRLObolMeshLodData *data);

BRLOBOL_EXPORT void
brlobol_mesh_lod_info_init(struct BRLObolMeshLodInfo *info);

BRLOBOL_EXPORT int
brlobol_mesh_lod_info_get(const struct BRLObolMeshLod *lod,
			  struct BRLObolMeshLodInfo *info);

BRLOBOL_EXPORT void
brlobol_mesh_lod_detail_init(struct BRLObolMeshLodDetail *detail);

BRLOBOL_EXPORT int
brlobol_mesh_lod_detail_callbacks_set(
	struct BRLObolMeshLod *lod,
	BRLObolMeshLodDetailSetupCallback setup_clbk,
	BRLObolMeshLodDetailClearCallback clear_clbk,
	BRLObolMeshLodDetailFreeCallback free_clbk,
	void *cb_data);

BRLOBOL_EXPORT void
brlobol_mesh_lod_detail_callbacks_clear(struct BRLObolMeshLod *lod);

BRLOBOL_EXPORT void
brlobol_mesh_lod_memshrink(struct BRLObolMeshLod *lod);

BRLOBOL_EXPORT void
brlobol_mesh_lod_destroy(struct BRLObolMeshLod *lod);

__END_DECLS

#endif /* BRLOBOL_MESH_LOD_CACHE_H */
