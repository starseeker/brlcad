/*               B M E S H L O D C A C H E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BMeshLodCache.h
 *
 * Persistent mesh LoD cache support for Obol drawing.  This API is owned by
 * libBObol because the cached mesh payloads are display realizations, not
 * core librt raytrace state.
 */

#ifndef BOBOL_BMESHLODCACHE_H
#define BOBOL_BMESHLODCACHE_H

#include "common.h"

#include "BObol/BDefines.h"

#include "bv/view.h"
#include "vmath.h"

#include <stddef.h>

__BEGIN_DECLS

struct db_i;
struct BObolMeshLod;

/* Borrowed active LoD arrays; valid until the LoD is reloaded or destroyed. */
struct BObolMeshLodData {
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
 * the array lifetimes; libBObol borrows them until the matching clear/free
 * callback. */
struct BObolMeshLodDetail {
    const int *faces;
    size_t face_count;
    const point_t *points;
    size_t point_count;
    const point_t *points_orig;
    size_t point_orig_count;
    const vect_t *normals;
    size_t normal_count;
};

typedef int (*BObolMeshLodDetailSetupCallback)(
	struct BObolMeshLodDetail *detail, void *cb_data);
typedef int (*BObolMeshLodDetailClearCallback)(void *cb_data);
typedef int (*BObolMeshLodDetailFreeCallback)(void *cb_data);

/* Stable summary of the active mesh LoD state. */
struct BObolMeshLodInfo {
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
struct BObolMeshLodCacheStatus {
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

#define BOBOL_MESH_LOD_INFO_INIT { -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, VINIT_ZERO, VINIT_ZERO }
#define BOBOL_MESH_LOD_CACHE_STATUS_INIT { 0, 0, 0, 0, 0, 0, 0, 0, 0 }

BOBOL_EXPORT void
bobol_mesh_lod_cache_init(struct db_i *dbip, int verbose);

BOBOL_EXPORT void
bobol_mesh_lod_cache_clear_database(struct db_i *dbip);

BOBOL_EXPORT void
bobol_mesh_lod_cache_clear_all(void);

BOBOL_EXPORT int
bobol_mesh_lod_cache_update(struct db_i *dbip, const char *name);

BOBOL_EXPORT void
bobol_mesh_lod_cache_status_init(
	struct BObolMeshLodCacheStatus *status);

BOBOL_EXPORT int
bobol_mesh_lod_cache_status(struct db_i *dbip,
			      const char *name,
			      struct BObolMeshLodCacheStatus *status);

BOBOL_EXPORT int
bobol_mesh_lod_cache_refresh(struct db_i *dbip,
			       const char *name,
			       struct BObolMeshLodCacheStatus *status);

BOBOL_EXPORT int
bobol_mesh_lod_cache_invalidate(struct db_i *dbip,
				  const char *name,
				  struct BObolMeshLodCacheStatus *status);

BOBOL_EXPORT int
bobol_mesh_lod_cache_store_mesh(
	struct db_i *dbip,
	const char *name,
	const point_t *vertices,
	size_t vertex_count,
	const vect_t *normals,
	const int *faces,
	size_t face_count,
	unsigned long long user_key,
	fastf_t fidelity_ratio,
	struct BObolMeshLodCacheStatus *status);

BOBOL_EXPORT struct BObolMeshLod *
bobol_mesh_lod_get(struct db_i *dbip, const char *name);

BOBOL_EXPORT int
bobol_mesh_lod_load_level(struct BObolMeshLod *lod,
			    int level,
			    int reset);

BOBOL_EXPORT int
bobol_mesh_lod_load_view(struct BObolMeshLod *lod,
			   const struct bv_view_info *info,
			   int reset);

BOBOL_EXPORT int
bobol_mesh_lod_current_level(const struct BObolMeshLod *lod);

BOBOL_EXPORT int
bobol_mesh_lod_has_active_data(const struct BObolMeshLod *lod);

BOBOL_EXPORT int
bobol_mesh_lod_data_get(const struct BObolMeshLod *lod,
			  struct BObolMeshLodData *data);

BOBOL_EXPORT void
bobol_mesh_lod_info_init(struct BObolMeshLodInfo *info);

BOBOL_EXPORT int
bobol_mesh_lod_info_get(const struct BObolMeshLod *lod,
			  struct BObolMeshLodInfo *info);

BOBOL_EXPORT void
bobol_mesh_lod_detail_init(struct BObolMeshLodDetail *detail);

BOBOL_EXPORT int
bobol_mesh_lod_detail_callbacks_set(
	struct BObolMeshLod *lod,
	BObolMeshLodDetailSetupCallback setup_clbk,
	BObolMeshLodDetailClearCallback clear_clbk,
	BObolMeshLodDetailFreeCallback free_clbk,
	void *cb_data);

BOBOL_EXPORT void
bobol_mesh_lod_detail_callbacks_clear(struct BObolMeshLod *lod);

BOBOL_EXPORT void
bobol_mesh_lod_memshrink(struct BObolMeshLod *lod);

BOBOL_EXPORT void
bobol_mesh_lod_destroy(struct BObolMeshLod *lod);

__END_DECLS

#endif /* BOBOL_BMESHLODCACHE_H */
