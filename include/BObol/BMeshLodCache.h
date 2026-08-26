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
#include <stdint.h>

__BEGIN_DECLS

struct db_i;
struct BObolMeshLod;
struct rt_bot_internal;

#define BOBOL_MESH_LOD_QUANTIZATION_BITS 16
#define BOBOL_MESH_LOD_CUT_COUNT_MAX \
    (1 + 3 * (BOBOL_MESH_LOD_QUANTIZATION_BITS - 1))

/* Display-provider contract version.  Bump this whenever the meaning or
 * completeness of a selected PoP cut changes. */
#define BOBOL_MESH_LOD_PROVIDER_VERSION "bobol-chunked-pop-v2"

/* Large logical leaves are divided into deterministic private cache pages.
 * This is a latency/working-set bound, not a visible tessellation rule. */
#define BOBOL_MESH_LOD_CHUNK_FACE_TARGET 65536

/* Spatial subresources are deliberately not display objects.  Sparse
 * occupied-cell records describe triangle-list ranges within the one
 * activation-ordered PoP stream owned by a logical leaf. */
#define BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION 8
#define BOBOL_MESH_LOD_CLUSTER_COUNT \
    (BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION * \
     BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION * \
     BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION)

struct BObolMeshLodClusterRange {
    uint32_t first_index;
    uint32_t index_count;
    uint8_t activation_cut;
};

struct BObolMeshLodClusterInfo {
    /* Uniform-grid cell number.  Only occupied cells are published, in
     * ascending order; an empty cell consumes no memory or cache record. */
    uint32_t cluster_id;
    point_t bmin;
    point_t bmax;
    const struct BObolMeshLodClusterRange *ranges;
    uint32_t range_count;
};

/* One independently resident private render subresource.  Chunks retain the
 * logical leaf's common cut schedule and quantization domain.  Counts are
 * cumulative within the chunk.  The records are borrowed from an opened
 * immutable cache handle. */
struct BObolMeshLodChunkCutInfo {
    uint32_t face_count;
    uint32_t point_count;
    uint64_t resident_bytes;
};

struct BObolMeshLodChunkInfo {
    uint32_t chunk_id;
    int min_cut;
    int max_cut;
    point_t bmin;
    point_t bmax;
    struct BObolMeshLodChunkCutInfo cuts[BOBOL_MESH_LOD_CUT_COUNT_MAX];
};

/* One borrowed chunk prefix.  Face indices are local to points and normals,
 * when present, contain three authored corner normals per face.  The pointers
 * are valid only for the duration of the callback. */
typedef int (*BObolMeshLodChunkCallback)(
	uint32_t chunk_id,
	int cut,
	const point_t *points,
	size_t point_count,
	const uint32_t *faces,
	size_t face_count,
	const vect_t *normals,
	size_t normal_count,
	void *cb_data);

/* Borrowed active LoD arrays; valid until the LoD is reloaded or destroyed. */
struct BObolMeshLodData {
    const uint32_t *faces;
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

/* One borrowed activation-cut suffix from an immutable PoP cache snapshot.
 * Points and faces contain only the records introduced at cut; face indices
 * address the cumulative activation-ordered point array.  normals is either
 * NULL or contains three corner normals per face.  All pointers are valid only
 * for the duration of the callback. */
typedef int (*BObolMeshLodSuffixCallback)(
	int cut,
	const point_t *points,
	size_t point_count,
	const uint32_t *faces,
	size_t face_count,
	const vect_t *normals,
	size_t normal_count,
	void *cb_data);

/* Stable summary of the active mesh LoD state. */
struct BObolMeshLodInfo {
    int active_cut;
    size_t face_count;
    size_t point_count;
    size_t point_orig_count;
    size_t normal_count;
    int has_faces;
    int has_points;
    int has_original_points;
    int has_snapped_points;
    int has_normals;
    int shaded_cull_backfaces;
    point_t bmin;
    point_t bmax;
};

/* Immutable facts for one independently drawable cumulative prefix. */
struct BObolMeshLodCutInfo {
    uint64_t face_count;
    uint64_t point_count;
    uint64_t resident_bytes;
    double object_error;
    uint8_t quantization_bits[3];
    uint8_t exact;
};

/* Immutable hierarchy metadata plus the currently resident cumulative cut.
 * Counts are cumulative: face_count[n] and point_count[n] describe the prefix
 * that may be drawn at cut n. */
struct BObolMeshLodHierarchyInfo {
    int min_cut;
    int max_cut;
    int resident_cut;
    uint32_t cut_count;
    /* Immutable representation facts needed before any prefix is loaded.
     * Resident-memory admission must not infer these from whichever arrays
     * happened to be active on an opened cache handle. */
    int has_normals;
    int shaded_cull_backfaces;
    /* The exact, tight quantization domain used when classifying and snapping
     * PoP vertices.  A retained renderer must use these values rather than
     * reconstructing them from separately rounded display bounds.  A
     * zero-extent axis is left unsnapped. */
    point_t quantization_min;
    point_t quantization_max;
    uint16_t cluster_grid_resolution;
    uint32_t cluster_count;
    /* Borrowed from the opened immutable cache handle. */
    const struct BObolMeshLodClusterInfo *clusters;
    uint32_t chunk_count;
    /* Borrowed from the opened immutable cache handle.  Unlike clusters,
     * these records describe independently readable storage. */
    const struct BObolMeshLodChunkInfo *chunks;
    struct BObolMeshLodCutInfo cuts[BOBOL_MESH_LOD_CUT_COUNT_MAX];
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

/* Database-wide cache coverage.  This is intentionally a name-map summary,
 * not a promise that every payload byte has been audited.  Cache generation
 * commits an immutable payload before publishing its name mapping, so the
 * summary is the inexpensive normal-operation readiness test.  Use the
 * per-object status API when diagnosing corruption or incomplete writes. */
struct BObolMeshLodCacheSummary {
    uint64_t database_bot_count;
    uint64_t mapped_bot_count;
    uint64_t missing_bot_count;
    int all_bots_mapped;
};

enum BObolMeshLodPreviewKind {
    /* A complete cumulative PoP prefix.  Its hierarchy metadata is valid for
     * normal mesh presentation, although durable cache work may continue. */
    BOBOL_MESH_LOD_PREVIEW_MESH_PREFIX = 1,
    /* A bounded point representation stratified over the whole serialized
     * source extent.  It is an overview representation, not a PoP cut: no
     * face data or hierarchy-level conclusion may be inferred from it. */
    BOBOL_MESH_LOD_PREVIEW_COVERAGE_POINTS = 2
};

enum {
    /* The coverage preview is a bounded visual occupancy summary.  It is not
     * cache hierarchy storage and never changes the durable PoP format. */
    BOBOL_MESH_LOD_COVERAGE_PREVIEW_CELL_AXIS = 24
};

/* A true-cold large mesh may spend substantial time constructing private
 * spatial pages and persisting them.  This callback borrows either a complete
 * PoP prefix or a separately identified coverage preview for the duration of
 * the call, allowing a background provider to publish useful content without
 * delaying durable cache work.  Implementations must copy retained data. */
typedef void (*BObolMeshLodPreviewCallback)(
	int preview_kind,
	unsigned long long cache_key,
	const struct BObolMeshLodData *data,
	const struct BObolMeshLodHierarchyInfo *hierarchy,
	void *callback_data);

/* One complete, independently drawable spatial cache page.  The page arrays
 * and hierarchy snapshot are borrowed for the duration of the callback.  A
 * receiver that keeps a page must make an immutable copy before returning;
 * the producer never exposes a partially written page. */
struct BObolMeshLodSpatialPage {
    uint32_t page_id;
    int cut;
    struct BObolMeshLodChunkInfo info;
    struct BObolMeshLodData data;
    struct BObolMeshLodHierarchyInfo hierarchy;
};

typedef void (*BObolMeshLodSpatialPageCallback)(
	unsigned long long cache_key,
	const struct BObolMeshLodSpatialPage *page,
	void *callback_data);

/* Cold cache construction can outlive the first useful presentation.  The
 * producer polls this optional callback at bounded source/page intervals so
 * an obsolete view may stop its private work without asynchronous teardown.
 * It runs on the source worker; callback_data must remain valid for the
 * synchronous refresh call and the callback must not re-enter this API. */
typedef int (*BObolMeshLodCancellationCallback)(void *callback_data);

struct BObolMeshLodPreviewRequest {
    int requested_cut;
    float projected_pixel_diameter;
    float target_pixel_error;
    /* Select the bounded serialized spatial producer for a large cold source.
     * This is request policy, not an environment-controlled cache format. */
    int spatial_leaf_producer;
    /* The caller's already-discovered local-space source extent permits a
     * coverage preview to avoid a duplicate full vertex scan.  Durable cache
     * generation still validates and persists independently scanned bounds. */
    int coverage_bounds_valid;
    point_t coverage_bmin;
    point_t coverage_bmax;
    /* Called only for complete, locally validated pages produced by the
     * serialized spatial-cache path.  It is distinct from the whole-source
     * preview callback above: a page never claims to cover the logical BoT. */
    BObolMeshLodSpatialPageCallback spatial_page_callback;
    void *spatial_page_data;
    BObolMeshLodCancellationCallback cancellation_callback;
    void *cancellation_data;
};

#define BOBOL_MESH_LOD_PREVIEW_REQUEST_INIT \
    { -1, 0.0f, 0.0f, 0, 0, VINIT_ZERO, VINIT_ZERO, NULL, NULL, NULL, NULL }

#define BOBOL_MESH_LOD_INFO_INIT { -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, VINIT_ZERO, VINIT_ZERO }
#define BOBOL_MESH_LOD_HIERARCHY_INFO_INIT { -1, -1, -1, 0, 0, 0, VINIT_ZERO, VINIT_ZERO, 0, 0, NULL, 0, NULL, {{0, 0, 0, 0.0, {0, 0, 0}, 0}} }
#define BOBOL_MESH_LOD_CACHE_STATUS_INIT { 0, 0, 0, 0, 0, 0, 0, 0, 0 }
#define BOBOL_MESH_LOD_CACHE_SUMMARY_INIT { 0, 0, 0, 0 }

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

/* Check all database BoT name mappings using one cache read snapshot.  This
 * is O(number of BoTs), but avoids one LMDB transaction and full hierarchy
 * materialization per object. */
BOBOL_EXPORT int
bobol_mesh_lod_cache_summary(struct db_i *dbip,
			       struct BObolMeshLodCacheSummary *summary);

BOBOL_EXPORT int
bobol_mesh_lod_cache_refresh(struct db_i *dbip,
	const char *name,
	struct BObolMeshLodCacheStatus *status);

/* Generate a missing cache from bot when supplied, or import name from dbip
 * otherwise, and return its opened prefix.  For a multi-page mesh preview is
 * called once after global prefixes are classified but before page
 * construction and persistence complete.  preview_request carries the same
 * ordinal/screen-error demand used by ordinary warm selection; generation
 * clamps it to the hierarchy and a transient-memory ceiling.
 */
BOBOL_EXPORT struct BObolMeshLod *
bobol_mesh_lod_cache_refresh_open(
	struct db_i *dbip,
	const char *name,
	const struct rt_bot_internal *bot,
	struct BObolMeshLodCacheStatus *status,
	const struct BObolMeshLodPreviewRequest *preview_request,
	BObolMeshLodPreviewCallback preview,
	void *preview_data);

/* Generate a missing V5 authored-BoT cache by decoding its serialized point
 * and face records directly.  This display-specific fallback avoids librt's
 * bu_bomb allocation semantics, validates all source records, and reports an
 * ordinary cache miss on capacity refusal.  It is intentionally not a
 * general replacement for rt_db_get_internal. */
BOBOL_EXPORT struct BObolMeshLod *
bobol_mesh_lod_cache_refresh_serialized_open(
	struct db_i *dbip,
	const char *name,
	struct BObolMeshLodCacheStatus *status,
	const struct BObolMeshLodPreviewRequest *preview_request,
	BObolMeshLodPreviewCallback preview,
	void *preview_data);

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
	int shaded_cull_backfaces,
	struct BObolMeshLodCacheStatus *status);

/* Store a representation variant without deleting another payload associated
 * with the same database object name.  Consumers must retain and reopen the
 * returned immutable cache key.  This is the BREP/tessellation-band path;
 * authored BoT replacement continues to use store_mesh above. */
BOBOL_EXPORT int
bobol_mesh_lod_cache_store_mesh_variant(
	struct db_i *dbip,
	const char *name,
	const point_t *vertices,
	size_t vertex_count,
	const vect_t *normals,
	const int *faces,
	size_t face_count,
	unsigned long long user_key,
	int shaded_cull_backfaces,
	struct BObolMeshLodCacheStatus *status);

BOBOL_EXPORT struct BObolMeshLod *
bobol_mesh_lod_get(struct db_i *dbip, const char *name);

/* Open the immutable PoP prefix currently associated with a name without
 * resolving or importing the database directory object.  This is the
 * retained-display path: the submitting source has already validated the
 * database object and its revision, and a stale/missing name mapping simply
 * returns NULL.  Unlike bobol_mesh_lod_get(), the returned handle deliberately
 * has no source-detail callbacks. */
BOBOL_EXPORT struct BObolMeshLod *
bobol_mesh_lod_get_named_cached_prefix(struct db_i *dbip,
	const char *name);

/* Open an already validated cache payload directly by immutable content key.
 * The hierarchy header and the immediately following resident-prefix load
 * share one read snapshot.  This is the high-throughput retained drawing path
 * after a source request has captured BObolMeshLodCacheStatus::cache_key. */
BOBOL_EXPORT struct BObolMeshLod *
bobol_mesh_lod_get_cached_prefix(struct db_i *dbip,
	unsigned long long cache_key);

/* Return the immutable cache identity backing an opened PoP handle.  This
 * lets warm-cache consumers avoid a second status/name-key transaction after
 * bobol_mesh_lod_get has already validated the same payload. */
BOBOL_EXPORT unsigned long long
bobol_mesh_lod_cache_key_get(const struct BObolMeshLod *lod);

/* Materialize exactly one producer-defined cumulative cut.  The requested
 * cut must lie in [min_cut, max_cut]; callers select it from hierarchy
 * metadata rather than relying on ordinal clamping. */
BOBOL_EXPORT int
bobol_mesh_lod_load_cut(struct BObolMeshLod *lod,
			  int cut,
			    int reset);

/* Append/trim exactly one resident cumulative prefix without materializing the
 * cut-snapped CPU point array.  Retained renderers consume points_orig and
 * select a draw prefix independently for each occurrence.  The same strict
 * [min_cut, max_cut] input contract applies. */
BOBOL_EXPORT int
bobol_mesh_lod_load_resident_cut(struct BObolMeshLod *lod,
				   int cut,
				     int reset);

/* Read only the activation records in (resident_cut, target_cut] from one
 * cache transaction.  Unlike bobol_mesh_lod_load_resident_cut(), this does
 * not materialize or retain a duplicate cumulative prefix in lod and does not
 * change its active cut.  It is the immutable-renderer append path after a
 * stable memory compaction has released the cache reader's working arrays. */
BOBOL_EXPORT int
bobol_mesh_lod_read_resident_suffix(struct BObolMeshLod *lod,
	int resident_cut,
	int target_cut,
	BObolMeshLodSuffixCallback callback,
	void *cb_data);

/* Read independently resident prefixes for the requested private chunks.
 * No unrequested chunk arrays are materialized.  chunk_ids must be strictly
 * increasing and each requested cut must be in the logical hierarchy range.
 * The callback is invoked once per chunk in input order. */
BOBOL_EXPORT int
bobol_mesh_lod_read_chunk_prefixes(struct BObolMeshLod *lod,
	const uint32_t *chunk_ids,
	size_t chunk_count,
	int cut,
	BObolMeshLodChunkCallback callback,
	void *cb_data);

/* Open an independent read-only handle for the same immutable cached
 * hierarchy.  The peer owns its cache transaction and active-prefix state,
 * allowing bounded callers to read disjoint private-page ranges in parallel.
 * It must be released with bobol_mesh_lod_destroy(). */
BOBOL_EXPORT struct BObolMeshLod *
bobol_mesh_lod_clone_reader(const struct BObolMeshLod *lod);

BOBOL_EXPORT int
bobol_mesh_lod_current_cut(const struct BObolMeshLod *lod);

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
bobol_mesh_lod_hierarchy_info_init(
	struct BObolMeshLodHierarchyInfo *info);

BOBOL_EXPORT int
bobol_mesh_lod_hierarchy_info_get(
	const struct BObolMeshLod *lod,
	struct BObolMeshLodHierarchyInfo *info);

/* True when this live asset retains only a bounded prefix because durable
 * cache publication was unavailable in its current capacity epoch.  It is
 * drawable and terminal at its advertised maximum cut, but the caller must
 * report the capacity limitation rather than mistake it for exact geometry. */
BOBOL_EXPORT int
bobol_mesh_lod_source_limited(const struct BObolMeshLod *lod);

/* Select the coarsest cut whose conservative object-space quantization error
 * projects to no more than target_pixel_error.  projected_pixel_diameter is
 * the occurrence's full projected bound; selection is therefore invariant to
 * the number and ordinal spacing of producer-defined cuts. */
BOBOL_EXPORT int
bobol_mesh_lod_select_cut(
	const struct BObolMeshLodHierarchyInfo *hierarchy,
	double projected_pixel_diameter,
	double target_pixel_error);

BOBOL_EXPORT void
bobol_mesh_lod_memshrink(struct BObolMeshLod *lod);

/* Bytes owned by the opened PoP handle.  The total includes fixed hierarchy
 * state and all retained prefix arrays; prefix_bytes reports only the
 * reloadable active arrays released by bobol_mesh_lod_memshrink(). */
BOBOL_EXPORT size_t
bobol_mesh_lod_resident_bytes(const struct BObolMeshLod *lod);

BOBOL_EXPORT size_t
bobol_mesh_lod_resident_prefix_bytes(const struct BObolMeshLod *lod);

BOBOL_EXPORT void
bobol_mesh_lod_destroy(struct BObolMeshLod *lod);

__END_DECLS

#endif /* BOBOL_BMESHLODCACHE_H */
