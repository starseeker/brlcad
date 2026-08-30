/*                  B D R A W C A C H E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BDrawCache.h
 *
 * Persistent cache support for Obol drawing data.  The cache is owned by
 * libBObol because it stores view/drawing realization inputs, not core
 * librt raytrace state.
 */

#ifndef BOBOL_BDRAWCACHE_H
#define BOBOL_BDRAWCACHE_H

#include "common.h"

#include "BObol/BDefines.h"

#include "vmath.h"

#include <stddef.h>
#include <stdint.h>

__BEGIN_DECLS

#define BOBOL_DRAW_CACHE_DIR ".DbDrawCache"
#define BOBOL_DRAW_CACHE_METADATA_NAME_MAX 128
#define BOBOL_DRAW_CACHE_METADATA_SHADER_MAX 256
#define BOBOL_DRAW_CACHE_LOD_ASSET_NAME_MAX 1024
#define BOBOL_DRAW_CACHE_PROXY_AABB 1
#define BOBOL_DRAW_CACHE_PROXY_OBB 2
#define BOBOL_DRAW_CACHE_MESH_ASSET_UNKNOWN 0
#define BOBOL_DRAW_CACHE_MESH_ASSET_BOT 1
#define BOBOL_DRAW_CACHE_MESH_ASSET_BREP 2

struct db_i;
struct db_tree_state;
struct directory;

struct BObolDrawCacheStatus {
    int directoryFound;
    int hasCachedPayload;
    int clearedCacheEntry;
    int generatedCacheEntry;
};

struct BObolDrawProxyRecord {
    int kind;
    size_t pointCount;
    point_t points[8];
};

/* A validated mapping from one baked database mesh to the canonical mesh
 * whose progressive payload it reuses.  assetToObject maps canonical asset
 * coordinates into this object's native coordinates. */
struct BObolDrawLodAssetRecord {
    char assetName[BOBOL_DRAW_CACHE_LOD_ASSET_NAME_MAX];
    uint64_t faceCount;
    uint64_t pointCount;
    point_t boundsMin;
    point_t boundsMax;
    point_t assetBoundsMin;
    point_t assetBoundsMax;
    mat_t assetToObject;
    /* Optional whole-asset PCA bounds in canonical asset coordinates.
     * These are presentation metadata for a terminal aggregate proxy; the
     * AABB fields above remain the authoritative coverage contract. */
    int assetOrientedBoundsValid;
    point_t assetOrientedBounds[8];
};

struct BObolDrawMetadataRecord {
    int directoryFound;
    int isPhony;
    int flags;
    int majorType;
    int minorType;
    int isSolid;
    int isComb;
    int isRegion;
    int isHidden;
    int hasRegionId;
    int regionId;
    int hasAircode;
    int aircode;
    int hasLos;
    int los;
    int hasMaterialId;
    int materialId;
    int hasInherit;
    int inherit;
    int hasColor;
    unsigned char color[3];
    int hasMaterialName;
    char materialName[BOBOL_DRAW_CACHE_METADATA_NAME_MAX];
    int hasShader;
    char shader[BOBOL_DRAW_CACHE_METADATA_SHADER_MAX];
};

/* A compact, cacheable structural occurrence used to rebuild the initial
 * drawing registry without recursively importing database objects. */
struct BObolDrawManifestOccurrence {
    char *path;
    char *sourceName;
    /* A completed leaf manifest carries enough immutable source-asset
     * identity to start view-LoD immediately.  Without this contract a warm
     * draw can publish every cached leaf box, but still has to repeat the
     * database hierarchy/import walk before any box is eligible for PoP
     * replacement.  Structural-only manifests leave sourceMeshRequestValid
     * false. */
    int sourceMeshRequestValid;
    int meshAssetKind;
    uint64_t meshAssetContentHash;
    double meshAssetTessellationAbsTol;
    double meshAssetTessellationRelTol;
    double meshAssetTessellationNormTol;
    char *meshAssetPath;
    char *meshAssetName;
    point_t meshAssetBoundsMin;
    point_t meshAssetBoundsMax;
    mat_t meshAssetMatrix;
    uint64_t sourceFaceCount;
    uint64_t sourcePointCount;
    mat_t localMatrix;
    point_t boundsMin;
    point_t boundsMax;
    /* Optional object-coordinate PCA proxy corners.  The AABB remains the
     * conservative coverage and view-framing contract. */
    int orientedBoundsValid;
    point_t orientedBounds[8];
    int booleanOperation;
    uint32_t occurrenceIndex;
    int metadataValid;
    struct BObolDrawMetadataRecord metadata;
};

struct BObolDrawManifest {
    /* Exact source-local extent of the complete draw target.  This is a
     * path-scoped fact, unlike the object-name proxy cache: a manifest for
     * a/b/c may include transforms and Boolean context which are not valid
     * for another occurrence of c. */
    int coverageBoundsValid;
    point_t coverageBoundsMin;
    point_t coverageBoundsMax;
    /* Complete, mesh-buffer-free population proof captured by cold discovery.
     * A warm consumer may use these aggregate facts to select its bounded
     * admission policy before occurrence chunks finish streaming.  Zero
     * values mean that this is a structural-only manifest. */
    uint64_t uniqueAssetCount;
    uint64_t encodedSourceBytes;
    uint64_t largestAssetBytes;
    size_t occurrenceCount;
    struct BObolDrawManifestOccurrence *occurrences;
};

BOBOL_EXPORT void
bobol_draw_cache_status_init(struct BObolDrawCacheStatus *status);

BOBOL_EXPORT void
bobol_draw_proxy_record_init(struct BObolDrawProxyRecord *record);

BOBOL_EXPORT void
bobol_draw_lod_asset_record_init(struct BObolDrawLodAssetRecord *record);

BOBOL_EXPORT int
bobol_draw_lod_asset_cache_store(struct db_i *dbip,
	const char *name,
	const struct BObolDrawLodAssetRecord *record);

BOBOL_EXPORT int
bobol_draw_lod_asset_cache_get(struct db_i *dbip,
	const char *name,
	struct BObolDrawLodAssetRecord *record);

BOBOL_EXPORT void
bobol_draw_metadata_record_init(struct BObolDrawMetadataRecord *record);

BOBOL_EXPORT void
bobol_draw_manifest_init(struct BObolDrawManifest *manifest);

BOBOL_EXPORT void
bobol_draw_manifest_free(struct BObolDrawManifest *manifest);

/* Supplies one borrowed occurrence for manifest persistence.  The store may
 * call this provider more than once for an index; pointers are needed only
 * until the callback returns. */
typedef int (*BObolDrawManifestOccurrenceProvider)(size_t occurrenceIndex,
    struct BObolDrawManifestOccurrence *occurrence, void *userData);

BOBOL_EXPORT int
bobol_draw_manifest_cache_store(struct db_i *dbip, const char *rootPath,
	const struct BObolDrawManifest *manifest);

/* Store a manifest from a bounded caller-owned occurrence source.  The
 * manifest supplies only the count and exact coverage extent; occurrences is
 * not read.  This avoids a second full C occurrence array during cold cache
 * creation for very large scenes. */
BOBOL_EXPORT int
bobol_draw_manifest_cache_store_visit(struct db_i *dbip,
    const char *rootPath, const struct BObolDrawManifest *manifest,
    BObolDrawManifestOccurrenceProvider occurrence, void *userData);

/* manifest must first be initialized with bobol_draw_manifest_init.  A
 * successful call replaces its contents; callers release them with free. */
BOBOL_EXPORT int
bobol_draw_manifest_cache_get(struct db_i *dbip, const char *rootPath,
	struct BObolDrawManifest *manifest);

/* Visit a manifest directly from its read-only cache mapping.  begin is
 * called once after the header and exact coverage extent are validated;
 * occurrence is then called in stored order.  Occurrence string pointers are
 * valid only for the duration of that callback.  Returning zero from either
 * callback stops the visit and reports BRLCAD_ERROR.  This API avoids copying
 * and retaining an entire large manifest before its first useful record can
 * be published. */
typedef int (*BObolDrawManifestBeginCallback)(
    const struct BObolDrawManifest *manifest, void *userData);
typedef int (*BObolDrawManifestOccurrenceCallback)(
    const struct BObolDrawManifestOccurrence *occurrence,
    size_t occurrenceIndex, void *userData);

BOBOL_EXPORT int
bobol_draw_manifest_cache_stream(struct db_i *dbip, const char *rootPath,
    BObolDrawManifestBeginCallback begin,
    BObolDrawManifestOccurrenceCallback occurrence, void *userData);

/* Read only the validated manifest description and exact coverage extent.
 * No occurrence strings or records are decoded or allocated, making this
 * suitable for an O(1) warm-start overview on the UI thread.  The returned
 * description never owns an occurrences array. */
BOBOL_EXPORT int
bobol_draw_manifest_cache_describe(struct db_i *dbip, const char *rootPath,
    struct BObolDrawManifest *description);

/* Manifests and transformed-LoD asset mappings describe relationships among
 * database objects, so an edit conservatively invalidates both record types. */
BOBOL_EXPORT int
bobol_draw_manifest_cache_invalidate_database(struct db_i *dbip);

BOBOL_EXPORT void
bobol_draw_metadata_record_from_tree_state(
    struct BObolDrawMetadataRecord *record,
    const struct db_tree_state *state,
    const struct directory *dp);

BOBOL_EXPORT void
bobol_draw_cache_clear_all(void);

BOBOL_EXPORT int
bobol_draw_cache_clear_database(struct db_i *dbip);

BOBOL_EXPORT int
bobol_draw_proxy_cache_status(struct db_i *dbip,
				const char *name,
				int kind,
				struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_proxy_cache_refresh(struct db_i *dbip,
				 const char *name,
				 int kind,
				 struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_proxy_cache_invalidate(struct db_i *dbip,
				    const char *name,
				    int kind,
				    struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_proxy_cache_store(struct db_i *dbip,
			       const char *name,
			       int kind,
			       const point_t *points,
			       size_t pointCount,
			       struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_proxy_cache_get(struct db_i *dbip,
			     const char *name,
			     int kind,
			     struct BObolDrawProxyRecord *record);

BOBOL_EXPORT int
bobol_draw_metadata_cache_status(struct db_i *dbip,
				   const char *name,
				   struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_metadata_cache_refresh(struct db_i *dbip,
				    const char *name,
				    struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_metadata_cache_invalidate(struct db_i *dbip,
				       const char *name,
				       struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_metadata_cache_store(struct db_i *dbip,
				  const char *name,
				  const struct BObolDrawMetadataRecord *record,
				  struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_metadata_cache_get(struct db_i *dbip,
				const char *name,
				struct BObolDrawMetadataRecord *record);

BOBOL_EXPORT int
bobol_draw_path_metadata_cache_status(struct db_i *dbip,
					const char *path,
					struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_path_metadata_cache_refresh(struct db_i *dbip,
					 const char *path,
					 struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_path_metadata_cache_invalidate(struct db_i *dbip,
					    const char *path,
					    struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_path_metadata_cache_invalidate_object(
					    struct db_i *dbip,
					    const char *name,
					    struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_path_metadata_cache_store(struct db_i *dbip,
				       const char *path,
				       const struct BObolDrawMetadataRecord *record,
				       struct BObolDrawCacheStatus *status);

BOBOL_EXPORT int
bobol_draw_path_metadata_cache_get(struct db_i *dbip,
				     const char *path,
				     struct BObolDrawMetadataRecord *record);

__END_DECLS

#endif /* BOBOL_BDRAWCACHE_H */
