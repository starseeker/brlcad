/*                    D R A W _ C A C H E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/draw_cache.h
 *
 * Persistent cache support for Obol drawing data.  The cache is owned by
 * libbrlobol because it stores view/drawing realization inputs, not core
 * librt raytrace state.
 */

#ifndef BRLOBOL_DRAW_CACHE_H
#define BRLOBOL_DRAW_CACHE_H

#include "common.h"

#include "brlobol/defines.h"

#include "vmath.h"

#include <stddef.h>
#include <stdint.h>

__BEGIN_DECLS

#define BRLOBOL_DRAW_CACHE_DIR ".DbDrawCache"
#define BRLOBOL_DRAW_CACHE_METADATA_NAME_MAX 128
#define BRLOBOL_DRAW_CACHE_METADATA_SHADER_MAX 256
#define BRLOBOL_DRAW_CACHE_PROXY_AABB 1
#define BRLOBOL_DRAW_CACHE_PROXY_OBB 2

struct db_i;
struct db_tree_state;
struct directory;

struct BRLObolDrawCacheStatus {
    int directoryFound;
    int hasCachedPayload;
    int clearedCacheEntry;
    int generatedCacheEntry;
};

struct BRLObolDrawProxyRecord {
    int kind;
    size_t pointCount;
    point_t points[8];
};

struct BRLObolDrawMetadataRecord {
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
    char materialName[BRLOBOL_DRAW_CACHE_METADATA_NAME_MAX];
    int hasShader;
    char shader[BRLOBOL_DRAW_CACHE_METADATA_SHADER_MAX];
};

/* A compact, cacheable structural occurrence used to rebuild the initial
 * drawing registry without recursively importing database objects. */
struct BRLObolDrawManifestOccurrence {
    char *path;
    char *sourceName;
    mat_t localMatrix;
    point_t boundsMin;
    point_t boundsMax;
    int booleanOperation;
    uint32_t occurrenceIndex;
    int metadataValid;
    struct BRLObolDrawMetadataRecord metadata;
};

struct BRLObolDrawManifest {
    size_t occurrenceCount;
    struct BRLObolDrawManifestOccurrence *occurrences;
};

BRLOBOL_EXPORT void
brlobol_draw_cache_status_init(struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT void
brlobol_draw_proxy_record_init(struct BRLObolDrawProxyRecord *record);

BRLOBOL_EXPORT void
brlobol_draw_metadata_record_init(struct BRLObolDrawMetadataRecord *record);

BRLOBOL_EXPORT void
brlobol_draw_manifest_init(struct BRLObolDrawManifest *manifest);

BRLOBOL_EXPORT void
brlobol_draw_manifest_free(struct BRLObolDrawManifest *manifest);

BRLOBOL_EXPORT int
brlobol_draw_manifest_cache_store(struct db_i *dbip, const char *rootPath,
	const struct BRLObolDrawManifest *manifest);

/* manifest must first be initialized with brlobol_draw_manifest_init.  A
 * successful call replaces its contents; callers release them with free. */
BRLOBOL_EXPORT int
brlobol_draw_manifest_cache_get(struct db_i *dbip, const char *rootPath,
	struct BRLObolDrawManifest *manifest);

/* A manifest describes an entire draw-root hierarchy, so an edit to any
 * object conservatively invalidates all manifests for that database. */
BRLOBOL_EXPORT int
brlobol_draw_manifest_cache_invalidate_database(struct db_i *dbip);

BRLOBOL_EXPORT void
brlobol_draw_metadata_record_from_tree_state(
    struct BRLObolDrawMetadataRecord *record,
    const struct db_tree_state *state,
    const struct directory *dp);

BRLOBOL_EXPORT void
brlobol_draw_cache_clear_all(void);

BRLOBOL_EXPORT int
brlobol_draw_cache_clear_database(struct db_i *dbip);

BRLOBOL_EXPORT int
brlobol_draw_proxy_cache_status(struct db_i *dbip,
				const char *name,
				int kind,
				struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_proxy_cache_refresh(struct db_i *dbip,
				 const char *name,
				 int kind,
				 struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_proxy_cache_invalidate(struct db_i *dbip,
				    const char *name,
				    int kind,
				    struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_proxy_cache_store(struct db_i *dbip,
			       const char *name,
			       int kind,
			       const point_t *points,
			       size_t pointCount,
			       struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_proxy_cache_get(struct db_i *dbip,
			     const char *name,
			     int kind,
			     struct BRLObolDrawProxyRecord *record);

BRLOBOL_EXPORT int
brlobol_draw_metadata_cache_status(struct db_i *dbip,
				   const char *name,
				   struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_metadata_cache_refresh(struct db_i *dbip,
				    const char *name,
				    struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_metadata_cache_invalidate(struct db_i *dbip,
				       const char *name,
				       struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_metadata_cache_store(struct db_i *dbip,
				  const char *name,
				  const struct BRLObolDrawMetadataRecord *record,
				  struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_metadata_cache_get(struct db_i *dbip,
				const char *name,
				struct BRLObolDrawMetadataRecord *record);

BRLOBOL_EXPORT int
brlobol_draw_path_metadata_cache_status(struct db_i *dbip,
					const char *path,
					struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_path_metadata_cache_refresh(struct db_i *dbip,
					 const char *path,
					 struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_path_metadata_cache_invalidate(struct db_i *dbip,
					    const char *path,
					    struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_path_metadata_cache_invalidate_object(
					    struct db_i *dbip,
					    const char *name,
					    struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_path_metadata_cache_store(struct db_i *dbip,
				       const char *path,
				       const struct BRLObolDrawMetadataRecord *record,
				       struct BRLObolDrawCacheStatus *status);

BRLOBOL_EXPORT int
brlobol_draw_path_metadata_cache_get(struct db_i *dbip,
				     const char *path,
				     struct BRLObolDrawMetadataRecord *record);

__END_DECLS

#endif /* BRLOBOL_DRAW_CACHE_H */
