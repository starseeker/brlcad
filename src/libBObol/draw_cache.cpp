/*                  D R A W _ C A C H E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file draw_cache.cpp */

#include "common.h"

#include "BObol/BDrawCache.h"
#include "BObol/BLodRealization.h"
#include "BObol/BMeshLodCache.h"

#include "draw_cache_private.h"

#include "raytrace.h"

#include "bu/app.h"
#include "bu/cache.h"
#include "bu/color.h"
#include "bu/file.h"
#include "bu/hash.h"
#include "bu/opt.h"
#include "bu/parallel.h"
#include "bu/path.h"
#include "bu/str.h"
#include "bu/vls.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include <cmath>
#include <array>
#include <map>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

#define BOBOL_DRAW_CACHE_FORMAT_FILE "draw_data.format"
/* Mesh LoD and draw records share BOBOL_DRAW_CACHE_DIR.  Keep the directory
 * format stable when only an independently versioned record changes; bumping
 * this value erases valid multi-gigabyte PoP caches as well as draw metadata. */
#define BOBOL_DRAW_CACHE_CURRENT_FORMAT 5
#define BOBOL_DRAW_CACHE_AABB "bb"
#define BOBOL_DRAW_CACHE_OBB "obb"
#define BOBOL_DRAW_CACHE_LOD_ASSET "lod_asset"
#define BOBOL_DRAW_CACHE_METADATA "meta"
#define BOBOL_DRAW_CACHE_PATH_METADATA "meta_path"
#define BOBOL_DRAW_CACHE_PATH_METADATA_INDEX "meta_path_index"
#define BOBOL_DRAW_CACHE_MANIFEST "manifest"
#define BOBOL_DRAW_METADATA_DISK_MAGIC 0x4f424d45u /* OBME */
#define BOBOL_DRAW_PROXY_DISK_MAGIC 0x4f425058u /* OBPX */
#define BOBOL_DRAW_PROXY_DISK_VERSION 1u
#define BOBOL_DRAW_LOD_ASSET_DISK_MAGIC 0x4f424c41u /* OBLA */
/* Version 3 adds optional canonical-asset OBB metadata.  Version 1 mappings
 * may also have been poisoned by an uncomposed reuse transform, so neither
 * older version is suitable for migration. */
#define BOBOL_DRAW_LOD_ASSET_DISK_VERSION 3u
#define BOBOL_DRAW_MANIFEST_CHUNKED_DISK_MAGIC 0x4f424d43u /* OBMC */
#define BOBOL_DRAW_MANIFEST_CHUNK_DISK_MAGIC 0x4f424d4bu /* OBMK */
#define BOBOL_DRAW_MANIFEST_CHUNK_DISK_VERSION 3u

/* Keep cache serialization bounded even for assemblies with hundreds of
 * thousands of leaves.  A descriptor is published only after every chunk has
 * been written, so an interrupted refresh is simply a cache miss. */
static const size_t manifestChunkTargetBytes = 8u * 1024u * 1024u;

struct BObolDrawCacheContext {
    bu_cache *cache;
    char *registryKey;
    db_i *validationDb;
    uint64_t validationFingerprint;
    std::unordered_set<std::string> *lodAssetKeys;
    std::unordered_set<std::string> *lodAssetOrientedBoundsKeys;
    bool lodAssetKeysLoaded;
};

struct BObolDrawCacheBinding {
    std::string filename;
    BObolDrawCacheContext *context;
};

struct BObolDrawCacheHandle {
    BObolDrawCacheContext *context;
    bu_cache *cache;
};

struct BObolDrawMetadataDiskRecord {
    uint32_t magic;
    uint32_t version;
    int32_t directoryFound;
    int32_t isPhony;
    int32_t flags;
    int32_t majorType;
    int32_t minorType;
    int32_t isSolid;
    int32_t isComb;
    int32_t isRegion;
    int32_t isHidden;
    int32_t hasRegionId;
    int32_t regionId;
    int32_t hasAircode;
    int32_t aircode;
    int32_t hasLos;
    int32_t los;
    int32_t hasMaterialId;
    int32_t materialId;
    int32_t hasInherit;
    int32_t inherit;
    int32_t hasColor;
    unsigned char color[3];
    unsigned char pad0;
    int32_t hasMaterialName;
    char materialName[BOBOL_DRAW_CACHE_METADATA_NAME_MAX];
    int32_t hasShader;
    char shader[BOBOL_DRAW_CACHE_METADATA_SHADER_MAX];
};

/* Cache keys are intentionally short hashes because bu_cache keys have a
 * fixed size.  The payload carries all identity fields so a hash collision or
 * a changed database object is always treated as a cache miss. */
struct BObolDrawProxyDiskHeader {
    uint32_t magic;
    uint32_t version;
    uint32_t kind;
    uint32_t pointCount;
    uint32_t nameLength;
    uint32_t directoryFlags;
    uint32_t directoryMajorType;
    uint32_t directoryMinorType;
    uint64_t databaseFingerprint;
    uint64_t directoryAddress;
    uint64_t directoryLength;
};

struct BObolDrawLodAssetDiskHeader {
    uint32_t magic;
    uint32_t version;
    uint32_t nameLength;
    uint32_t assetNameLength;
    uint32_t directoryFlags;
    uint32_t directoryMajorType;
    uint32_t directoryMinorType;
    uint32_t assetDirectoryFlags;
    uint32_t assetDirectoryMajorType;
    uint32_t assetDirectoryMinorType;
    uint32_t assetOrientedBoundsValid;
    uint64_t databaseFingerprint;
    uint64_t directoryAddress;
    uint64_t directoryLength;
    uint64_t assetDirectoryAddress;
    uint64_t assetDirectoryLength;
    uint64_t faceCount;
    uint64_t pointCount;
    point_t boundsMin;
    point_t boundsMax;
    point_t assetBoundsMin;
    point_t assetBoundsMax;
    mat_t assetToObject;
    point_t assetOrientedBounds[8];
};

struct BObolDrawManifestChunkedDiskHeader {
    uint32_t magic;
    uint32_t version;
    uint64_t databaseFingerprint;
    uint64_t occurrenceCount;
    uint64_t uniqueAssetCount;
    uint64_t encodedSourceBytes;
    uint64_t largestAssetBytes;
    uint32_t rootPathLength;
    uint32_t coverageBoundsValid;
    point_t coverageBoundsMin;
    point_t coverageBoundsMax;
    uint32_t chunkCount;
    uint32_t reserved;
};

struct BObolDrawManifestChunkDiskHeader {
    uint32_t magic;
    uint32_t version;
    uint64_t databaseFingerprint;
    uint32_t chunkIndex;
    uint32_t occurrenceCount;
};

struct BObolDrawManifestOccurrenceDiskHeader {
    uint32_t pathLength;
    uint32_t sourceNameLength;
    uint32_t meshAssetPathLength;
    uint32_t meshAssetNameLength;
    int32_t booleanOperation;
    uint32_t occurrenceIndex;
    uint32_t metadataValid;
    uint32_t sourceMeshRequestValid;
    uint32_t meshAssetKind;
    uint32_t orientedBoundsValid;
    uint64_t meshAssetContentHash;
    uint64_t sourceFaceCount;
    uint64_t sourcePointCount;
    double meshAssetTessellationAbsTol;
    double meshAssetTessellationRelTol;
    double meshAssetTessellationNormTol;
    fastf_t localMatrix[16];
    fastf_t meshAssetMatrix[16];
    point_t boundsMin;
    point_t boundsMax;
    point_t meshAssetBoundsMin;
    point_t meshAssetBoundsMax;
    BObolDrawMetadataDiskRecord metadata;
};

static int
bobol_draw_cache_semaphore(void)
{
    static int sem = 0;
    if (!sem)
	sem = bu_semaphore_register("BOBOL_DRAW_CACHE");
    return sem;
}

static std::mutex &
bobol_draw_cache_registry_mutex(void)
{
    static std::mutex m;
    return m;
}

static std::map<std::string, BObolDrawCacheContext *> &
bobol_draw_cache_registry(void)
{
    static std::map<std::string, BObolDrawCacheContext *> registry;
    return registry;
}

/* File-backed db_i instances issue large bursts of cache probes from coverage
 * workers.  Retain the already-resolved cache context for that live database
 * handle instead of repeating directory creation checks, path normalization,
 * and registry-key construction for every leaf.  In-memory databases are not
 * bound because an allocator may reuse their pointer with no filename by
 * which to distinguish the new database. */
static std::map<const db_i *, BObolDrawCacheBinding> &
bobol_draw_cache_bindings(void)
{
    static std::map<const db_i *, BObolDrawCacheBinding> bindings;
    return bindings;
}

static std::array<std::mutex, 127> &
bobol_draw_proxy_refresh_locks(void)
{
    static std::array<std::mutex, 127> refreshLocks;
    return refreshLocks;
}

/*
 * Source-realization workers may be active until the process-wide coordinator
 * is destroyed.  Its constructor calls this routine before starting workers,
 * which deliberately constructs every non-trivial function-local static a
 * cache operation can touch.  C++ then destroys the coordinator (and joins
 * its workers) before destroying these registries.  Do not make a new lazy
 * cache-global available to realization workers without adding it here.
 *
 * This is an internal lifetime barrier, not part of the public cache API.
 */
void
bobol_draw_cache_runtime_prepare(void)
{
    (void)bobol_draw_cache_registry_mutex();
    (void)bobol_draw_cache_registry();
    (void)bobol_draw_cache_bindings();
    (void)bobol_draw_proxy_refresh_locks();
}

/* A compact assembly may contain thousands of occurrences of one leaf.  The
 * persistent record is shared, but cold-cache refresh used to duplicate its
 * realization for every occurrence. */
static std::mutex &
bobol_draw_proxy_refresh_mutex(struct db_i *dbip, const char *name, int kind)
{
    std::array<std::mutex, 127> &refreshLocks =
	bobol_draw_proxy_refresh_locks();
    uint64_t hash = bu_data_hash(name ? name : "", name ? strlen(name) : 0);

    hash ^= static_cast<uint64_t>(reinterpret_cast<uintptr_t>(dbip));
    hash *= 1099511628211ULL;
    hash ^= static_cast<uint64_t>(kind);
    hash *= 1099511628211ULL;
    return refreshLocks[hash % refreshLocks.size()];
}

static void
bobol_draw_cache_context_close(BObolDrawCacheContext *context)
{
    if (!context)
	return;
    if (context->validationDb)
	db_close(context->validationDb);
    if (context->cache)
	bu_cache_close(context->cache);
    if (context->registryKey)
	bu_free(context->registryKey, "bobol draw cache registry key");
    delete context->lodAssetKeys;
    delete context->lodAssetOrientedBoundsKeys;
    BU_PUT(context, BObolDrawCacheContext);
}

static void
bobol_draw_cache_registry_close_all(void)
{
    std::map<std::string, BObolDrawCacheContext *> closeMap;

    {
	std::lock_guard<std::mutex> guard(bobol_draw_cache_registry_mutex());
	closeMap.swap(bobol_draw_cache_registry());
	bobol_draw_cache_bindings().clear();
    }

    for (std::map<std::string, BObolDrawCacheContext *>::iterator it =
	     closeMap.begin(); it != closeMap.end(); ++it)
	bobol_draw_cache_context_close(it->second);
}

static const char *
bobol_draw_proxy_component(int kind)
{
    switch (kind) {
	case BOBOL_LOD_PROXY_AABB:
	    return BOBOL_DRAW_CACHE_AABB;
	case BOBOL_LOD_PROXY_OBB:
	    return BOBOL_DRAW_CACHE_OBB;
	default:
	    return NULL;
    }
}

static size_t
bobol_draw_proxy_point_count(int kind)
{
    switch (kind) {
	case BOBOL_LOD_PROXY_AABB:
	    return 2;
	case BOBOL_LOD_PROXY_OBB:
	    return 8;
	default:
	    return 0;
    }
}

static int
bobol_draw_proxy_leaf_solid_directory(const directory *dp)
{
    if (!dp || dp == RT_DIR_NULL)
	return 0;
    if (dp->d_addr == RT_DIR_PHONY_ADDR)
	return 0;
    if (dp->d_flags & RT_DIR_COMB)
	return 0;
    return (dp->d_flags & RT_DIR_SOLID) ? 1 : 0;
}


static int
bobol_draw_proxy_leaf_solid_name(db_i *dbip, const char *name,
				   directory **dp_out)
{
    if (dp_out)
	*dp_out = RT_DIR_NULL;
    if (!dbip || !name || !name[0])
	return 0;

    directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp_out)
	*dp_out = dp;
    return bobol_draw_proxy_leaf_solid_directory(dp);
}

static int
bobol_draw_proxy_cache_directory(db_i *dbip, const char *name,
				  directory **dp_out)
{
    if (dp_out)
	*dp_out = RT_DIR_NULL;
    if (!dbip || !name || !name[0])
	return 0;

    directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp_out)
	*dp_out = dp;
    return dp != RT_DIR_NULL && dp->d_addr != RT_DIR_PHONY_ADDR &&
	!(dp->d_flags & RT_DIR_NON_GEOM) ? 1 : 0;
}

static uint64_t
bobol_draw_cache_mix(uint64_t hash, uint64_t value)
{
    hash ^= value;
    return hash * 1099511628211ULL;
}

static uint64_t
bobol_draw_cache_database_fingerprint(const db_i *dbip)
{
    if (!dbip)
	return 0;

    char canonical[MAXPATHLEN] = {0};
    const char *name = dbip->dbi_filename ? dbip->dbi_filename : "";
    if (name[0]) {
	(void)bu_file_realpath(name, canonical);
	name = canonical;
    }
    uint64_t fingerprint = 1469598103934665603ULL;
    fingerprint = bobol_draw_cache_mix(fingerprint,
	bu_data_hash(name, strlen(name)));
    struct stat sb;
    if (name[0] && stat(name, &sb) == 0) {
	fingerprint = bobol_draw_cache_mix(fingerprint,
	    static_cast<uint64_t>(sb.st_dev));
	fingerprint = bobol_draw_cache_mix(fingerprint,
	    static_cast<uint64_t>(sb.st_ino));
	fingerprint = bobol_draw_cache_mix(fingerprint,
	    static_cast<uint64_t>(sb.st_size));
	fingerprint = bobol_draw_cache_mix(fingerprint,
	    static_cast<uint64_t>(sb.st_mtime));
	fingerprint = bobol_draw_cache_mix(fingerprint,
	    static_cast<uint64_t>(sb.st_ctime));
    }
    return fingerprint;
}

static int bobol_draw_proxy_bbox_valid(const point_t bmin,
	const point_t bmax);

static int
bobol_draw_lod_asset_oriented_bounds_valid(
	const BObolDrawLodAssetRecord *record)
{
    if (!record || record->assetOrientedBoundsValid == 0)
	return record ? 1 : 0;
    if (record->assetOrientedBoundsValid != 1)
	return 0;

    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    VMOVE(hierarchy.quantization_min, record->assetBoundsMin);
    VMOVE(hierarchy.quantization_max, record->assetBoundsMax);
    hierarchy.oriented_bounds_valid = 1;
    for (size_t corner = 0; corner < 8; ++corner)
	VMOVE(hierarchy.oriented_bounds[corner],
	    record->assetOrientedBounds[corner]);
    return bobol_mesh_lod_oriented_bounds_validate(&hierarchy);
}

static size_t
bobol_draw_proxy_disk_size(size_t name_length, size_t point_count)
{
    if (name_length > SIZE_MAX - sizeof(BObolDrawProxyDiskHeader) ||
	point_count > (SIZE_MAX - sizeof(BObolDrawProxyDiskHeader) -
	name_length) / sizeof(point_t))
	return 0;
    return sizeof(BObolDrawProxyDiskHeader) + name_length +
	point_count * sizeof(point_t);
}

static int
bobol_draw_proxy_disk_pack(std::vector<unsigned char> &buffer,
	const db_i *dbip, const directory *dp, const char *name, int kind,
	const point_t *points, size_t point_count)
{
    if (!dbip || !dp || !name || !points || !point_count)
	return 0;
    const size_t name_length = strlen(name);
    const size_t disk_size = bobol_draw_proxy_disk_size(name_length,
	point_count);
    if (!disk_size || name_length > UINT32_MAX || point_count > UINT32_MAX)
	return 0;

    buffer.assign(disk_size, 0);
    BObolDrawProxyDiskHeader header;
    memset(&header, 0, sizeof(header));
    header.magic = BOBOL_DRAW_PROXY_DISK_MAGIC;
    header.version = BOBOL_DRAW_PROXY_DISK_VERSION;
    header.kind = static_cast<uint32_t>(kind);
    header.pointCount = static_cast<uint32_t>(point_count);
    header.nameLength = static_cast<uint32_t>(name_length);
    header.directoryFlags = static_cast<uint32_t>(dp->d_flags);
    header.directoryMajorType = static_cast<uint32_t>(dp->d_major_type);
    header.directoryMinorType = static_cast<uint32_t>(dp->d_minor_type);
    header.databaseFingerprint = bobol_draw_cache_database_fingerprint(dbip);
    header.directoryAddress = static_cast<uint64_t>(dp->d_addr);
    header.directoryLength = static_cast<uint64_t>(dp->d_len);
    memcpy(buffer.data(), &header, sizeof(header));
    memcpy(buffer.data() + sizeof(header), name, name_length);
    memcpy(buffer.data() + sizeof(header) + name_length, points,
	point_count * sizeof(point_t));
    return 1;
}

static int
bobol_draw_proxy_disk_unpack(const void *data, size_t data_size,
	const db_i *dbip, const directory *dp, const char *name, int kind,
	point_t *points, size_t point_count)
{
    if (!data || !dbip || !dp || !name || !points || !point_count ||
	data_size < sizeof(BObolDrawProxyDiskHeader))
	return 0;

    BObolDrawProxyDiskHeader header;
    memcpy(&header, data, sizeof(header));
    const size_t name_length = strlen(name);
    const size_t expected_size = bobol_draw_proxy_disk_size(name_length,
	point_count);
    if (!expected_size || data_size != expected_size ||
	header.magic != BOBOL_DRAW_PROXY_DISK_MAGIC ||
	header.version != BOBOL_DRAW_PROXY_DISK_VERSION ||
	header.kind != static_cast<uint32_t>(kind) ||
	header.pointCount != point_count ||
	header.nameLength != name_length ||
	header.directoryFlags != static_cast<uint32_t>(dp->d_flags) ||
	header.directoryMajorType != static_cast<uint32_t>(dp->d_major_type) ||
	header.directoryMinorType != static_cast<uint32_t>(dp->d_minor_type) ||
	header.databaseFingerprint != bobol_draw_cache_database_fingerprint(dbip) ||
	header.directoryAddress != static_cast<uint64_t>(dp->d_addr) ||
	header.directoryLength != static_cast<uint64_t>(dp->d_len))
	return 0;

    const unsigned char *bytes = static_cast<const unsigned char *>(data);
    if (memcmp(bytes + sizeof(header), name, name_length) != 0)
	return 0;
    memcpy(points, bytes + sizeof(header) + name_length,
	point_count * sizeof(point_t));
    return 1;
}

static size_t
bobol_draw_lod_asset_disk_size(size_t nameLength, size_t assetNameLength)
{
    if (nameLength > SIZE_MAX - sizeof(BObolDrawLodAssetDiskHeader) ||
	assetNameLength > SIZE_MAX - sizeof(BObolDrawLodAssetDiskHeader) -
	    nameLength)
	return 0;
    return sizeof(BObolDrawLodAssetDiskHeader) + nameLength +
	assetNameLength;
}

static int
bobol_draw_lod_asset_disk_pack(std::vector<unsigned char> &buffer,
	const db_i *dbip, db_i *validationDb, const char *name,
	const BObolDrawLodAssetRecord *record)
{
    if (!dbip || !validationDb || !name || !record || !record->faceCount ||
	!record->pointCount ||
	!bobol_draw_proxy_bbox_valid(record->boundsMin, record->boundsMax) ||
	!bobol_draw_proxy_bbox_valid(record->assetBoundsMin,
	    record->assetBoundsMax) ||
	!bobol_draw_lod_asset_oriented_bounds_valid(record))
	return 0;

    const size_t nameLength = strlen(name);
    const size_t assetNameLength = strnlen(record->assetName,
	BOBOL_DRAW_CACHE_LOD_ASSET_NAME_MAX);
    if (!nameLength || !assetNameLength ||
	assetNameLength >= BOBOL_DRAW_CACHE_LOD_ASSET_NAME_MAX)
	return 0;
    /* A canonical asset never needs a persisted reference to itself.  More
     * importantly, a non-identity self mapping is impossible: applying it
     * moves the authored arrays away from their own database coordinates.
     * Reject this at the cache boundary so a future reuse-planner regression
     * cannot poison every subsequent warm draw. */
    if (nameLength == assetNameLength &&
	memcmp(name, record->assetName, nameLength) == 0) {
	for (size_t i = 0; i < 16; ++i) {
	    const fastf_t expected = (i == 0 || i == 5 || i == 10 ||
		i == 15) ? 1.0 : 0.0;
	    if (fabs(record->assetToObject[i] - expected) > VUNITIZE_TOL)
		return 0;
	}
    }
    directory *dp = db_lookup(validationDb, name, LOOKUP_QUIET);
    directory *assetDp = db_lookup(validationDb,
	record->assetName, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL || assetDp == RT_DIR_NULL ||
	dp->d_addr == RT_DIR_PHONY_ADDR ||
	assetDp->d_addr == RT_DIR_PHONY_ADDR ||
	dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT ||
	assetDp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
	return 0;
    for (size_t i = 0; i < 16; i++)
	if (!std::isfinite(record->assetToObject[i]))
	    return 0;

    const size_t diskSize =
	bobol_draw_lod_asset_disk_size(nameLength, assetNameLength);
    if (!diskSize || nameLength > UINT32_MAX ||
	assetNameLength > UINT32_MAX)
	return 0;

    buffer.assign(diskSize, 0);
    BObolDrawLodAssetDiskHeader header;
    memset(&header, 0, sizeof(header));
    header.magic = BOBOL_DRAW_LOD_ASSET_DISK_MAGIC;
    header.version = BOBOL_DRAW_LOD_ASSET_DISK_VERSION;
    header.nameLength = static_cast<uint32_t>(nameLength);
    header.assetNameLength = static_cast<uint32_t>(assetNameLength);
    header.directoryFlags = static_cast<uint32_t>(dp->d_flags);
    header.directoryMajorType = static_cast<uint32_t>(dp->d_major_type);
    header.directoryMinorType = static_cast<uint32_t>(dp->d_minor_type);
    header.assetDirectoryFlags = static_cast<uint32_t>(assetDp->d_flags);
    header.assetDirectoryMajorType =
	static_cast<uint32_t>(assetDp->d_major_type);
    header.assetDirectoryMinorType =
	static_cast<uint32_t>(assetDp->d_minor_type);
    header.databaseFingerprint = bobol_draw_cache_database_fingerprint(dbip);
    header.directoryAddress = static_cast<uint64_t>(dp->d_addr);
    header.directoryLength = static_cast<uint64_t>(dp->d_len);
    header.assetDirectoryAddress = static_cast<uint64_t>(assetDp->d_addr);
    header.assetDirectoryLength = static_cast<uint64_t>(assetDp->d_len);
    header.faceCount = record->faceCount;
    header.pointCount = record->pointCount;
    VMOVE(header.boundsMin, record->boundsMin);
    VMOVE(header.boundsMax, record->boundsMax);
    VMOVE(header.assetBoundsMin, record->assetBoundsMin);
    VMOVE(header.assetBoundsMax, record->assetBoundsMax);
    MAT_COPY(header.assetToObject, record->assetToObject);
    header.assetOrientedBoundsValid =
	static_cast<uint32_t>(record->assetOrientedBoundsValid);
    if (record->assetOrientedBoundsValid)
	for (size_t corner = 0; corner < 8; ++corner)
	    VMOVE(header.assetOrientedBounds[corner],
		record->assetOrientedBounds[corner]);
    memcpy(buffer.data(), &header, sizeof(header));
    memcpy(buffer.data() + sizeof(header), name, nameLength);
    memcpy(buffer.data() + sizeof(header) + nameLength, record->assetName,
	assetNameLength);
    return 1;
}

static int
bobol_draw_lod_asset_disk_unpack(const void *data, size_t dataSize,
	const db_i *dbip, db_i *validationDb, const char *name,
	BObolDrawLodAssetRecord *record)
{
    if (!data || !dbip || !validationDb || !name || !record ||
	dataSize < sizeof(BObolDrawLodAssetDiskHeader))
	return 0;
    directory *dp = db_lookup(validationDb, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
	return 0;

    BObolDrawLodAssetDiskHeader header;
    memcpy(&header, data, sizeof(header));
    const size_t nameLength = strlen(name);
    const size_t expectedSize = bobol_draw_lod_asset_disk_size(
	header.nameLength, header.assetNameLength);
    if (!expectedSize || dataSize != expectedSize ||
	header.magic != BOBOL_DRAW_LOD_ASSET_DISK_MAGIC ||
	header.version != BOBOL_DRAW_LOD_ASSET_DISK_VERSION ||
	header.nameLength != nameLength || !header.assetNameLength ||
	header.assetNameLength >= BOBOL_DRAW_CACHE_LOD_ASSET_NAME_MAX ||
	header.directoryFlags != static_cast<uint32_t>(dp->d_flags) ||
	header.directoryMajorType != static_cast<uint32_t>(dp->d_major_type) ||
	header.directoryMinorType != static_cast<uint32_t>(dp->d_minor_type) ||
	header.databaseFingerprint != bobol_draw_cache_database_fingerprint(dbip) ||
	header.directoryAddress != static_cast<uint64_t>(dp->d_addr) ||
	header.directoryLength != static_cast<uint64_t>(dp->d_len) ||
	!header.faceCount || !header.pointCount ||
	header.assetOrientedBoundsValid > 1 ||
	!bobol_draw_proxy_bbox_valid(header.boundsMin, header.boundsMax) ||
	!bobol_draw_proxy_bbox_valid(header.assetBoundsMin,
	    header.assetBoundsMax))
	return 0;
    for (size_t i = 0; i < 16; i++)
	if (!std::isfinite(header.assetToObject[i]))
	    return 0;

    const unsigned char *bytes = static_cast<const unsigned char *>(data);
    if (memcmp(bytes + sizeof(header), name, nameLength) != 0)
	return 0;
    char assetName[BOBOL_DRAW_CACHE_LOD_ASSET_NAME_MAX] = {0};
    memcpy(assetName, bytes + sizeof(header) + nameLength,
	header.assetNameLength);
    directory *assetDp = db_lookup(validationDb, assetName,
	LOOKUP_QUIET);
    if (assetDp == RT_DIR_NULL ||
	header.assetDirectoryFlags !=
	    static_cast<uint32_t>(assetDp->d_flags) ||
	header.assetDirectoryMajorType !=
	    static_cast<uint32_t>(assetDp->d_major_type) ||
	header.assetDirectoryMinorType !=
	    static_cast<uint32_t>(assetDp->d_minor_type) ||
	header.assetDirectoryAddress !=
	    static_cast<uint64_t>(assetDp->d_addr) ||
	header.assetDirectoryLength != static_cast<uint64_t>(assetDp->d_len))
	return 0;

    bobol_draw_lod_asset_record_init(record);
    bu_strlcpy(record->assetName, assetName, sizeof(record->assetName));
    record->faceCount = header.faceCount;
    record->pointCount = header.pointCount;
    VMOVE(record->boundsMin, header.boundsMin);
    VMOVE(record->boundsMax, header.boundsMax);
    VMOVE(record->assetBoundsMin, header.assetBoundsMin);
    VMOVE(record->assetBoundsMax, header.assetBoundsMax);
    MAT_COPY(record->assetToObject, header.assetToObject);
    record->assetOrientedBoundsValid =
	static_cast<int>(header.assetOrientedBoundsValid);
    if (record->assetOrientedBoundsValid)
	for (size_t corner = 0; corner < 8; ++corner)
	    VMOVE(record->assetOrientedBounds[corner],
		header.assetOrientedBounds[corner]);
    if (!bobol_draw_lod_asset_oriented_bounds_valid(record))
	return 0;
    return 1;
}


static int
bobol_draw_proxy_point_finite(const point_t p)
{
    return std::isfinite(p[X]) && std::isfinite(p[Y]) &&
	   std::isfinite(p[Z]);
}


static int
bobol_draw_proxy_bbox_valid(const point_t bmin,
			      const point_t bmax)
{
    return bobol_draw_proxy_point_finite(bmin) &&
	   bobol_draw_proxy_point_finite(bmax) &&
	   bmin[X] <= bmax[X] && bmin[Y] <= bmax[Y] &&
	   bmin[Z] <= bmax[Z];
}


static void
bobol_draw_cache_key(char *key, const char *name, const char *component)
{
    if (!key)
	return;
    key[0] = '\0';
    if (!name || !component)
	return;

    /* This byte sequence is persistent cache-key identity.  Preserve the
     * historical thirteen-G suffix for short names, but do not allocate and
     * format a bu_vls for every probe.  Large discovery streams call this
     * once per leaf and allocator contention was measurable at 150k parts. */
    const size_t nameLength = strlen(name);
    const char shortSuffix[] = "GGGGGGGGGGGGG";
    char shortName[10 + sizeof(shortSuffix)] = {0};
    const void *hashBytes = name;
    size_t hashLength = nameLength;
    if (nameLength < 10) {
	memcpy(shortName, name, nameLength);
	memcpy(shortName + nameLength, shortSuffix,
	    sizeof(shortSuffix) - 1);
	hashBytes = shortName;
	hashLength += sizeof(shortSuffix) - 1;
    }
    const unsigned long long hash = bu_data_hash(hashBytes, hashLength);
    snprintf(key, BU_CACHE_KEY_MAXLEN, "%llu:%s", hash, component);
}

static void
bobol_draw_proxy_cache_key(char *key, const char *name, int kind)
{
    bobol_draw_cache_key(key, name, bobol_draw_proxy_component(kind));
}

static void
bobol_draw_path_metadata_index_key(char *key, const char *path)
{
    bobol_draw_cache_key(key, path, BOBOL_DRAW_CACHE_PATH_METADATA_INDEX);
}

static int
bobol_draw_cache_key_component_is(const char *key, const char *component)
{
    if (!key || !component)
	return 0;

    const char *sep = strrchr(key, ':');
    return (sep && BU_STR_EQUAL(sep + 1, component)) ? 1 : 0;
}

/* Called with BOBOL_DRAW_CACHE held.  A process-local negative index avoids
 * opening one LMDB read transaction for every cold leaf.  A concurrent writer
 * in another process can only make this process conservatively miss a newly
 * added hint; it can never make an invalid record appear valid. */
static void
bobol_draw_cache_load_lod_asset_keys(BObolDrawCacheContext *context)
{
    if (!context || !context->cache || context->lodAssetKeysLoaded)
	return;

    char **keys = NULL;
    const int nkeys = bu_cache_keys(&keys, context->cache);
    for (int i = 0; i < nkeys; ++i) {
	if (keys[i] && bobol_draw_cache_key_component_is(keys[i],
		BOBOL_DRAW_CACHE_LOD_ASSET))
	    context->lodAssetKeys->insert(keys[i]);
    }
    if (nkeys > 0)
	bu_argv_free(static_cast<size_t>(nkeys), keys);
    context->lodAssetKeysLoaded = true;
}

static int
bobol_draw_path_component_is_instance_of(const char *component,
	size_t componentLen,
	const char *name,
	size_t nameLen)
{
    if (!component || !componentLen || !name || !nameLen)
	return 0;
    if (componentLen == nameLen &&
	bu_strncmp(component, name, nameLen) == 0)
	return 1;
    if (componentLen <= nameLen + 1 ||
	bu_strncmp(component, name, nameLen) != 0 ||
	component[nameLen] != '@')
	return 0;

    for (size_t i = nameLen + 1; i < componentLen; i++) {
	if (component[i] < '0' || component[i] > '9')
	    return 0;
    }
    return 1;
}

static int
bobol_draw_path_contains_object_name(const char *path, const char *name)
{
    if (!path || !path[0] || !name || !name[0])
	return 0;

    size_t nameLen = strlen(name);
    const char *start = path;
    while (*start) {
	const char *end = strchr(start, '/');
	size_t componentLen = end ? (size_t)(end - start) : strlen(start);
	if (bobol_draw_path_component_is_instance_of(start, componentLen,
		name, nameLen))
	    return 1;
	if (!end)
	    break;
	start = end + 1;
    }

    return 0;
}

static int
bobol_draw_cache_ensure_root(void)
{
    char dir[MAXPATHLEN];

    bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, NULL);
    if (dir[0] && !bu_file_exists(dir, NULL))
	bu_mkdir(dir);
    if (dir[0] && !bu_file_exists(dir, NULL))
	return 0;

    bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, BOBOL_DRAW_CACHE_DIR, NULL);
    if (!bu_file_exists(dir, NULL))
	bu_mkdir(dir);
    if (!bu_file_exists(dir, NULL))
	return 0;

    char formatPath[MAXPATHLEN];
    bu_dir(formatPath, MAXPATHLEN, BU_DIR_CACHE, BOBOL_DRAW_CACHE_DIR,
	   BOBOL_DRAW_CACHE_FORMAT_FILE, NULL);
    long diskFormatVersion = -1;
    FILE *fp = fopen(formatPath, "r");
    if (fp) {
	if (fscanf(fp, "%ld", &diskFormatVersion) != 1)
	    diskFormatVersion = -1;
	fclose(fp);
    }
    if (diskFormatVersion > 0 &&
	diskFormatVersion != BOBOL_DRAW_CACHE_CURRENT_FORMAT) {
	bobol_draw_cache_registry_close_all();
	bu_cache_erase(BOBOL_DRAW_CACHE_DIR);
	bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, BOBOL_DRAW_CACHE_DIR, NULL);
	if (!bu_file_exists(dir, NULL))
	    bu_mkdir(dir);
	if (!bu_file_exists(dir, NULL))
	    return 0;
    }
    fp = fopen(formatPath, "w");
    if (fp) {
	fprintf(fp, "%d\n", BOBOL_DRAW_CACHE_CURRENT_FORMAT);
	fclose(fp);
    }

    return 1;
}

static int
bobol_draw_cache_db_name(bu_vls *fname, db_i *dbip)
{
    if (!fname || !dbip)
	return 0;

    const char *ctxName = dbip->dbi_filename;
    char inmemCtxName[128] = {0};
    char canonicalCtxName[MAXPATHLEN] = {0};
    if (!ctxName || !ctxName[0]) {
	snprintf(inmemCtxName, sizeof(inmemCtxName),
		 "bobol_inmem_draw_cache_%p", (void *)dbip);
	ctxName = inmemCtxName;
    } else {
	(void)bu_file_realpath(ctxName, canonicalCtxName);
	ctxName = canonicalCtxName;
    }

    bu_vls_sprintf(fname, "%s", bu_path_normalize(ctxName));
    if (bu_vls_strlen(fname) < 10)
	bu_vls_printf(fname, "GGGGGGGGGGGGG");
    unsigned long long hash = bu_data_hash(bu_vls_cstr(fname),
					   bu_vls_strlen(fname) * sizeof(char));
    bu_path_component(fname, bu_path_normalize(ctxName),
		      BU_PATH_BASENAME_EXTLESS);
    bu_vls_printf(fname, "_%llu", hash);
    return 1;
}

static void
bobol_draw_cache_validation_refresh(BObolDrawCacheContext *context,
	db_i *dbip)
{
    if (!context || !dbip)
	return;
    const uint64_t fingerprint = bobol_draw_cache_database_fingerprint(dbip);
    if (context->validationDb &&
	context->validationFingerprint == fingerprint)
	return;
    if (context->validationDb) {
	db_close(context->validationDb);
	context->validationDb = NULL;
    }
    context->validationFingerprint = fingerprint;
    const char *filename = dbip->dbi_filename;
    if (!filename || !filename[0] || !bu_file_exists(filename, NULL))
	return;
    db_i *validationDb = db_open(filename, DB_OPEN_READONLY);
    if (!validationDb)
	return;
    if (db_dirbuild(validationDb) < 0) {
	db_close(validationDb);
	return;
    }
    context->validationDb = validationDb;
}

static int
bobol_draw_cache_open(BObolDrawCacheHandle *handle, db_i *dbip)
{
    bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls cpath = BU_VLS_INIT_ZERO;
    char cacheRoot[MAXPATHLEN] = {0};

    if (!handle || !dbip)
	return 0;
    handle->context = NULL;
    handle->cache = NULL;

    const char *filename = dbip->dbi_filename ? dbip->dbi_filename : "";
    if (filename[0]) {
	std::lock_guard<std::mutex> guard(bobol_draw_cache_registry_mutex());
	auto bound = bobol_draw_cache_bindings().find(dbip);
	if (bound != bobol_draw_cache_bindings().end()) {
	    if (bound->second.filename == filename &&
		bound->second.context && bound->second.context->cache) {
		handle->context = bound->second.context;
		handle->cache = bound->second.context->cache;
		return 1;
	    }
	    bobol_draw_cache_bindings().erase(bound);
	}
    }

    if (!bobol_draw_cache_ensure_root())
	return 0;
    if (!bobol_draw_cache_db_name(&fname, dbip))
	return 0;

    bu_vls_sprintf(&cpath, "%s/%s_draw", BOBOL_DRAW_CACHE_DIR,
		   bu_vls_cstr(&fname));
    bu_dir(cacheRoot, MAXPATHLEN, BU_DIR_CACHE, NULL);
    std::string registryKey(cacheRoot);
    registryKey.append("/");
    registryKey.append(bu_vls_cstr(&cpath));

    {
	std::lock_guard<std::mutex> guard(bobol_draw_cache_registry_mutex());
	std::map<std::string, BObolDrawCacheContext *>::iterator it =
	    bobol_draw_cache_registry().find(registryKey);
	if (it != bobol_draw_cache_registry().end()) {
	    handle->context = it->second;
	    handle->cache = it->second->cache;
	    if (filename[0])
		bobol_draw_cache_bindings()[dbip] = {filename, it->second};
	    bu_vls_free(&cpath);
	    bu_vls_free(&fname);
	    return handle->cache ? 1 : 0;
	}

	BObolDrawCacheContext *context;
	BU_GET(context, BObolDrawCacheContext);
	context->cache = bu_cache_open(bu_vls_cstr(&cpath), 1, 0);
	context->registryKey = bu_strdup(registryKey.c_str());
	context->validationDb = NULL;
	context->validationFingerprint = 0;
	context->lodAssetKeys = new std::unordered_set<std::string>;
	context->lodAssetOrientedBoundsKeys =
	    new std::unordered_set<std::string>;
	context->lodAssetKeysLoaded = false;
	if (!context->cache) {
	    bobol_draw_cache_context_close(context);
	    bu_vls_free(&cpath);
	    bu_vls_free(&fname);
	    return 0;
	}
	bobol_draw_cache_registry()[registryKey] = context;
	if (filename[0])
	    bobol_draw_cache_bindings()[dbip] = {filename, context};
	handle->context = context;
	handle->cache = context->cache;
    }

    bu_vls_free(&cpath);
    bu_vls_free(&fname);
    return handle->cache ? 1 : 0;
}

static void
bobol_draw_cache_close(BObolDrawCacheHandle *handle)
{
    if (!handle)
	return;
    handle->context = NULL;
    handle->cache = NULL;
}

extern "C" void
bobol_draw_cache_status_init(BObolDrawCacheStatus *status)
{
    if (!status)
	return;
    memset(status, 0, sizeof(*status));
}

extern "C" void
bobol_draw_proxy_record_init(BObolDrawProxyRecord *record)
{
    if (!record)
	return;
    memset(record, 0, sizeof(*record));
}

extern "C" void
bobol_draw_lod_asset_record_init(BObolDrawLodAssetRecord *record)
{
    if (!record)
	return;
    memset(record, 0, sizeof(*record));
    MAT_IDN(record->assetToObject);
}

extern "C" void
bobol_draw_metadata_record_init(BObolDrawMetadataRecord *record)
{
    if (!record)
	return;
    memset(record, 0, sizeof(*record));
}

extern "C" void
bobol_draw_cache_clear_all(void)
{
    char dir[MAXPATHLEN];
    int sem = bobol_draw_cache_semaphore();

    bu_semaphore_acquire(sem);
    bobol_draw_cache_registry_close_all();
    bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, BOBOL_DRAW_CACHE_DIR, NULL);
    bu_dirclear(dir);
    bu_semaphore_release(sem);
}

extern "C" int
bobol_draw_cache_clear_database(db_i *dbip)
{
    BObolDrawCacheHandle handle;

    if (!dbip)
	return BRLCAD_ERROR;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    int opened = bobol_draw_cache_open(&handle, dbip);
    if (opened) {
	opened = bu_cache_clear_all(handle.cache) == BRLCAD_OK ? 1 : 0;
	if (opened && handle.context) {
	    handle.context->lodAssetKeys->clear();
	    handle.context->lodAssetOrientedBoundsKeys->clear();
	    handle.context->lodAssetKeysLoaded = true;
	}
    }
    bobol_draw_cache_close(&handle);
    bu_semaphore_release(sem);

    return opened ? BRLCAD_OK : BRLCAD_ERROR;
}

extern "C" int
bobol_draw_manifest_cache_invalidate_database(db_i *dbip)
{
    BObolDrawCacheHandle handle;
    char **keys = NULL;
    int nkeys = 0;

    if (!dbip)
	return BRLCAD_ERROR;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    int opened = bobol_draw_cache_open(&handle, dbip);
    if (opened) {
	nkeys = bu_cache_keys(&keys, handle.cache);
	struct bu_cache_txn *writeTxn = NULL;
	for (int i = 0; i < nkeys; i++) {
	    if (bobol_draw_cache_key_component_is(keys[i],
		    BOBOL_DRAW_CACHE_MANIFEST) ||
		bobol_draw_cache_key_component_is(keys[i],
		    BOBOL_DRAW_CACHE_LOD_ASSET))
		bu_cache_clear(keys[i], handle.cache, &writeTxn);
	}
	if (writeTxn &&
	    bu_cache_write_commit(handle.cache, &writeTxn) != BRLCAD_OK)
	    opened = 0;
	if (opened && handle.context) {
	    handle.context->lodAssetKeys->clear();
	    handle.context->lodAssetOrientedBoundsKeys->clear();
	    handle.context->lodAssetKeysLoaded = true;
	}
    }
    bobol_draw_cache_close(&handle);
    bu_semaphore_release(sem);

    if (nkeys)
	bu_argv_free((size_t)nkeys, keys);
    return opened ? BRLCAD_OK : BRLCAD_ERROR;
}

static void
bobol_draw_proxy_status_current(db_i *dbip,
				  const char *name,
				  int kind,
				  BObolDrawCacheStatus *status)
{
    void *data = NULL;
    char key[BU_CACHE_KEY_MAXLEN] = {0};
    const size_t pointCount = bobol_draw_proxy_point_count(kind);

    if (!status)
	return;
    bobol_draw_cache_status_init(status);
    if (!dbip || !name || !pointCount)
	return;

    directory *dp = RT_DIR_NULL;
    int cacheable = bobol_draw_proxy_cache_directory(dbip, name, &dp);
    status->directoryFound = (dp != RT_DIR_NULL) ? 1 : 0;
    if (!cacheable)
	return;

    bobol_draw_proxy_cache_key(key, name, kind);
    if (!key[0])
	return;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
	if (bobol_draw_cache_open(&handle, dbip)) {
	size_t dsize = bu_cache_get(&data, key, handle.cache, NULL);
	point_t points[8];
	status->hasCachedPayload = bobol_draw_proxy_disk_unpack(data, dsize,
	    dbip, dp, name, kind, points, pointCount);
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    if (data)
	bu_free(data, "bobol draw proxy status data");
}

extern "C" int
bobol_draw_proxy_cache_status(db_i *dbip,
				const char *name,
				int kind,
				BObolDrawCacheStatus *status)
{
    if (!status || !dbip || !name || !bobol_draw_proxy_point_count(kind))
	return BRLCAD_ERROR;
    bobol_draw_proxy_status_current(dbip, name, kind, status);
    return BRLCAD_OK;
}

extern "C" int
bobol_draw_proxy_cache_invalidate(db_i *dbip,
				    const char *name,
				    int kind,
				    BObolDrawCacheStatus *status)
{
    BObolDrawCacheStatus current;
    char key[BU_CACHE_KEY_MAXLEN] = {0};

    if (status)
	bobol_draw_cache_status_init(status);
    if (!dbip || !name || !bobol_draw_proxy_point_count(kind))
	return BRLCAD_ERROR;

    bobol_draw_proxy_status_current(dbip, name, kind, &current);
    bobol_draw_proxy_cache_key(key, name, kind);
    if (!key[0])
	return BRLCAD_ERROR;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    if (bobol_draw_cache_open(&handle, dbip)) {
	bu_cache_clear(key, handle.cache, NULL);
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    current.clearedCacheEntry = 1;
    current.hasCachedPayload = 0;
    if (status)
	*status = current;
    return BRLCAD_OK;
}

extern "C" int
bobol_draw_proxy_cache_store(db_i *dbip,
			       const char *name,
			       int kind,
			       const point_t *points,
			       size_t pointCount,
			       BObolDrawCacheStatus *status)
{
    BObolDrawCacheStatus current;
    char key[BU_CACHE_KEY_MAXLEN] = {0};
    size_t expectedCount = bobol_draw_proxy_point_count(kind);
    std::vector<unsigned char> disk;

    if (status)
	bobol_draw_cache_status_init(status);
    if (!dbip || !name || !points || !expectedCount ||
	pointCount != expectedCount)
	return BRLCAD_ERROR;
    bobol_draw_cache_status_init(&current);

    directory *dp = RT_DIR_NULL;
    int cacheable = bobol_draw_proxy_cache_directory(dbip, name, &dp);
    current.directoryFound = (dp != RT_DIR_NULL) ? 1 : 0;
    if (!cacheable || !bobol_draw_proxy_disk_pack(disk, dbip, dp, name,
	kind, points, pointCount)) {
	if (status)
	    *status = current;
	return BRLCAD_ERROR;
    }

    bobol_draw_proxy_cache_key(key, name, kind);
    if (!key[0])
	return BRLCAD_ERROR;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    size_t wsize = 0;
    if (bobol_draw_cache_open(&handle, dbip)) {
	wsize = bu_cache_write(disk.data(), disk.size(), key, handle.cache,
	    NULL);
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    current.hasCachedPayload = (wsize == disk.size()) ? 1 : 0;
    current.generatedCacheEntry = current.hasCachedPayload ? 1 : 0;
    if (status)
	*status = current;
    return current.hasCachedPayload ? BRLCAD_OK : BRLCAD_ERROR;
}

extern "C" int
bobol_draw_proxy_cache_get(db_i *dbip,
			     const char *name,
			     int kind,
			     BObolDrawProxyRecord *record)
{
    void *data = NULL;
    char key[BU_CACHE_KEY_MAXLEN] = {0};
    size_t expectedCount = bobol_draw_proxy_point_count(kind);

    if (!dbip || !name || !record || !expectedCount)
	return BRLCAD_ERROR;
    bobol_draw_proxy_record_init(record);
    directory *dp = RT_DIR_NULL;
    if (!bobol_draw_proxy_cache_directory(dbip, name, &dp))
	return BRLCAD_ERROR;

    bobol_draw_proxy_cache_key(key, name, kind);
    if (!key[0])
	return BRLCAD_ERROR;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    size_t dsize = 0;
    if (bobol_draw_cache_open(&handle, dbip)) {
	dsize = bu_cache_get(&data, key, handle.cache, NULL);
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    if (!bobol_draw_proxy_disk_unpack(data, dsize, dbip, dp, name, kind,
	record->points, expectedCount)) {
	if (data)
	    bu_free(data, "bobol draw proxy get data");
	return BRLCAD_ERROR;
    }

    record->kind = kind;
    record->pointCount = expectedCount;
    bu_free(data, "bobol draw proxy get data");
    return BRLCAD_OK;
}

extern "C" int
bobol_draw_lod_asset_cache_store(db_i *dbip, const char *name,
	const BObolDrawLodAssetRecord *record)
{
    if (!dbip || !name || !record)
	return BRLCAD_ERROR;

    char key[BU_CACHE_KEY_MAXLEN] = {0};
    bobol_draw_cache_key(key, name, BOBOL_DRAW_CACHE_LOD_ASSET);
    if (!key[0])
	return BRLCAD_ERROR;

    std::vector<unsigned char> disk;
    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    size_t written = 0;
    if (bobol_draw_cache_open(&handle, dbip)) {
	bobol_draw_cache_validation_refresh(handle.context, dbip);
	db_i *validationDb = handle.context &&
	    handle.context->validationDb ?
	    handle.context->validationDb : dbip;
	BObolDrawLodAssetRecord merged = *record;
	/* A cold reuse proof may finish after the resident PoP worker has
	 * published its OBB.  Preserve that optional field when the canonical
	 * identity is unchanged.  The process-local positive hint prevents this
	 * merge from adding a read transaction to ordinary 50k-scale stores. */
	if (!merged.assetOrientedBoundsValid && handle.context &&
	    handle.context->lodAssetOrientedBoundsKeys->find(key) !=
		handle.context->lodAssetOrientedBoundsKeys->end()) {
	    void *priorData = NULL;
	    const size_t priorSize = bu_cache_get(
		&priorData, key, handle.cache, NULL);
	    BObolDrawLodAssetRecord prior;
	    if (priorData && bobol_draw_lod_asset_disk_unpack(
		    priorData, priorSize, dbip, validationDb, name, &prior) &&
		BU_STR_EQUAL(prior.assetName, merged.assetName) &&
		prior.assetOrientedBoundsValid == 1) {
		merged.assetOrientedBoundsValid = 1;
		for (size_t corner = 0; corner < 8; ++corner)
		    VMOVE(merged.assetOrientedBounds[corner],
			prior.assetOrientedBounds[corner]);
	    }
	    if (priorData)
		bu_free(priorData, "prior draw LoD asset data");
	}
	if (bobol_draw_lod_asset_disk_pack(disk, dbip, validationDb, name,
		&merged))
	    written = bu_cache_write(disk.data(), disk.size(), key,
		handle.cache, NULL);
	if (!disk.empty() && written == disk.size()) {
	    bobol_draw_cache_load_lod_asset_keys(handle.context);
	    handle.context->lodAssetKeys->insert(key);
	    if (merged.assetOrientedBoundsValid)
		handle.context->lodAssetOrientedBoundsKeys->insert(key);
	    else
		handle.context->lodAssetOrientedBoundsKeys->erase(key);
	}
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);
    return !disk.empty() && written == disk.size() ?
	BRLCAD_OK : BRLCAD_ERROR;
}

extern "C" int
bobol_draw_lod_asset_cache_get(db_i *dbip, const char *name,
	BObolDrawLodAssetRecord *record)
{
    if (!dbip || !name || !record)
	return BRLCAD_ERROR;
    bobol_draw_lod_asset_record_init(record);

    char key[BU_CACHE_KEY_MAXLEN] = {0};
    bobol_draw_cache_key(key, name, BOBOL_DRAW_CACHE_LOD_ASSET);
    if (!key[0])
	return BRLCAD_ERROR;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    void *data = NULL;
    size_t dataSize = 0;
    int valid = 0;
    if (bobol_draw_cache_open(&handle, dbip)) {
	bobol_draw_cache_load_lod_asset_keys(handle.context);
	const bool keyPresent = handle.context->lodAssetKeys->find(key) !=
	    handle.context->lodAssetKeys->end();
	if (keyPresent)
	    dataSize = bu_cache_get(&data, key, handle.cache, NULL);
	/*
	 * Validation protects a persistent hit from a concurrently edited live
	 * database by comparing it with a read-only directory snapshot.  A miss
	 * has no payload to trust and therefore needs no snapshot at all.  Opening
	 * and directory-building a multi-gigabyte database before testing the
	 * cache key made the first of 150k cold probes serialize every coverage
	 * worker behind several seconds of work whose answer was necessarily
	 * "not cached".
	 */
	if (data && dataSize >= sizeof(BObolDrawLodAssetDiskHeader)) {
	    bobol_draw_cache_validation_refresh(handle.context, dbip);
	    db_i *validationDb = handle.context &&
		handle.context->validationDb ?
		handle.context->validationDb : dbip;
	    valid = bobol_draw_lod_asset_disk_unpack(data, dataSize,
		dbip, validationDb, name, record);
	    if (valid && record->assetOrientedBoundsValid == 1)
		handle.context->lodAssetOrientedBoundsKeys->insert(key);
	}
	/* A key deleted by another process after our snapshot is safe to retire
	 * locally.  An invalid nonempty payload remains indexed: database edits
	 * may be followed by undo, and validation must continue rejecting or
	 * accepting it from authoritative disk identity rather than a miss hint. */
	if (keyPresent && !dataSize)
	    handle.context->lodAssetKeys->erase(key);
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    if (data)
	bu_free(data, "bobol draw LoD asset data");
    return valid ? BRLCAD_OK : BRLCAD_ERROR;
}

int
bobol_draw_lod_asset_oriented_bounds_publish(
    struct db_i *dbip, const char *name, uint64_t faceCount,
    uint64_t pointCount, const point_t boundsMin, const point_t boundsMax,
    const point_t orientedBounds[8])
{
    if (!dbip || !name || !name[0] || !faceCount || !pointCount ||
	!orientedBounds || !bobol_draw_proxy_bbox_valid(boundsMin, boundsMax))
	return BRLCAD_ERROR;

    BObolDrawLodAssetRecord record;
    if (bobol_draw_lod_asset_cache_get(dbip, name, &record) != BRLCAD_OK) {
	bobol_draw_lod_asset_record_init(&record);
	bu_strlcpy(record.assetName, name, sizeof(record.assetName));
	record.faceCount = faceCount;
	record.pointCount = pointCount;
	VMOVE(record.boundsMin, boundsMin);
	VMOVE(record.boundsMax, boundsMax);
	VMOVE(record.assetBoundsMin, boundsMin);
	VMOVE(record.assetBoundsMax, boundsMax);
    } else if (!BU_STR_EQUAL(record.assetName, name)) {
	/* orientedBounds describe name itself, not an alias target. */
	return BRLCAD_ERROR;
    } else if (record.assetOrientedBoundsValid == 1) {
	return BRLCAD_OK;
    }

    record.assetOrientedBoundsValid = 1;
    for (size_t corner = 0; corner < 8; ++corner)
	VMOVE(record.assetOrientedBounds[corner], orientedBounds[corner]);
    return bobol_draw_lod_asset_cache_store(dbip, name, &record);
}

extern "C" int
bobol_draw_proxy_cache_refresh(db_i *dbip,
				 const char *name,
				 int kind,
				 BObolDrawCacheStatus *status)
{
    point_t points[8];
    rt_db_internal intern;
    const bn_tol tol = BN_TOL_INIT_TOL;

    if (status)
	bobol_draw_cache_status_init(status);
    if (!dbip || !name || !bobol_draw_proxy_point_count(kind))
	return BRLCAD_ERROR;

    std::lock_guard<std::mutex> refreshGuard(
	bobol_draw_proxy_refresh_mutex(dbip, name, kind));

    /* A competing occurrence may have filled this persistent record while
     * this request waited for its leaf/kind refresh lock. */
    BObolDrawCacheStatus cached;
    if (bobol_draw_proxy_cache_status(dbip, name, kind, &cached) ==
	BRLCAD_OK && cached.hasCachedPayload) {
	if (status)
	    *status = cached;
	return BRLCAD_OK;
    }

    directory *dp = RT_DIR_NULL;
    int leaf_solid = bobol_draw_proxy_leaf_solid_name(dbip, name, &dp);
    if (!leaf_solid) {
	if (status)
	    status->directoryFound = (dp != RT_DIR_NULL) ? 1 : 0;
	return BRLCAD_ERROR;
    }

    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
	return BRLCAD_ERROR;

    int ret = BRLCAD_ERROR;
    if (kind == BOBOL_LOD_PROXY_AABB) {
	point_t bmin, bmax;
	VSETALL(bmin, INFINITY);
	VSETALL(bmax, -INFINITY);
	if (intern.idb_meth && intern.idb_meth->ft_bbox &&
	    intern.idb_meth->ft_bbox(&intern, &bmin, &bmax,
				     &tol) == 0 &&
	    bobol_draw_proxy_bbox_valid(bmin, bmax)) {
	    VMOVE(points[0], bmin);
	    VMOVE(points[1], bmax);
	    ret = bobol_draw_proxy_cache_store(dbip, name, kind, points,
						 2, status);
	}
    } else if (kind == BOBOL_LOD_PROXY_OBB) {
	rt_arb_internal arb;
	arb.magic = RT_ARB_INTERNAL_MAGIC;
	if (intern.idb_meth && intern.idb_meth->ft_oriented_bbox &&
	    intern.idb_meth->ft_oriented_bbox(&arb, &intern,
					      BN_TOL_DIST) == 0) {
	    for (size_t i = 0; i < 8; i++)
		VMOVE(points[i], arb.pt[i]);
	    ret = bobol_draw_proxy_cache_store(dbip, name, kind, points,
						 8, status);
	}
    }

    rt_db_free_internal(&intern);
    if (ret != BRLCAD_OK)
	(void)bobol_draw_proxy_cache_invalidate(dbip, name, kind, status);
    return ret;
}

static int
bobol_draw_metadata_int_attr(const bu_attribute_value_set *avs,
			       const char *name,
			       int *out)
{
    const char *val;
    int parsed = 0;

    if (!avs || !name || !out)
	return 0;
    val = bu_avs_get(avs, name);
    if (!val || !val[0])
	return 0;
    if (bu_opt_int(NULL, 1, &val, (void *)&parsed) != 1)
	return 0;
    *out = parsed;
    return 1;
}

static int
bobol_draw_metadata_color_attr(const bu_attribute_value_set *avs,
				 BObolDrawMetadataRecord *record)
{
    bu_color color = BU_COLOR_INIT_ZERO;
    const char *val;

    if (!avs || !record)
	return 0;
    val = bu_avs_get(avs, "color");
    if (!val || !val[0])
	val = bu_avs_get(avs, "rgb");
    if (!val || !val[0])
	return 0;
    if (bu_opt_color(NULL, 1, &val, (void *)&color) != 1)
	return 0;
    bu_color_to_rgb_chars(&color, record->color);
    record->hasColor = 1;
    return 1;
}

static void
bobol_draw_metadata_string_attr(const bu_attribute_value_set *avs,
				  const char *name,
				  char *out,
				  size_t outSize,
				  int *hasOut)
{
    const char *val;

    if (hasOut)
	*hasOut = 0;
    if (!avs || !name || !out || !outSize)
	return;
    out[0] = '\0';
    val = bu_avs_get(avs, name);
    if (!val || !val[0])
	return;
    bu_strlcpy(out, val, outSize);
    if (hasOut)
	*hasOut = 1;
}

static void
bobol_draw_metadata_to_disk(BObolDrawMetadataDiskRecord *disk,
			      const BObolDrawMetadataRecord *record)
{
    memset(disk, 0, sizeof(*disk));
    disk->magic = BOBOL_DRAW_METADATA_DISK_MAGIC;
    disk->version = BOBOL_DRAW_CACHE_CURRENT_FORMAT;
    disk->directoryFound = record->directoryFound;
    disk->isPhony = record->isPhony;
    disk->flags = record->flags;
    disk->majorType = record->majorType;
    disk->minorType = record->minorType;
    disk->isSolid = record->isSolid;
    disk->isComb = record->isComb;
    disk->isRegion = record->isRegion;
    disk->isHidden = record->isHidden;
    disk->hasRegionId = record->hasRegionId;
    disk->regionId = record->regionId;
    disk->hasAircode = record->hasAircode;
    disk->aircode = record->aircode;
    disk->hasLos = record->hasLos;
    disk->los = record->los;
    disk->hasMaterialId = record->hasMaterialId;
    disk->materialId = record->materialId;
    disk->hasInherit = record->hasInherit;
    disk->inherit = record->inherit;
    disk->hasColor = record->hasColor;
    memcpy(disk->color, record->color, sizeof(disk->color));
    disk->hasMaterialName = record->hasMaterialName;
    bu_strlcpy(disk->materialName, record->materialName,
	       sizeof(disk->materialName));
    disk->hasShader = record->hasShader;
    bu_strlcpy(disk->shader, record->shader, sizeof(disk->shader));
}

static void
bobol_draw_metadata_from_disk(BObolDrawMetadataRecord *record,
				const BObolDrawMetadataDiskRecord *disk)
{
    bobol_draw_metadata_record_init(record);
    record->directoryFound = disk->directoryFound;
    record->isPhony = disk->isPhony;
    record->flags = disk->flags;
    record->majorType = disk->majorType;
    record->minorType = disk->minorType;
    record->isSolid = disk->isSolid;
    record->isComb = disk->isComb;
    record->isRegion = disk->isRegion;
    record->isHidden = disk->isHidden;
    record->hasRegionId = disk->hasRegionId;
    record->regionId = disk->regionId;
    record->hasAircode = disk->hasAircode;
    record->aircode = disk->aircode;
    record->hasLos = disk->hasLos;
    record->los = disk->los;
    record->hasMaterialId = disk->hasMaterialId;
    record->materialId = disk->materialId;
    record->hasInherit = disk->hasInherit;
    record->inherit = disk->inherit;
    record->hasColor = disk->hasColor;
    memcpy(record->color, disk->color, sizeof(record->color));
    record->hasMaterialName = disk->hasMaterialName;
    bu_strlcpy(record->materialName, disk->materialName,
	       sizeof(record->materialName));
    record->hasShader = disk->hasShader;
    bu_strlcpy(record->shader, disk->shader, sizeof(record->shader));
}

static int
bobol_draw_metadata_disk_valid(const void *data, size_t dsize)
{
    const BObolDrawMetadataDiskRecord *disk =
	(const BObolDrawMetadataDiskRecord *)data;
    return data && dsize == sizeof(*disk) &&
	   disk->magic == BOBOL_DRAW_METADATA_DISK_MAGIC &&
	   disk->version == BOBOL_DRAW_CACHE_CURRENT_FORMAT;
}

static int
bobol_draw_metadata_component_found(db_i *dbip, const char *component)
{
    if (!dbip || !component || !component[0])
	return 0;

    if (db_lookup(dbip, component, LOOKUP_QUIET) != RT_DIR_NULL)
	return 1;

    const char *at = strrchr(component, '@');
    if (!at || at == component || !at[1])
	return 0;
    for (const char *cp = at + 1; *cp; cp++) {
	if (*cp < '0' || *cp > '9')
	    return 0;
    }

    std::string base(component, static_cast<size_t>(at - component));
    return db_lookup(dbip, base.c_str(), LOOKUP_QUIET) != RT_DIR_NULL;
}

static int
bobol_draw_metadata_path_found(db_i *dbip, const char *name)
{
    if (!dbip || !name || !name[0])
	return 0;

    const char *start = name;
    while (*start == '/')
	start++;
    if (!start[0])
	return 1;

    char *copy = bu_strdup(start);
    size_t len = strlen(copy);
    while (len > 0 && copy[len - 1] == '/') {
	copy[len - 1] = '\0';
	len--;
    }

    int found = len > 0 ? 1 : 0;
    char *component = copy;
    while (found && component && component[0]) {
	char *slash = strchr(component, '/');
	if (slash)
	    *slash = '\0';
	if (!bobol_draw_metadata_component_found(dbip, component))
	    found = 0;
	component = slash ? slash + 1 : NULL;
    }

    bu_free(copy, "bobol draw metadata path probe");
    return found;
}

static int
bobol_draw_metadata_target_found(db_i *dbip, const char *name, int pathAware)
{
    if (!dbip || !name || !name[0])
	return 0;

    if (!pathAware)
	return db_lookup(dbip, name, LOOKUP_QUIET) != RT_DIR_NULL;

    return bobol_draw_metadata_path_found(dbip, name);
}

static void
bobol_draw_metadata_record_apply_directory(
    BObolDrawMetadataRecord *record,
    const directory *dp)
{
    if (!record || !dp)
	return;

    record->directoryFound = 1;
    record->isPhony = (dp->d_addr == RT_DIR_PHONY_ADDR) ? 1 : 0;
    record->flags = dp->d_flags;
    record->majorType = dp->d_major_type;
    record->minorType = dp->d_minor_type;
    record->isSolid = (dp->d_flags & RT_DIR_SOLID) ? 1 : 0;
    record->isComb = (dp->d_flags & RT_DIR_COMB) ? 1 : 0;
    record->isRegion = (dp->d_flags & RT_DIR_REGION) ? 1 : 0;
    record->isHidden = (dp->d_flags & RT_DIR_HIDDEN) ? 1 : 0;
}

static int
bobol_draw_metadata_color_channel(float value)
{
    long color = lround(value * 255.0f);
    if (color < 0)
	color = 0;
    if (color > 255)
	color = 255;
    return (int)color;
}

static void
bobol_draw_metadata_record_set_color_from_mater(
    BObolDrawMetadataRecord *record,
    const mater_info *mater)
{
    if (!record || !mater || !mater->ma_color_valid)
	return;

    record->hasColor = 1;
    record->color[0] =
	(unsigned char)bobol_draw_metadata_color_channel(
	    mater->ma_color[0]);
    record->color[1] =
	(unsigned char)bobol_draw_metadata_color_channel(
	    mater->ma_color[1]);
    record->color[2] =
	(unsigned char)bobol_draw_metadata_color_channel(
	    mater->ma_color[2]);
}

static void
bobol_draw_metadata_record_apply_tree_state(
    BObolDrawMetadataRecord *record,
    const db_tree_state *state)
{
    if (!record || !state)
	return;

    if (state->ts_sofar & TS_SOFAR_REGION) {
	record->isRegion = 1;
	record->hasRegionId = 1;
	record->regionId = state->ts_regionid;
	record->hasAircode = 1;
	record->aircode = state->ts_aircode;
	record->hasLos = 1;
	record->los = state->ts_los;
	record->hasMaterialId = 1;
	record->materialId = state->ts_gmater;
    }

    bobol_draw_metadata_record_set_color_from_mater(record,
						      &state->ts_mater);
    if (!record->hasColor && record->hasRegionId) {
	struct region rp;
	memset(&rp, 0, sizeof(rp));
	rp.reg_regionid = record->regionId;
	db_mater_color_region(state->ts_dbip, &rp);
	bobol_draw_metadata_record_set_color_from_mater(record,
							  &rp.reg_mater);
    }
    if (state->ts_mater.ma_shader && state->ts_mater.ma_shader[0]) {
	record->hasShader = 1;
	bu_strlcpy(record->shader, state->ts_mater.ma_shader,
		   sizeof(record->shader));
    }
    if (state->ts_mater.ma_cinherit == DB_INH_HIGHER ||
	state->ts_mater.ma_minherit == DB_INH_HIGHER) {
	record->hasInherit = 1;
	record->inherit = 1;
    }
}

extern "C" void
bobol_draw_metadata_record_from_tree_state(
    BObolDrawMetadataRecord *record,
    const db_tree_state *state,
    const directory *dp)
{
    if (!record)
	return;

    bobol_draw_metadata_record_init(record);
    bobol_draw_metadata_record_apply_directory(record, dp);
    bobol_draw_metadata_record_apply_tree_state(record, state);
}

static void
bobol_draw_metadata_status_current_for_component(
    db_i *dbip,
    const char *name,
    const char *component,
    int pathAware,
    BObolDrawCacheStatus *status)
{
    void *data = NULL;
    char key[BU_CACHE_KEY_MAXLEN] = {0};

    if (!status)
	return;
    bobol_draw_cache_status_init(status);
    if (!dbip || !name)
	return;

    status->directoryFound =
	bobol_draw_metadata_target_found(dbip, name, pathAware);

    bobol_draw_cache_key(key, name, component);
    if (!key[0])
	return;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    if (bobol_draw_cache_open(&handle, dbip)) {
	size_t dsize = bu_cache_get(&data, key, handle.cache, NULL);
	status->hasCachedPayload =
	    bobol_draw_metadata_disk_valid(data, dsize) ? 1 : 0;
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    if (data)
	bu_free(data, "bobol draw metadata status data");
}

static int
bobol_draw_metadata_cache_invalidate_for_component(
    db_i *dbip,
    const char *name,
    const char *component,
    int pathAware,
    BObolDrawCacheStatus *status)
{
    BObolDrawCacheStatus current;
    char key[BU_CACHE_KEY_MAXLEN] = {0};
    char pathIndexKey[BU_CACHE_KEY_MAXLEN] = {0};

    if (status)
	bobol_draw_cache_status_init(status);
    if (!dbip || !name || !component)
	return BRLCAD_ERROR;
    bobol_draw_metadata_status_current_for_component(dbip, name, component,
	    pathAware, &current);
    bobol_draw_cache_key(key, name, component);
    if (!key[0])
	return BRLCAD_ERROR;
    if (BU_STR_EQUAL(component, BOBOL_DRAW_CACHE_PATH_METADATA))
	bobol_draw_path_metadata_index_key(pathIndexKey, name);

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    if (bobol_draw_cache_open(&handle, dbip)) {
	bu_cache_clear(key, handle.cache, NULL);
	if (pathIndexKey[0])
	    bu_cache_clear(pathIndexKey, handle.cache, NULL);
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    current.clearedCacheEntry = 1;
    current.hasCachedPayload = 0;
    if (status)
	*status = current;
    return BRLCAD_OK;
}

static int
bobol_draw_metadata_cache_store_for_component(
    db_i *dbip,
    const char *name,
    const char *component,
    int pathAware,
    const BObolDrawMetadataRecord *record,
    BObolDrawCacheStatus *status)
{
    BObolDrawCacheStatus current;
    BObolDrawMetadataDiskRecord disk;
    char key[BU_CACHE_KEY_MAXLEN] = {0};
    char pathIndexKey[BU_CACHE_KEY_MAXLEN] = {0};

    if (status)
	bobol_draw_cache_status_init(status);
    if (!dbip || !name || !component || !record)
	return BRLCAD_ERROR;
    bobol_draw_cache_status_init(&current);

    current.directoryFound =
	bobol_draw_metadata_target_found(dbip, name, pathAware);
    if (!current.directoryFound) {
	if (status)
	    *status = current;
	return BRLCAD_ERROR;
    }

    bobol_draw_cache_key(key, name, component);
    if (!key[0])
	return BRLCAD_ERROR;
    if (BU_STR_EQUAL(component, BOBOL_DRAW_CACHE_PATH_METADATA))
	bobol_draw_path_metadata_index_key(pathIndexKey, name);
    bobol_draw_metadata_to_disk(&disk, record);

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    size_t wsize = 0;
    size_t iwsize = 0;
    if (bobol_draw_cache_open(&handle, dbip)) {
	wsize = bu_cache_write((void *)&disk, sizeof(disk), key,
			       handle.cache, NULL);
	if (pathIndexKey[0] && wsize == sizeof(disk))
	    iwsize = bu_cache_write((void *)name, strlen(name) + 1,
				    pathIndexKey, handle.cache, NULL);
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    current.hasCachedPayload =
	(wsize == sizeof(disk) &&
	 (!pathIndexKey[0] || iwsize == strlen(name) + 1)) ? 1 : 0;
    current.generatedCacheEntry = current.hasCachedPayload ? 1 : 0;
    if (status)
	*status = current;
    return current.hasCachedPayload ? BRLCAD_OK : BRLCAD_ERROR;
}

static int
bobol_draw_metadata_cache_get_for_component(
    db_i *dbip,
    const char *name,
    const char *component,
    BObolDrawMetadataRecord *record)
{
    void *data = NULL;
    char key[BU_CACHE_KEY_MAXLEN] = {0};

    if (!dbip || !name || !component || !record)
	return BRLCAD_ERROR;
    bobol_draw_metadata_record_init(record);

    bobol_draw_cache_key(key, name, component);
    if (!key[0])
	return BRLCAD_ERROR;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    size_t dsize = 0;
    if (bobol_draw_cache_open(&handle, dbip)) {
	dsize = bu_cache_get(&data, key, handle.cache, NULL);
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    if (!bobol_draw_metadata_disk_valid(data, dsize)) {
	if (data)
	    bu_free(data, "bobol draw metadata get data");
	return BRLCAD_ERROR;
    }

    bobol_draw_metadata_from_disk(record,
				    (const BObolDrawMetadataDiskRecord *)data);
    bu_free(data, "bobol draw metadata get data");
    return BRLCAD_OK;
}

extern "C" int
bobol_draw_metadata_cache_status(db_i *dbip,
				   const char *name,
				   BObolDrawCacheStatus *status)
{
    if (!status || !dbip || !name)
	return BRLCAD_ERROR;
    bobol_draw_metadata_status_current_for_component(dbip, name,
	    BOBOL_DRAW_CACHE_METADATA, 0, status);
    return BRLCAD_OK;
}

extern "C" int
bobol_draw_metadata_cache_invalidate(db_i *dbip,
				       const char *name,
				       BObolDrawCacheStatus *status)
{
    return bobol_draw_metadata_cache_invalidate_for_component(
	       dbip, name, BOBOL_DRAW_CACHE_METADATA, 0, status);
}

extern "C" int
bobol_draw_metadata_cache_store(db_i *dbip,
				  const char *name,
				  const BObolDrawMetadataRecord *record,
				  BObolDrawCacheStatus *status)
{
    return bobol_draw_metadata_cache_store_for_component(
	       dbip, name, BOBOL_DRAW_CACHE_METADATA, 0, record, status);
}

extern "C" int
bobol_draw_metadata_cache_get(db_i *dbip,
				const char *name,
				BObolDrawMetadataRecord *record)
{
    return bobol_draw_metadata_cache_get_for_component(
	       dbip, name, BOBOL_DRAW_CACHE_METADATA, record);
}

extern "C" int
bobol_draw_metadata_cache_refresh(db_i *dbip,
				    const char *name,
				    BObolDrawCacheStatus *status)
{
    BObolDrawMetadataRecord record;
    bu_attribute_value_set avs = BU_AVS_INIT_ZERO;
    int haveAttrs = 0;
    int attrInt = 0;

    if (status)
	bobol_draw_cache_status_init(status);
    if (!dbip || !name)
	return BRLCAD_ERROR;

    directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
	if (status)
	    status->directoryFound = 0;
	return BRLCAD_ERROR;
    }

    bobol_draw_metadata_record_init(&record);
    bobol_draw_metadata_record_apply_directory(&record, dp);

    if (db5_get_attributes(dbip, &avs, dp) == 0)
	haveAttrs = 1;

    if (haveAttrs) {
	const char *region = bu_avs_get(&avs, "region");
	if (region && bu_str_true(region))
	    record.isRegion = 1;

	if (bobol_draw_metadata_int_attr(&avs, "region_id", &attrInt)) {
	    record.hasRegionId = 1;
	    record.regionId = attrInt;
	} else if (record.isRegion) {
	    record.hasRegionId = 1;
	    record.regionId = 0;
	}
	if (bobol_draw_metadata_int_attr(&avs, "air", &attrInt)) {
	    record.hasAircode = 1;
	    record.aircode = attrInt;
	}
	if (bobol_draw_metadata_int_attr(&avs, "los", &attrInt)) {
	    record.hasLos = 1;
	    record.los = attrInt;
	}
	if (bobol_draw_metadata_int_attr(&avs, "material_id", &attrInt)) {
	    record.hasMaterialId = 1;
	    record.materialId = attrInt;
	}

	const char *inherit = bu_avs_get(&avs, "inherit");
	if (inherit && inherit[0]) {
	    record.hasInherit = 1;
	    record.inherit = bu_str_true(inherit) ? 1 : 0;
	}

	(void)bobol_draw_metadata_color_attr(&avs, &record);
	bobol_draw_metadata_string_attr(&avs, "material_name",
					  record.materialName, sizeof(record.materialName),
					  &record.hasMaterialName);
	bobol_draw_metadata_string_attr(&avs, "shader",
					  record.shader, sizeof(record.shader), &record.hasShader);
    } else if (record.isRegion) {
	record.hasRegionId = 1;
	record.regionId = 0;
    }

    if (haveAttrs)
	bu_avs_free(&avs);

    return bobol_draw_metadata_cache_store(dbip, name, &record, status);
}

extern "C" int
bobol_draw_path_metadata_cache_status(db_i *dbip,
					const char *path,
					BObolDrawCacheStatus *status)
{
    if (!status || !dbip || !path)
	return BRLCAD_ERROR;
    bobol_draw_metadata_status_current_for_component(dbip, path,
	    BOBOL_DRAW_CACHE_PATH_METADATA, 1, status);
    return BRLCAD_OK;
}

extern "C" int
bobol_draw_path_metadata_cache_invalidate(db_i *dbip,
	const char *path,
	BObolDrawCacheStatus *status)
{
    return bobol_draw_metadata_cache_invalidate_for_component(
	       dbip, path, BOBOL_DRAW_CACHE_PATH_METADATA, 1, status);
}

extern "C" int
bobol_draw_path_metadata_cache_invalidate_object(
    db_i *dbip,
    const char *name,
    BObolDrawCacheStatus *status)
{
    char **keys = NULL;
    int nkeys = 0;
    int cleared = 0;
    int opened = 0;

    if (status)
	bobol_draw_cache_status_init(status);
    if (!dbip || !name || !name[0])
	return BRLCAD_ERROR;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    opened = bobol_draw_cache_open(&handle, dbip);
    if (opened) {
	nkeys = bu_cache_keys(&keys, handle.cache);
	std::vector<std::string> clearKeys;
	struct bu_cache_txn *readTxn = NULL;
	for (int i = 0; i < nkeys; i++) {
	    if (!bobol_draw_cache_key_component_is(keys[i],
		BOBOL_DRAW_CACHE_PATH_METADATA_INDEX))
		continue;

	    void *data = NULL;
	    size_t dsize = bu_cache_get(&data, keys[i], handle.cache,
		&readTxn);
	    const char *path = (const char *)data;
	    if (data && dsize > 0 && ((const char *)data)[dsize - 1] == '\0' &&
		bobol_draw_path_contains_object_name(path, name)) {
		char pathMetadataKey[BU_CACHE_KEY_MAXLEN] = {0};
		bobol_draw_cache_key(pathMetadataKey, path,
				       BOBOL_DRAW_CACHE_PATH_METADATA);
		if (pathMetadataKey[0])
		    clearKeys.emplace_back(pathMetadataKey);
		clearKeys.emplace_back(keys[i]);
		cleared = 1;
	    }
	}
	bu_cache_get_done(&readTxn);
	struct bu_cache_txn *writeTxn = NULL;
	for (const std::string &clearKey : clearKeys)
	    bu_cache_clear(clearKey.c_str(), handle.cache, &writeTxn);
	if (writeTxn &&
	    bu_cache_write_commit(handle.cache, &writeTxn) != BRLCAD_OK) {
	    opened = 0;
	    cleared = 0;
	}
	bobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    if (nkeys)
	bu_argv_free((size_t)nkeys, keys);

    if (status) {
	status->directoryFound = opened ? 1 : 0;
	status->clearedCacheEntry = cleared;
	status->hasCachedPayload = 0;
    }
    return opened ? BRLCAD_OK : BRLCAD_ERROR;
}

extern "C" int
bobol_draw_path_metadata_cache_store(db_i *dbip,
				       const char *path,
				       const BObolDrawMetadataRecord *record,
				       BObolDrawCacheStatus *status)
{
    return bobol_draw_metadata_cache_store_for_component(
	       dbip, path, BOBOL_DRAW_CACHE_PATH_METADATA, 1, record, status);
}

extern "C" int
bobol_draw_path_metadata_cache_get(db_i *dbip,
				     const char *path,
				     BObolDrawMetadataRecord *record)
{
    return bobol_draw_metadata_cache_get_for_component(
	       dbip, path, BOBOL_DRAW_CACHE_PATH_METADATA, record);
}

extern "C" int
bobol_draw_path_metadata_cache_refresh(db_i *dbip,
	const char *path,
	BObolDrawCacheStatus *status)
{
    BObolDrawMetadataRecord record;
    struct db_tree_state state;
    struct db_full_path fullPath;

    if (status)
	bobol_draw_cache_status_init(status);
    if (!dbip || !path)
	return BRLCAD_ERROR;

    db_init_db_tree_state(&state, dbip);
    db_full_path_init(&fullPath);
    if (db_follow_path_for_state(&state, &fullPath, path, LOOKUP_QUIET) < 0 ||
	!DB_FULL_PATH_CUR_DIR(&fullPath)) {
	if (status)
	    status->directoryFound = 0;
	db_free_full_path(&fullPath);
	db_free_db_tree_state(&state);
	return BRLCAD_ERROR;
    }

    bobol_draw_metadata_record_from_tree_state(&record, &state,
	    DB_FULL_PATH_CUR_DIR(&fullPath));

    db_free_full_path(&fullPath);
    db_free_db_tree_state(&state);

    return bobol_draw_path_metadata_cache_store(dbip, path, &record,
	    status);
}

static int
bobol_draw_manifest_matrix_valid(const fastf_t matrix[16])
{
    if (!matrix)
	return 0;
    for (size_t i = 0; i < 16; i++) {
	if (!std::isfinite(matrix[i]))
	    return 0;
    }
    return 1;
}

static int
bobol_draw_manifest_boolean_valid(int operation)
{
    return operation == DB_OP_UNION || operation == DB_OP_SUBTRACT ||
	operation == DB_OP_INTERSECT;
}

static int
bobol_draw_manifest_profile_valid_or_empty(uint64_t occurrenceCount,
	uint64_t uniqueAssetCount, uint64_t encodedSourceBytes,
	uint64_t largestAssetBytes)
{
    const bool empty = !uniqueAssetCount && !encodedSourceBytes &&
	!largestAssetBytes;
    if (empty)
	return 1;
    return occurrenceCount > 0 && uniqueAssetCount > 0 &&
	uniqueAssetCount <= occurrenceCount && encodedSourceBytes > 0 &&
	largestAssetBytes > 0 && largestAssetBytes <= encodedSourceBytes;
}

static int
bobol_draw_manifest_oriented_bounds_valid(
	const BObolDrawManifestOccurrence &occurrence)
{
    if (!occurrence.orientedBoundsValid)
	return occurrence.orientedBoundsValid == 0;
    if (occurrence.orientedBoundsValid != 1)
	return 0;

    BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    VMOVE(hierarchy.quantization_min, occurrence.boundsMin);
    VMOVE(hierarchy.quantization_max, occurrence.boundsMax);
    hierarchy.oriented_bounds_valid = 1;
    for (size_t corner = 0; corner < 8; ++corner)
	VMOVE(hierarchy.oriented_bounds[corner],
	    occurrence.orientedBounds[corner]);
    return bobol_mesh_lod_oriented_bounds_validate(&hierarchy);
}

extern "C" void
bobol_draw_manifest_init(BObolDrawManifest *manifest)
{
    if (!manifest)
	return;
    manifest->coverageBoundsValid = 0;
    VSETALL(manifest->coverageBoundsMin, 0.0);
    VSETALL(manifest->coverageBoundsMax, 0.0);
    manifest->uniqueAssetCount = 0;
    manifest->encodedSourceBytes = 0;
    manifest->largestAssetBytes = 0;
    manifest->occurrenceCount = 0;
    manifest->occurrences = NULL;
}

extern "C" void
bobol_draw_manifest_free(BObolDrawManifest *manifest)
{
    if (!manifest)
	return;

    if (manifest->occurrences) {
	for (size_t i = 0; i < manifest->occurrenceCount; i++) {
	    if (manifest->occurrences[i].path)
		bu_free(manifest->occurrences[i].path,
		    "bobol draw manifest path");
	    if (manifest->occurrences[i].sourceName)
		bu_free(manifest->occurrences[i].sourceName,
		    "bobol draw manifest source name");
	    if (manifest->occurrences[i].meshAssetPath)
		bu_free(manifest->occurrences[i].meshAssetPath,
		    "bobol draw manifest mesh asset path");
	    if (manifest->occurrences[i].meshAssetName)
		bu_free(manifest->occurrences[i].meshAssetName,
		    "bobol draw manifest mesh asset name");
	}
    }
    if (manifest->occurrences)
	bu_free(manifest->occurrences, "bobol draw manifest occurrences");
    bobol_draw_manifest_init(manifest);
}

static int
bobol_draw_manifest_cache_store_with_provider(db_i *dbip,
	const char *rootPath, const BObolDrawManifest *manifest,
	BObolDrawManifestOccurrenceProvider provider, void *userData)
{
    if (!dbip || !rootPath || !rootPath[0] || !manifest ||
	!manifest->occurrenceCount || !provider)
	return BRLCAD_ERROR;

    if ((manifest->coverageBoundsValid != 0 &&
	 manifest->coverageBoundsValid != 1) ||
	(manifest->coverageBoundsValid &&
	 !bobol_draw_proxy_bbox_valid(manifest->coverageBoundsMin,
	     manifest->coverageBoundsMax)) ||
	!bobol_draw_manifest_profile_valid_or_empty(
	    static_cast<uint64_t>(manifest->occurrenceCount),
	    manifest->uniqueAssetCount, manifest->encodedSourceBytes,
	    manifest->largestAssetBytes))
	return BRLCAD_ERROR;

    const size_t rootLength = strlen(rootPath);
    if (rootLength > UINT32_MAX)
	return BRLCAD_ERROR;
    for (size_t i = 0; i < manifest->occurrenceCount; i++) {
	BObolDrawManifestOccurrence occurrence;
	memset(&occurrence, 0, sizeof(occurrence));
	if (!provider(i, &occurrence, userData))
	    return BRLCAD_ERROR;
	if (!occurrence.path || !occurrence.path[0] || !occurrence.sourceName)
	    return BRLCAD_ERROR;
	const size_t pathLength = strlen(occurrence.path);
	const size_t sourceNameLength = strlen(occurrence.sourceName);
	const size_t meshAssetPathLength = occurrence.meshAssetPath ?
	    strlen(occurrence.meshAssetPath) : 0;
	const size_t meshAssetNameLength = occurrence.meshAssetName ?
	    strlen(occurrence.meshAssetName) : 0;
	if (pathLength > UINT32_MAX || sourceNameLength > UINT32_MAX ||
	    meshAssetPathLength > UINT32_MAX ||
	    meshAssetNameLength > UINT32_MAX ||
	    (occurrence.metadataValid != 0 && occurrence.metadataValid != 1) ||
	    (occurrence.sourceMeshRequestValid != 0 &&
	     occurrence.sourceMeshRequestValid != 1) ||
	    occurrence.meshAssetKind < BOBOL_DRAW_CACHE_MESH_ASSET_UNKNOWN ||
	    occurrence.meshAssetKind > BOBOL_DRAW_CACHE_MESH_ASSET_BREP ||
	    !std::isfinite(occurrence.meshAssetTessellationAbsTol) ||
	    !std::isfinite(occurrence.meshAssetTessellationRelTol) ||
	    !std::isfinite(occurrence.meshAssetTessellationNormTol) ||
	    occurrence.meshAssetTessellationAbsTol < 0.0 ||
	    occurrence.meshAssetTessellationRelTol < 0.0 ||
	    occurrence.meshAssetTessellationNormTol < 0.0 ||
	    !bobol_draw_manifest_oriented_bounds_valid(occurrence) ||
	    !bobol_draw_manifest_boolean_valid(occurrence.booleanOperation) ||
	    !bobol_draw_manifest_matrix_valid(occurrence.localMatrix) ||
	    !bobol_draw_proxy_bbox_valid(occurrence.boundsMin,
		occurrence.boundsMax) ||
	    (occurrence.sourceMeshRequestValid &&
	     (!meshAssetPathLength || !meshAssetNameLength ||
	      !bobol_draw_manifest_matrix_valid(
		  occurrence.meshAssetMatrix) ||
	      !bobol_draw_proxy_bbox_valid(occurrence.meshAssetBoundsMin,
		  occurrence.meshAssetBoundsMax))) ||
	    pathLength > SIZE_MAX - sizeof(BObolDrawManifestOccurrenceDiskHeader) ||
	    sourceNameLength > SIZE_MAX - sizeof(BObolDrawManifestOccurrenceDiskHeader) ||
	    meshAssetPathLength > SIZE_MAX - sizeof(BObolDrawManifestOccurrenceDiskHeader) ||
	    meshAssetNameLength > SIZE_MAX - sizeof(BObolDrawManifestOccurrenceDiskHeader))
	    return BRLCAD_ERROR;
    }

    char key[BU_CACHE_KEY_MAXLEN] = {0};
    bobol_draw_cache_key(key, rootPath, BOBOL_DRAW_CACHE_MANIFEST);
    if (!key[0])
	return BRLCAD_ERROR;

    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    if (!bobol_draw_cache_open(&handle, dbip)) {
	bu_semaphore_release(sem);
	return BRLCAD_ERROR;
    }

    const uint64_t fingerprint = bobol_draw_cache_database_fingerprint(dbip);
    std::vector<unsigned char> chunk;
    uint32_t chunkIndex = 0;
    uint32_t chunkOccurrenceCount = 0;
    int status = BRLCAD_OK;
    auto startChunk = [&]() {
	chunk.clear();
	chunk.resize(sizeof(BObolDrawManifestChunkDiskHeader), 0);
	chunkOccurrenceCount = 0;
    };
    auto finishChunk = [&]() -> int {
	if (!chunkOccurrenceCount)
	    return BRLCAD_OK;
	BObolDrawManifestChunkDiskHeader header;
	memset(&header, 0, sizeof(header));
	header.magic = BOBOL_DRAW_MANIFEST_CHUNK_DISK_MAGIC;
	header.version = BOBOL_DRAW_MANIFEST_CHUNK_DISK_VERSION;
	header.databaseFingerprint = fingerprint;
	header.chunkIndex = chunkIndex;
	header.occurrenceCount = chunkOccurrenceCount;
	memcpy(chunk.data(), &header, sizeof(header));
	std::string chunkIdentity = std::string(rootPath) +
	    "|manifest-chunk-v3|" + std::to_string(chunkIndex);
	char chunkKey[BU_CACHE_KEY_MAXLEN] = {0};
	bobol_draw_cache_key(chunkKey, chunkIdentity.c_str(),
	    BOBOL_DRAW_CACHE_MANIFEST);
	if (!chunkKey[0] ||
	    bu_cache_write(chunk.data(), chunk.size(), chunkKey, handle.cache,
		NULL) != chunk.size())
	    return BRLCAD_ERROR;
	chunkIndex++;
	return BRLCAD_OK;
    };
    try {
	startChunk();
    } catch (const std::bad_alloc &) {
	status = BRLCAD_ERROR;
    }
    for (size_t i = 0; status == BRLCAD_OK &&
	i < manifest->occurrenceCount; i++) {
	BObolDrawManifestOccurrence occurrence;
	memset(&occurrence, 0, sizeof(occurrence));
	if (!provider(i, &occurrence, userData)) {
	    status = BRLCAD_ERROR;
	    break;
	}
	const size_t pathLength = strlen(occurrence.path);
	const size_t sourceNameLength = strlen(occurrence.sourceName);
	const size_t meshAssetPathLength = occurrence.meshAssetPath ?
	    strlen(occurrence.meshAssetPath) : 0;
	const size_t meshAssetNameLength = occurrence.meshAssetName ?
	    strlen(occurrence.meshAssetName) : 0;
	const size_t orientedBoundsBytes = occurrence.orientedBoundsValid ?
	    sizeof(occurrence.orientedBounds) : 0;
	size_t occurrenceSize = sizeof(BObolDrawManifestOccurrenceDiskHeader);
	const size_t fieldSizes[] = {
	    pathLength, sourceNameLength, meshAssetPathLength,
	    meshAssetNameLength, orientedBoundsBytes
	};
	bool sizeValid = true;
	for (const size_t fieldSize : fieldSizes) {
	    if (fieldSize > SIZE_MAX - occurrenceSize) {
		sizeValid = false;
		break;
	    }
	    occurrenceSize += fieldSize;
	}
	if (!sizeValid || chunk.size() > SIZE_MAX - occurrenceSize) {
	    status = BRLCAD_ERROR;
	    break;
	}
	if (chunkOccurrenceCount &&
	    chunk.size() + occurrenceSize > manifestChunkTargetBytes) {
	    if (finishChunk() != BRLCAD_OK) {
		status = BRLCAD_ERROR;
		break;
	    }
	    try {
		startChunk();
	    } catch (const std::bad_alloc &) {
		status = BRLCAD_ERROR;
		break;
	    }
	}
	BObolDrawManifestOccurrenceDiskHeader occurrenceHeader;
	memset(&occurrenceHeader, 0, sizeof(occurrenceHeader));
	occurrenceHeader.pathLength = static_cast<uint32_t>(pathLength);
	occurrenceHeader.sourceNameLength =
	    static_cast<uint32_t>(sourceNameLength);
	occurrenceHeader.meshAssetPathLength =
	    static_cast<uint32_t>(meshAssetPathLength);
	occurrenceHeader.meshAssetNameLength =
	    static_cast<uint32_t>(meshAssetNameLength);
	occurrenceHeader.booleanOperation = occurrence.booleanOperation;
	occurrenceHeader.occurrenceIndex = occurrence.occurrenceIndex;
	occurrenceHeader.metadataValid = occurrence.metadataValid ? 1u : 0u;
	occurrenceHeader.sourceMeshRequestValid =
	    occurrence.sourceMeshRequestValid ? 1u : 0u;
	occurrenceHeader.meshAssetKind =
	    static_cast<uint32_t>(occurrence.meshAssetKind);
	occurrenceHeader.orientedBoundsValid =
	    occurrence.orientedBoundsValid ? 1u : 0u;
	occurrenceHeader.meshAssetContentHash =
	    occurrence.meshAssetContentHash;
	occurrenceHeader.sourceFaceCount = occurrence.sourceFaceCount;
	occurrenceHeader.sourcePointCount = occurrence.sourcePointCount;
	occurrenceHeader.meshAssetTessellationAbsTol =
	    occurrence.meshAssetTessellationAbsTol;
	occurrenceHeader.meshAssetTessellationRelTol =
	    occurrence.meshAssetTessellationRelTol;
	occurrenceHeader.meshAssetTessellationNormTol =
	    occurrence.meshAssetTessellationNormTol;
	memcpy(occurrenceHeader.localMatrix, occurrence.localMatrix,
	    sizeof(occurrenceHeader.localMatrix));
	memcpy(occurrenceHeader.meshAssetMatrix,
	    occurrence.meshAssetMatrix,
	    sizeof(occurrenceHeader.meshAssetMatrix));
	VMOVE(occurrenceHeader.boundsMin, occurrence.boundsMin);
	VMOVE(occurrenceHeader.boundsMax, occurrence.boundsMax);
	if (occurrence.sourceMeshRequestValid) {
	    VMOVE(occurrenceHeader.meshAssetBoundsMin,
		occurrence.meshAssetBoundsMin);
	    VMOVE(occurrenceHeader.meshAssetBoundsMax,
		occurrence.meshAssetBoundsMax);
	}
	if (occurrence.metadataValid)
	    bobol_draw_metadata_to_disk(&occurrenceHeader.metadata,
		&occurrence.metadata);
	const size_t appendOffset = chunk.size();
	try {
	    chunk.resize(appendOffset + occurrenceSize);
	} catch (const std::bad_alloc &) {
	    status = BRLCAD_ERROR;
	    break;
	}
	unsigned char *write = chunk.data() + appendOffset;
	memcpy(write, &occurrenceHeader, sizeof(occurrenceHeader));
	write += sizeof(occurrenceHeader);
	memcpy(write, occurrence.path, pathLength);
	write += pathLength;
	memcpy(write, occurrence.sourceName, sourceNameLength);
	write += sourceNameLength;
	if (meshAssetPathLength) {
	    memcpy(write, occurrence.meshAssetPath, meshAssetPathLength);
	    write += meshAssetPathLength;
	}
	if (meshAssetNameLength) {
	    memcpy(write, occurrence.meshAssetName, meshAssetNameLength);
	    write += meshAssetNameLength;
	}
	if (orientedBoundsBytes)
	    memcpy(write, occurrence.orientedBounds, orientedBoundsBytes);
	chunkOccurrenceCount++;
    }

    if (status == BRLCAD_OK &&
	(finishChunk() != BRLCAD_OK || !chunkIndex))
	status = BRLCAD_ERROR;

    BObolDrawManifestChunkedDiskHeader descriptor;
    memset(&descriptor, 0, sizeof(descriptor));
    descriptor.magic = BOBOL_DRAW_MANIFEST_CHUNKED_DISK_MAGIC;
    descriptor.version = BOBOL_DRAW_MANIFEST_CHUNK_DISK_VERSION;
    descriptor.databaseFingerprint = fingerprint;
    descriptor.occurrenceCount = static_cast<uint64_t>(manifest->occurrenceCount);
    descriptor.uniqueAssetCount = manifest->uniqueAssetCount;
    descriptor.encodedSourceBytes = manifest->encodedSourceBytes;
    descriptor.largestAssetBytes = manifest->largestAssetBytes;
    descriptor.rootPathLength = static_cast<uint32_t>(rootLength);
    descriptor.coverageBoundsValid = manifest->coverageBoundsValid ? 1u : 0u;
    descriptor.chunkCount = chunkIndex;
    if (manifest->coverageBoundsValid) {
	VMOVE(descriptor.coverageBoundsMin, manifest->coverageBoundsMin);
	VMOVE(descriptor.coverageBoundsMax, manifest->coverageBoundsMax);
    }
    if (status == BRLCAD_OK) {
	std::vector<unsigned char> disk(sizeof(descriptor) + rootLength);
	memcpy(disk.data(), &descriptor, sizeof(descriptor));
	memcpy(disk.data() + sizeof(descriptor), rootPath, rootLength);
	if (bu_cache_write(disk.data(), disk.size(), key, handle.cache, NULL) !=
	    disk.size())
	    status = BRLCAD_ERROR;
    }
	bobol_draw_cache_close(&handle);
    bu_semaphore_release(sem);
    return status;
}

struct BObolDrawManifestArrayProvider {
    const BObolDrawManifest *manifest = nullptr;
};

static int
bobol_draw_manifest_array_provider(size_t occurrenceIndex,
	BObolDrawManifestOccurrence *occurrence, void *userData)
{
    const BObolDrawManifestArrayProvider *provider =
	static_cast<const BObolDrawManifestArrayProvider *>(userData);
    if (!provider || !provider->manifest || !occurrence ||
	!provider->manifest->occurrences ||
	occurrenceIndex >= provider->manifest->occurrenceCount)
	return 0;
    *occurrence = provider->manifest->occurrences[occurrenceIndex];
    return 1;
}

extern "C" int
bobol_draw_manifest_cache_store(db_i *dbip, const char *rootPath,
	const BObolDrawManifest *manifest)
{
    if (!manifest || !manifest->occurrences)
	return BRLCAD_ERROR;
    BObolDrawManifestArrayProvider provider;
    provider.manifest = manifest;
    return bobol_draw_manifest_cache_store_with_provider(dbip, rootPath,
	manifest, bobol_draw_manifest_array_provider, &provider);
}

extern "C" int
bobol_draw_manifest_cache_store_visit(db_i *dbip, const char *rootPath,
	const BObolDrawManifest *manifest,
	BObolDrawManifestOccurrenceProvider occurrence, void *userData)
{
	return bobol_draw_manifest_cache_store_with_provider(dbip, rootPath,
	manifest, occurrence, userData);
}

static int
bobol_draw_manifest_stream_chunk(const unsigned char *bytes, size_t dataSize,
	size_t occurrenceCount, size_t *visited,
	BObolDrawManifestOccurrenceCallback visit, void *userData)
{
    if (!bytes || !visited || !visit)
	return BRLCAD_ERROR;
    size_t offset = 0;
    for (size_t i = 0; i < occurrenceCount; i++) {
	if (offset > dataSize ||
	    dataSize - offset < sizeof(BObolDrawManifestOccurrenceDiskHeader))
	    return BRLCAD_ERROR;
	BObolDrawManifestOccurrenceDiskHeader header;
	memcpy(&header, bytes + offset, sizeof(header));
	offset += sizeof(header);
	const size_t pathLength = header.pathLength;
	const size_t sourceNameLength = header.sourceNameLength;
	const size_t assetPathLength = header.meshAssetPathLength;
	const size_t assetNameLength = header.meshAssetNameLength;
	const size_t orientedBoundsBytes = header.orientedBoundsValid ?
	    sizeof(point_t) * 8 : 0;
	if (!pathLength || pathLength > dataSize - offset ||
	    sourceNameLength > dataSize - offset - pathLength ||
	    assetPathLength > dataSize - offset - pathLength - sourceNameLength ||
	    assetNameLength > dataSize - offset - pathLength - sourceNameLength -
		assetPathLength ||
	    orientedBoundsBytes > dataSize - offset - pathLength -
		sourceNameLength - assetPathLength - assetNameLength ||
	    header.metadataValid > 1 ||
	    header.sourceMeshRequestValid > 1 ||
	    header.orientedBoundsValid > 1 ||
	    header.meshAssetKind > BOBOL_DRAW_CACHE_MESH_ASSET_BREP ||
	    !std::isfinite(header.meshAssetTessellationAbsTol) ||
	    !std::isfinite(header.meshAssetTessellationRelTol) ||
	    !std::isfinite(header.meshAssetTessellationNormTol) ||
	    header.meshAssetTessellationAbsTol < 0.0 ||
	    header.meshAssetTessellationRelTol < 0.0 ||
	    header.meshAssetTessellationNormTol < 0.0 ||
	    (header.sourceMeshRequestValid &&
	     (!assetPathLength || !assetNameLength ||
	      !bobol_draw_manifest_matrix_valid(header.meshAssetMatrix) ||
	      !bobol_draw_proxy_bbox_valid(header.meshAssetBoundsMin,
		  header.meshAssetBoundsMax))) ||
	    !bobol_draw_manifest_boolean_valid(header.booleanOperation) ||
	    !bobol_draw_manifest_matrix_valid(header.localMatrix) ||
	    !bobol_draw_proxy_bbox_valid(header.boundsMin, header.boundsMax) ||
	    memchr(bytes + offset, '\0', pathLength) ||
	    memchr(bytes + offset + pathLength, '\0', sourceNameLength) ||
	    (assetPathLength && memchr(bytes + offset + pathLength +
		 sourceNameLength, '\0', assetPathLength)) ||
	    (assetNameLength && memchr(bytes + offset + pathLength +
		 sourceNameLength + assetPathLength, '\0', assetNameLength)) ||
	    (header.metadataValid && !bobol_draw_metadata_disk_valid(
		&header.metadata, sizeof(header.metadata))))
	    return BRLCAD_ERROR;
	std::string path(reinterpret_cast<const char *>(bytes + offset), pathLength);
	offset += pathLength;
	std::string sourceName(reinterpret_cast<const char *>(bytes + offset),
	    sourceNameLength);
	offset += sourceNameLength;
	std::string assetPath;
	std::string assetName;
	if (assetPathLength) {
	    assetPath.assign(reinterpret_cast<const char *>(bytes + offset),
		assetPathLength);
	    offset += assetPathLength;
	}
	if (assetNameLength) {
	    assetName.assign(reinterpret_cast<const char *>(bytes + offset),
		assetNameLength);
	    offset += assetNameLength;
	}
	BObolDrawManifestOccurrence occurrence;
	memset(&occurrence, 0, sizeof(occurrence));
	occurrence.path = const_cast<char *>(path.c_str());
	occurrence.sourceName = const_cast<char *>(sourceName.c_str());
	occurrence.meshAssetPath = assetPathLength ?
	    const_cast<char *>(assetPath.c_str()) : NULL;
	occurrence.meshAssetName = assetNameLength ?
	    const_cast<char *>(assetName.c_str()) : NULL;
	occurrence.booleanOperation = header.booleanOperation;
	occurrence.occurrenceIndex = header.occurrenceIndex;
	occurrence.metadataValid = header.metadataValid ? 1 : 0;
	occurrence.sourceMeshRequestValid = header.sourceMeshRequestValid ? 1 : 0;
	occurrence.orientedBoundsValid = header.orientedBoundsValid ? 1 : 0;
	occurrence.meshAssetKind = static_cast<int>(header.meshAssetKind);
	occurrence.meshAssetContentHash = header.meshAssetContentHash;
	occurrence.sourceFaceCount = header.sourceFaceCount;
	occurrence.sourcePointCount = header.sourcePointCount;
	occurrence.meshAssetTessellationAbsTol = header.meshAssetTessellationAbsTol;
	occurrence.meshAssetTessellationRelTol = header.meshAssetTessellationRelTol;
	occurrence.meshAssetTessellationNormTol = header.meshAssetTessellationNormTol;
	memcpy(occurrence.localMatrix, header.localMatrix, sizeof(occurrence.localMatrix));
	memcpy(occurrence.meshAssetMatrix, header.meshAssetMatrix,
	    sizeof(occurrence.meshAssetMatrix));
	VMOVE(occurrence.boundsMin, header.boundsMin);
	VMOVE(occurrence.boundsMax, header.boundsMax);
	if (occurrence.sourceMeshRequestValid) {
	    VMOVE(occurrence.meshAssetBoundsMin, header.meshAssetBoundsMin);
	    VMOVE(occurrence.meshAssetBoundsMax, header.meshAssetBoundsMax);
	}
	if (occurrence.orientedBoundsValid) {
	    memcpy(occurrence.orientedBounds, bytes + offset,
		orientedBoundsBytes);
	    offset += orientedBoundsBytes;
	    if (!bobol_draw_manifest_oriented_bounds_valid(occurrence))
		return BRLCAD_ERROR;
	}
	bobol_draw_metadata_from_disk(&occurrence.metadata, &header.metadata);
	if (!visit(&occurrence, *visited, userData))
	    return BRLCAD_ERROR;
	(*visited)++;
    }
    return offset == dataSize ? BRLCAD_OK : BRLCAD_ERROR;
}

extern "C" int
bobol_draw_manifest_cache_stream(db_i *dbip, const char *rootPath,
	BObolDrawManifestBeginCallback begin,
	BObolDrawManifestOccurrenceCallback visit, void *userData)
{
    if (!dbip || !rootPath || !rootPath[0])
	return BRLCAD_ERROR;

    char key[BU_CACHE_KEY_MAXLEN] = {0};
    bobol_draw_cache_key(key, rootPath, BOBOL_DRAW_CACHE_MANIFEST);
    if (!key[0])
	return BRLCAD_ERROR;

    void *data = NULL;
    size_t dataSize = 0;
    bu_cache_txn *transaction = NULL;
    int sem = bobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BObolDrawCacheHandle handle;
    if (bobol_draw_cache_open(&handle, dbip)) {
	/* A retained read transaction returns LMDB's read-only mapping rather
	 * than allocating and copying a potentially hundreds-of-megabytes value.
	 * Keep the draw-cache lifecycle lock until the callbacks finish so an
	 * explicit cache clear cannot retire the mapped environment underneath
	 * them.  PoP mesh data uses its own cache and remains fully concurrent. */
	dataSize = bu_cache_get(&data, key, handle.cache, &transaction);
    }
    if (!data || dataSize < sizeof(BObolDrawManifestChunkedDiskHeader)) {
	if (getenv("BOBOL_DRAW_TIMING"))
	    bu_log("[obol-timing] manifest cache record unavailable for %s "
		"(bytes=%zu)\n", rootPath, dataSize);
	bu_cache_get_done(&transaction);
	bobol_draw_cache_close(&handle);
	bu_semaphore_release(sem);
	return BRLCAD_ERROR;
    }

    /* New records are a compact descriptor plus independently mapped,
     * bounded occurrence chunks.  Do not reassemble them: a 150k assembly
     * must remain streamable under the same memory policy as its geometry. */
    if (dataSize >= sizeof(BObolDrawManifestChunkedDiskHeader)) {
	BObolDrawManifestChunkedDiskHeader descriptor;
	memcpy(&descriptor, data, sizeof(descriptor));
	if (descriptor.magic == BOBOL_DRAW_MANIFEST_CHUNKED_DISK_MAGIC) {
	    const size_t rootLength = strlen(rootPath);
	    const uint64_t fingerprint = bobol_draw_cache_database_fingerprint(dbip);
	    const unsigned char *descriptorBytes =
		static_cast<const unsigned char *>(data);
	    const int descriptorValid = rootLength <= UINT32_MAX &&
		descriptor.version == BOBOL_DRAW_MANIFEST_CHUNK_DISK_VERSION &&
		descriptor.databaseFingerprint == fingerprint &&
		descriptor.rootPathLength == rootLength &&
		descriptor.chunkCount && descriptor.occurrenceCount &&
		bobol_draw_manifest_profile_valid_or_empty(
		    descriptor.occurrenceCount,
		    descriptor.uniqueAssetCount,
		    descriptor.encodedSourceBytes,
		    descriptor.largestAssetBytes) &&
		descriptor.coverageBoundsValid <= 1 &&
		(!descriptor.coverageBoundsValid ||
		 bobol_draw_proxy_bbox_valid(descriptor.coverageBoundsMin,
		     descriptor.coverageBoundsMax)) &&
		dataSize == sizeof(descriptor) + rootLength &&
		!memcmp(descriptorBytes + sizeof(descriptor), rootPath, rootLength);
	    if (!descriptorValid) {
		bu_cache_get_done(&transaction);
		bobol_draw_cache_close(&handle);
		bu_semaphore_release(sem);
		return BRLCAD_ERROR;
	    }
	    BObolDrawManifest description;
	    bobol_draw_manifest_init(&description);
	    description.coverageBoundsValid = descriptor.coverageBoundsValid ? 1 : 0;
	    if (description.coverageBoundsValid) {
		VMOVE(description.coverageBoundsMin, descriptor.coverageBoundsMin);
		VMOVE(description.coverageBoundsMax, descriptor.coverageBoundsMax);
	    }
	    description.occurrenceCount = static_cast<size_t>(descriptor.occurrenceCount);
	    description.uniqueAssetCount = descriptor.uniqueAssetCount;
	    description.encodedSourceBytes = descriptor.encodedSourceBytes;
	    description.largestAssetBytes = descriptor.largestAssetBytes;
	    if (begin && !begin(&description, userData)) {
		bu_cache_get_done(&transaction);
		bobol_draw_cache_close(&handle);
		bu_semaphore_release(sem);
		return BRLCAD_ERROR;
	    }
	    if (!visit) {
		bu_cache_get_done(&transaction);
		bobol_draw_cache_close(&handle);
		bu_semaphore_release(sem);
		return BRLCAD_OK;
	    }
	    bu_cache_get_done(&transaction);
	    size_t visited = 0;
	    int valid = 1;
	    for (uint32_t i = 0; valid && i < descriptor.chunkCount; i++) {
		std::string chunkIdentity = std::string(rootPath) +
		    "|manifest-chunk-v3|" + std::to_string(i);
		char chunkKey[BU_CACHE_KEY_MAXLEN] = {0};
		bobol_draw_cache_key(chunkKey, chunkIdentity.c_str(),
		    BOBOL_DRAW_CACHE_MANIFEST);
		void *chunkData = NULL;
		size_t chunkSize = 0;
		bu_cache_txn *chunkTransaction = NULL;
		if (!chunkKey[0] ||
		    !(chunkSize = bu_cache_get(&chunkData, chunkKey, handle.cache,
			&chunkTransaction)) ||
		    chunkSize < sizeof(BObolDrawManifestChunkDiskHeader)) {
		    valid = 0;
		} else {
		    BObolDrawManifestChunkDiskHeader chunkHeader;
		    memcpy(&chunkHeader, chunkData, sizeof(chunkHeader));
		    if (chunkHeader.magic != BOBOL_DRAW_MANIFEST_CHUNK_DISK_MAGIC ||
			chunkHeader.version != BOBOL_DRAW_MANIFEST_CHUNK_DISK_VERSION ||
			chunkHeader.databaseFingerprint != fingerprint ||
			chunkHeader.chunkIndex != i || !chunkHeader.occurrenceCount ||
			chunkHeader.occurrenceCount > descriptor.occurrenceCount - visited ||
			bobol_draw_manifest_stream_chunk(
			    static_cast<const unsigned char *>(chunkData) +
				sizeof(chunkHeader),
			    chunkSize - sizeof(chunkHeader),
			    chunkHeader.occurrenceCount, &visited, visit,
			    userData) != BRLCAD_OK)
			valid = 0;
		}
		bu_cache_get_done(&chunkTransaction);
	    }
	    bobol_draw_cache_close(&handle);
	    bu_semaphore_release(sem);
	    return valid && visited == description.occurrenceCount ?
		BRLCAD_OK : BRLCAD_ERROR;
	}
    }

    bu_cache_get_done(&transaction);
    bobol_draw_cache_close(&handle);
    bu_semaphore_release(sem);
    return BRLCAD_ERROR;
}

struct BObolDrawManifestDescriptionContext {
    BObolDrawManifest *description = NULL;
    int captured = 0;
};

static int
bobol_draw_manifest_describe_begin(const BObolDrawManifest *description,
	void *userData)
{
    BObolDrawManifestDescriptionContext *context =
	static_cast<BObolDrawManifestDescriptionContext *>(userData);
    if (!context || !context->description || !description)
	return 0;
    context->description->coverageBoundsValid =
	description->coverageBoundsValid;
    VMOVE(context->description->coverageBoundsMin,
	description->coverageBoundsMin);
    VMOVE(context->description->coverageBoundsMax,
	description->coverageBoundsMax);
    context->description->occurrenceCount = description->occurrenceCount;
    context->description->uniqueAssetCount = description->uniqueAssetCount;
    context->description->encodedSourceBytes = description->encodedSourceBytes;
    context->description->largestAssetBytes = description->largestAssetBytes;
    context->description->occurrences = NULL;
    context->captured = 1;
    return 1;
}

extern "C" int
bobol_draw_manifest_cache_describe(db_i *dbip, const char *rootPath,
	BObolDrawManifest *description)
{
    if (!description)
	return BRLCAD_ERROR;
    bobol_draw_manifest_init(description);
    BObolDrawManifestDescriptionContext context;
    context.description = description;
    const int status = bobol_draw_manifest_cache_stream(dbip, rootPath,
	bobol_draw_manifest_describe_begin, NULL, &context);
    return status == BRLCAD_OK && context.captured ?
	BRLCAD_OK : BRLCAD_ERROR;
}

struct BObolDrawManifestCollector {
    BObolDrawManifest manifest;
};

static int
bobol_draw_manifest_collect_begin(const BObolDrawManifest *description,
	void *userData)
{
    BObolDrawManifestCollector *collector =
	static_cast<BObolDrawManifestCollector *>(userData);
    if (!collector || !description || !description->occurrenceCount)
	return 0;
    collector->manifest.coverageBoundsValid =
	description->coverageBoundsValid;
    VMOVE(collector->manifest.coverageBoundsMin,
	description->coverageBoundsMin);
    VMOVE(collector->manifest.coverageBoundsMax,
	description->coverageBoundsMax);
    collector->manifest.occurrenceCount = description->occurrenceCount;
    collector->manifest.uniqueAssetCount = description->uniqueAssetCount;
    collector->manifest.encodedSourceBytes = description->encodedSourceBytes;
    collector->manifest.largestAssetBytes = description->largestAssetBytes;
    collector->manifest.occurrences =
	static_cast<BObolDrawManifestOccurrence *>(bu_calloc(
	    description->occurrenceCount,
	    sizeof(*collector->manifest.occurrences),
	    "bobol draw manifest occurrences"));
    return collector->manifest.occurrences ? 1 : 0;
}

static int
bobol_draw_manifest_collect_occurrence(
	const BObolDrawManifestOccurrence *occurrence, size_t index,
	void *userData)
{
    BObolDrawManifestCollector *collector =
	static_cast<BObolDrawManifestCollector *>(userData);
    if (!collector || !occurrence ||
	index >= collector->manifest.occurrenceCount)
	return 0;
    BObolDrawManifestOccurrence &copy =
	collector->manifest.occurrences[index];
    copy = *occurrence;
    copy.path = occurrence->path ? bu_strdup(occurrence->path) : NULL;
    copy.sourceName = occurrence->sourceName ?
	bu_strdup(occurrence->sourceName) : NULL;
    copy.meshAssetPath = occurrence->meshAssetPath ?
	bu_strdup(occurrence->meshAssetPath) : NULL;
    copy.meshAssetName = occurrence->meshAssetName ?
	bu_strdup(occurrence->meshAssetName) : NULL;
    return copy.path && copy.sourceName &&
	(!occurrence->meshAssetPath || copy.meshAssetPath) &&
	(!occurrence->meshAssetName || copy.meshAssetName);
}

extern "C" int
bobol_draw_manifest_cache_get(db_i *dbip, const char *rootPath,
	BObolDrawManifest *manifest)
{
    if (!manifest)
	return BRLCAD_ERROR;
    BObolDrawManifestCollector collector;
    bobol_draw_manifest_init(&collector.manifest);
    const int status = bobol_draw_manifest_cache_stream(dbip, rootPath,
	bobol_draw_manifest_collect_begin,
	bobol_draw_manifest_collect_occurrence, &collector);
    if (status != BRLCAD_OK) {
	bobol_draw_manifest_free(&collector.manifest);
	return BRLCAD_ERROR;
    }
    bobol_draw_manifest_free(manifest);
    *manifest = collector.manifest;
    return BRLCAD_OK;
}
