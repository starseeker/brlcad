/*                  D R A W _ C A C H E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file draw_cache.cpp */

#include "common.h"

#include "brlobol/draw_cache.h"
#include "brlobol/lod_realization.h"

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
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <cmath>
#include <map>
#include <mutex>
#include <string>

#define BRLOBOL_DRAW_CACHE_FORMAT_FILE "draw_data.format"
#define BRLOBOL_DRAW_CACHE_CURRENT_FORMAT 3
#define BRLOBOL_DRAW_CACHE_AABB "bb"
#define BRLOBOL_DRAW_CACHE_OBB "obb"
#define BRLOBOL_DRAW_CACHE_METADATA "meta"
#define BRLOBOL_DRAW_CACHE_PATH_METADATA "meta_path"
#define BRLOBOL_DRAW_CACHE_PATH_METADATA_INDEX "meta_path_index"
#define BRLOBOL_DRAW_METADATA_DISK_MAGIC 0x4f424d45u /* OBME */

struct BRLObolDrawCacheContext {
    bu_cache *cache;
    char *registryKey;
};

struct BRLObolDrawCacheHandle {
    BRLObolDrawCacheContext *context;
    bu_cache *cache;
};

struct BRLObolDrawMetadataDiskRecord {
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
    char materialName[BRLOBOL_DRAW_CACHE_METADATA_NAME_MAX];
    int32_t hasShader;
    char shader[BRLOBOL_DRAW_CACHE_METADATA_SHADER_MAX];
};

static int
brlobol_draw_cache_semaphore(void)
{
    static int sem = 0;
    if (!sem)
	sem = bu_semaphore_register("BRLOBOL_DRAW_CACHE");
    return sem;
}

static std::mutex &
brlobol_draw_cache_registry_mutex(void)
{
    static std::mutex m;
    return m;
}

static std::map<std::string, BRLObolDrawCacheContext *> &
brlobol_draw_cache_registry(void)
{
    static std::map<std::string, BRLObolDrawCacheContext *> registry;
    return registry;
}

static void
brlobol_draw_cache_context_close(BRLObolDrawCacheContext *context)
{
    if (!context)
	return;
    if (context->cache)
	bu_cache_close(context->cache);
    if (context->registryKey)
	bu_free(context->registryKey, "brlobol draw cache registry key");
    BU_PUT(context, BRLObolDrawCacheContext);
}

static void
brlobol_draw_cache_registry_close_all(void)
{
    std::map<std::string, BRLObolDrawCacheContext *> closeMap;

    {
	std::lock_guard<std::mutex> guard(brlobol_draw_cache_registry_mutex());
	closeMap.swap(brlobol_draw_cache_registry());
    }

    for (std::map<std::string, BRLObolDrawCacheContext *>::iterator it =
	     closeMap.begin(); it != closeMap.end(); ++it)
	brlobol_draw_cache_context_close(it->second);
}

static const char *
brlobol_draw_proxy_component(int kind)
{
    switch (kind) {
	case BRLOBOL_LOD_PROXY_AABB:
	    return BRLOBOL_DRAW_CACHE_AABB;
	case BRLOBOL_LOD_PROXY_OBB:
	    return BRLOBOL_DRAW_CACHE_OBB;
	default:
	    return NULL;
    }
}

static size_t
brlobol_draw_proxy_point_count(int kind)
{
    switch (kind) {
	case BRLOBOL_LOD_PROXY_AABB:
	    return 2;
	case BRLOBOL_LOD_PROXY_OBB:
	    return 8;
	default:
	    return 0;
    }
}

static union tree *
    brlobol_draw_proxy_make_nop_tree(void)
{
    union tree *curtree = TREE_NULL;

    BU_GET(curtree, union tree);
    RT_TREE_INIT(curtree);
    curtree->tr_op = OP_NOP;
    return curtree;
}


static int
brlobol_draw_proxy_point_finite(const point_t p)
{
    return std::isfinite(p[X]) && std::isfinite(p[Y]) &&
	   std::isfinite(p[Z]);
}


static int
brlobol_draw_proxy_bbox_valid(const point_t bmin,
			      const point_t bmax)
{
    return brlobol_draw_proxy_point_finite(bmin) &&
	   brlobol_draw_proxy_point_finite(bmax) &&
	   bmin[X] <= bmax[X] && bmin[Y] <= bmax[Y] &&
	   bmin[Z] <= bmax[Z];
}


struct brlobol_draw_proxy_aabb_walk_ctx {
    struct db_i *dbip;
    point_t bmin;
    point_t bmax;
    int haveBounds;
};


static union tree *
    brlobol_draw_proxy_aabb_leaf_cb(struct db_tree_state *tsp,
				    const struct db_full_path *UNUSED(pathp),
				    struct directory *dp,
				    void *client_data)
{
    struct brlobol_draw_proxy_aabb_walk_ctx *ctx =
	(struct brlobol_draw_proxy_aabb_walk_ctx *)client_data;
    if (!ctx || !ctx->dbip || !tsp || !dp)
	return brlobol_draw_proxy_make_nop_tree();

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, ctx->dbip, tsp->ts_mat) < 0)
	return brlobol_draw_proxy_make_nop_tree();

    point_t bmin, bmax;
    const bn_tol tol = BN_TOL_INIT_TOL;
    VSETALL(bmin, INFINITY);
    VSETALL(bmax, -INFINITY);
    if (intern.idb_meth && intern.idb_meth->ft_bbox &&
	intern.idb_meth->ft_bbox(&intern, &bmin, &bmax, &tol) == 0 &&
	brlobol_draw_proxy_bbox_valid(bmin, bmax)) {
	VMINMAX(ctx->bmin, ctx->bmax, bmin);
	VMINMAX(ctx->bmin, ctx->bmax, bmax);
	ctx->haveBounds = 1;
    }

    rt_db_free_internal(&intern);
    return brlobol_draw_proxy_make_nop_tree();
}


static int
brlobol_draw_proxy_generate_fast_aabb(db_i *dbip,
				      const char *name,
				      point_t bmin,
				      point_t bmax)
{
    if (!dbip || !name || !name[0])
	return 0;

    struct brlobol_draw_proxy_aabb_walk_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.dbip = dbip;
    VSETALL(ctx.bmin, INFINITY);
    VSETALL(ctx.bmax, -INFINITY);

    struct db_tree_state init_state;
    db_init_db_tree_state(&init_state, dbip);
    init_state.ts_stop_at_regions = 0;

    const char *av[1] = {name};
    int ret = db_walk_tree_leaf_instances(dbip, 1, av, 1, &init_state,
					  NULL, NULL,
					  brlobol_draw_proxy_aabb_leaf_cb,
					  (void *)&ctx);
    db_free_db_tree_state(&init_state);
    if (ret < 0 || !ctx.haveBounds ||
	!brlobol_draw_proxy_bbox_valid(ctx.bmin, ctx.bmax))
	return 0;

    VMOVE(bmin, ctx.bmin);
    VMOVE(bmax, ctx.bmax);
    return 1;
}


static void
brlobol_draw_cache_key(char *key, const char *name, const char *component)
{
    bu_vls keystr = BU_VLS_INIT_ZERO;

    if (!key)
	return;
    key[0] = '\0';
    if (!name || !component)
	return;

    bu_vls_sprintf(&keystr, "%s", name);
    if (bu_vls_strlen(&keystr) < 10)
	bu_vls_printf(&keystr, "GGGGGGGGGGGGG");
    unsigned long long hash = bu_data_hash(bu_vls_cstr(&keystr),
					   bu_vls_strlen(&keystr) * sizeof(char));
    snprintf(key, BU_CACHE_KEY_MAXLEN, "%llu:%s", hash, component);
    bu_vls_free(&keystr);
}

static void
brlobol_draw_proxy_cache_key(char *key, const char *name, int kind)
{
    brlobol_draw_cache_key(key, name, brlobol_draw_proxy_component(kind));
}

static void
brlobol_draw_path_metadata_index_key(char *key, const char *path)
{
    brlobol_draw_cache_key(key, path, BRLOBOL_DRAW_CACHE_PATH_METADATA_INDEX);
}

static int
brlobol_draw_cache_key_component_is(const char *key, const char *component)
{
    if (!key || !component)
	return 0;

    const char *sep = strrchr(key, ':');
    return (sep && BU_STR_EQUAL(sep + 1, component)) ? 1 : 0;
}

static int
brlobol_draw_path_component_is_instance_of(const char *component,
					   size_t componentLen,
					   const char *name,
					   size_t nameLen)
{
    if (!component || !componentLen || !name || !nameLen)
	return 0;
    if (componentLen == nameLen &&
	strncmp(component, name, nameLen) == 0)
	return 1;
    if (componentLen <= nameLen + 1 ||
	strncmp(component, name, nameLen) != 0 ||
	component[nameLen] != '@')
	return 0;

    for (size_t i = nameLen + 1; i < componentLen; i++) {
	if (component[i] < '0' || component[i] > '9')
	    return 0;
    }
    return 1;
}

static int
brlobol_draw_path_contains_object_name(const char *path, const char *name)
{
    if (!path || !path[0] || !name || !name[0])
	return 0;

    size_t nameLen = strlen(name);
    const char *start = path;
    while (*start) {
	const char *end = strchr(start, '/');
	size_t componentLen = end ? (size_t)(end - start) : strlen(start);
	if (brlobol_draw_path_component_is_instance_of(start, componentLen,
		name, nameLen))
	    return 1;
	if (!end)
	    break;
	start = end + 1;
    }

    return 0;
}

static int
brlobol_draw_cache_ensure_root(void)
{
    char dir[MAXPATHLEN];

    bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, NULL);
    if (!bu_file_exists(dir, NULL))
	bu_mkdir(dir);
    if (!bu_file_exists(dir, NULL))
	return 0;

    bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, BRLOBOL_DRAW_CACHE_DIR, NULL);
    if (!bu_file_exists(dir, NULL))
	bu_mkdir(dir);
    if (!bu_file_exists(dir, NULL))
	return 0;

    char formatPath[MAXPATHLEN];
    bu_dir(formatPath, MAXPATHLEN, BU_DIR_CACHE, BRLOBOL_DRAW_CACHE_DIR,
	   BRLOBOL_DRAW_CACHE_FORMAT_FILE, NULL);
    long diskFormatVersion = -1;
    FILE *fp = fopen(formatPath, "r");
    if (fp) {
	(void)fscanf(fp, "%ld", &diskFormatVersion);
	fclose(fp);
    }
    if (diskFormatVersion > 0 &&
	diskFormatVersion != BRLOBOL_DRAW_CACHE_CURRENT_FORMAT) {
	brlobol_draw_cache_registry_close_all();
	bu_cache_erase(BRLOBOL_DRAW_CACHE_DIR);
	bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, BRLOBOL_DRAW_CACHE_DIR, NULL);
	if (!bu_file_exists(dir, NULL))
	    bu_mkdir(dir);
	if (!bu_file_exists(dir, NULL))
	    return 0;
    }
    fp = fopen(formatPath, "w");
    if (fp) {
	fprintf(fp, "%d\n", BRLOBOL_DRAW_CACHE_CURRENT_FORMAT);
	fclose(fp);
    }

    return 1;
}

static int
brlobol_draw_cache_db_name(bu_vls *fname, db_i *dbip)
{
    if (!fname || !dbip)
	return 0;

    const char *ctxName = dbip->dbi_filename;
    char inmemCtxName[128] = {0};
    if (!ctxName || !ctxName[0]) {
	snprintf(inmemCtxName, sizeof(inmemCtxName),
		 "brlobol_inmem_draw_cache_%p", (void *)dbip);
	ctxName = inmemCtxName;
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

static int
brlobol_draw_cache_open(BRLObolDrawCacheHandle *handle, db_i *dbip)
{
    bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls cpath = BU_VLS_INIT_ZERO;
    char cacheRoot[MAXPATHLEN] = {0};

    if (!handle || !dbip)
	return 0;
    handle->context = NULL;
    handle->cache = NULL;

    if (!brlobol_draw_cache_ensure_root())
	return 0;
    if (!brlobol_draw_cache_db_name(&fname, dbip))
	return 0;

    bu_vls_sprintf(&cpath, "%s/%s_draw", BRLOBOL_DRAW_CACHE_DIR,
		   bu_vls_cstr(&fname));
    bu_dir(cacheRoot, MAXPATHLEN, BU_DIR_CACHE, NULL);
    std::string registryKey(cacheRoot);
    registryKey.append("/");
    registryKey.append(bu_vls_cstr(&cpath));

    {
	std::lock_guard<std::mutex> guard(brlobol_draw_cache_registry_mutex());
	std::map<std::string, BRLObolDrawCacheContext *>::iterator it =
	    brlobol_draw_cache_registry().find(registryKey);
	if (it != brlobol_draw_cache_registry().end()) {
	    handle->context = it->second;
	    handle->cache = it->second->cache;
	    bu_vls_free(&cpath);
	    bu_vls_free(&fname);
	    return handle->cache ? 1 : 0;
	}

	BRLObolDrawCacheContext *context;
	BU_GET(context, BRLObolDrawCacheContext);
	context->cache = bu_cache_open(bu_vls_cstr(&cpath), 1, 0);
	context->registryKey = bu_strdup(registryKey.c_str());
	if (!context->cache) {
	    brlobol_draw_cache_context_close(context);
	    bu_vls_free(&cpath);
	    bu_vls_free(&fname);
	    return 0;
	}

	brlobol_draw_cache_registry()[registryKey] = context;
	handle->context = context;
	handle->cache = context->cache;
    }

    bu_vls_free(&cpath);
    bu_vls_free(&fname);
    return handle->cache ? 1 : 0;
}

static void
brlobol_draw_cache_close(BRLObolDrawCacheHandle *handle)
{
    if (!handle)
	return;
    handle->context = NULL;
    handle->cache = NULL;
}

extern "C" void
brlobol_draw_cache_status_init(BRLObolDrawCacheStatus *status)
{
    if (!status)
	return;
    memset(status, 0, sizeof(*status));
}

extern "C" void
brlobol_draw_proxy_record_init(BRLObolDrawProxyRecord *record)
{
    if (!record)
	return;
    memset(record, 0, sizeof(*record));
}

extern "C" void
brlobol_draw_metadata_record_init(BRLObolDrawMetadataRecord *record)
{
    if (!record)
	return;
    memset(record, 0, sizeof(*record));
}

extern "C" void
brlobol_draw_cache_clear_all(void)
{
    char dir[MAXPATHLEN];
    int sem = brlobol_draw_cache_semaphore();

    bu_semaphore_acquire(sem);
    brlobol_draw_cache_registry_close_all();
    bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, BRLOBOL_DRAW_CACHE_DIR, NULL);
    bu_dirclear(dir);
    bu_semaphore_release(sem);
}

extern "C" int
brlobol_draw_cache_clear_database(db_i *dbip)
{
    BRLObolDrawCacheHandle handle;
    char **keys = NULL;
    int nkeys = 0;

    if (!dbip)
	return BRLCAD_ERROR;

    int sem = brlobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    int opened = brlobol_draw_cache_open(&handle, dbip);
    if (opened) {
	nkeys = bu_cache_keys(&keys, handle.cache);
	for (int i = 0; i < nkeys; i++)
	    bu_cache_clear(keys[i], handle.cache, NULL);
    }
    brlobol_draw_cache_close(&handle);
    bu_semaphore_release(sem);

    if (nkeys)
	bu_argv_free((size_t)nkeys, keys);
    return opened ? BRLCAD_OK : BRLCAD_ERROR;
}

static void
brlobol_draw_proxy_status_current(db_i *dbip,
				  const char *name,
				  int kind,
				  BRLObolDrawCacheStatus *status)
{
    void *data = NULL;
    char key[BU_CACHE_KEY_MAXLEN] = {0};
    size_t expectedSize = brlobol_draw_proxy_point_count(kind) *
			  sizeof(point_t);

    if (!status)
	return;
    brlobol_draw_cache_status_init(status);
    if (!dbip || !name || !expectedSize)
	return;

    directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    status->directoryFound = (dp != RT_DIR_NULL) ? 1 : 0;

    brlobol_draw_proxy_cache_key(key, name, kind);
    if (!key[0])
	return;

    int sem = brlobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BRLObolDrawCacheHandle handle;
    if (brlobol_draw_cache_open(&handle, dbip)) {
	size_t dsize = bu_cache_get(&data, key, handle.cache, NULL);
	status->hasCachedPayload = (data && dsize == expectedSize) ? 1 : 0;
	brlobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    if (data)
	bu_free(data, "brlobol draw proxy status data");
}

extern "C" int
brlobol_draw_proxy_cache_status(db_i *dbip,
				const char *name,
				int kind,
				BRLObolDrawCacheStatus *status)
{
    if (!status || !dbip || !name || !brlobol_draw_proxy_point_count(kind))
	return BRLCAD_ERROR;
    brlobol_draw_proxy_status_current(dbip, name, kind, status);
    return BRLCAD_OK;
}

extern "C" int
brlobol_draw_proxy_cache_invalidate(db_i *dbip,
				    const char *name,
				    int kind,
				    BRLObolDrawCacheStatus *status)
{
    BRLObolDrawCacheStatus current;
    char key[BU_CACHE_KEY_MAXLEN] = {0};

    if (status)
	brlobol_draw_cache_status_init(status);
    if (!dbip || !name || !brlobol_draw_proxy_point_count(kind))
	return BRLCAD_ERROR;

    brlobol_draw_proxy_status_current(dbip, name, kind, &current);
    brlobol_draw_proxy_cache_key(key, name, kind);
    if (!key[0])
	return BRLCAD_ERROR;

    int sem = brlobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BRLObolDrawCacheHandle handle;
    if (brlobol_draw_cache_open(&handle, dbip)) {
	bu_cache_clear(key, handle.cache, NULL);
	brlobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    current.clearedCacheEntry = 1;
    current.hasCachedPayload = 0;
    if (status)
	*status = current;
    return BRLCAD_OK;
}

extern "C" int
brlobol_draw_proxy_cache_store(db_i *dbip,
			       const char *name,
			       int kind,
			       const point_t *points,
			       size_t pointCount,
			       BRLObolDrawCacheStatus *status)
{
    BRLObolDrawCacheStatus current;
    char key[BU_CACHE_KEY_MAXLEN] = {0};
    size_t expectedCount = brlobol_draw_proxy_point_count(kind);

    if (status)
	brlobol_draw_cache_status_init(status);
    if (!dbip || !name || !points || !expectedCount ||
	pointCount != expectedCount)
	return BRLCAD_ERROR;
    brlobol_draw_cache_status_init(&current);

    directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    current.directoryFound = (dp != RT_DIR_NULL) ? 1 : 0;
    if (dp == RT_DIR_NULL) {
	if (status)
	    *status = current;
	return BRLCAD_ERROR;
    }

    brlobol_draw_proxy_cache_key(key, name, kind);
    if (!key[0])
	return BRLCAD_ERROR;

    int sem = brlobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BRLObolDrawCacheHandle handle;
    size_t wsize = 0;
    if (brlobol_draw_cache_open(&handle, dbip)) {
	wsize = bu_cache_write((void *)points,
			       expectedCount * sizeof(point_t), key, handle.cache, NULL);
	brlobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    current.hasCachedPayload = (wsize == expectedCount *
				sizeof(point_t)) ? 1 : 0;
    current.generatedCacheEntry = current.hasCachedPayload ? 1 : 0;
    if (status)
	*status = current;
    return current.hasCachedPayload ? BRLCAD_OK : BRLCAD_ERROR;
}

extern "C" int
brlobol_draw_proxy_cache_get(db_i *dbip,
			     const char *name,
			     int kind,
			     BRLObolDrawProxyRecord *record)
{
    void *data = NULL;
    char key[BU_CACHE_KEY_MAXLEN] = {0};
    size_t expectedCount = brlobol_draw_proxy_point_count(kind);
    size_t expectedSize = expectedCount * sizeof(point_t);

    if (!dbip || !name || !record || !expectedCount)
	return BRLCAD_ERROR;
    brlobol_draw_proxy_record_init(record);

    brlobol_draw_proxy_cache_key(key, name, kind);
    if (!key[0])
	return BRLCAD_ERROR;

    int sem = brlobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BRLObolDrawCacheHandle handle;
    size_t dsize = 0;
    if (brlobol_draw_cache_open(&handle, dbip)) {
	dsize = bu_cache_get(&data, key, handle.cache, NULL);
	brlobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    if (!data || dsize != expectedSize) {
	if (data)
	    bu_free(data, "brlobol draw proxy get data");
	return BRLCAD_ERROR;
    }

    record->kind = kind;
    record->pointCount = expectedCount;
    memcpy(record->points, data, expectedSize);
    bu_free(data, "brlobol draw proxy get data");
    return BRLCAD_OK;
}

extern "C" int
brlobol_draw_proxy_cache_refresh(db_i *dbip,
				 const char *name,
				 int kind,
				 BRLObolDrawCacheStatus *status)
{
    point_t points[8];
    rt_db_internal intern;
    const bn_tol tol = BN_TOL_INIT_TOL;

    if (status)
	brlobol_draw_cache_status_init(status);
    if (!dbip || !name || !brlobol_draw_proxy_point_count(kind))
	return BRLCAD_ERROR;

    directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
	if (status)
	    status->directoryFound = 0;
	return BRLCAD_ERROR;
    }

    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
	return BRLCAD_ERROR;

    int ret = BRLCAD_ERROR;
    if (kind == BRLOBOL_LOD_PROXY_AABB) {
	point_t bmin, bmax;
	VSETALL(bmin, INFINITY);
	VSETALL(bmax, -INFINITY);
	if (intern.idb_meth && intern.idb_meth->ft_bbox &&
	    intern.idb_meth->ft_bbox(&intern, &bmin, &bmax,
				     &tol) == 0 &&
	    brlobol_draw_proxy_bbox_valid(bmin, bmax)) {
	    VMOVE(points[0], bmin);
	    VMOVE(points[1], bmax);
	    ret = brlobol_draw_proxy_cache_store(dbip, name, kind, points,
						 2, status);
	} else if (brlobol_draw_proxy_generate_fast_aabb(dbip, name,
		   bmin, bmax)) {
	    VMOVE(points[0], bmin);
	    VMOVE(points[1], bmax);
	    ret = brlobol_draw_proxy_cache_store(dbip, name, kind, points,
						 2, status);
	}
    } else if (kind == BRLOBOL_LOD_PROXY_OBB) {
	rt_arb_internal arb;
	arb.magic = RT_ARB_INTERNAL_MAGIC;
	if (intern.idb_meth && intern.idb_meth->ft_oriented_bbox &&
	    intern.idb_meth->ft_oriented_bbox(&arb, &intern,
					      BN_TOL_DIST) == 0) {
	    for (size_t i = 0; i < 8; i++)
		VMOVE(points[i], arb.pt[i]);
	    ret = brlobol_draw_proxy_cache_store(dbip, name, kind, points,
						 8, status);
	}
    }

    rt_db_free_internal(&intern);
    if (ret != BRLCAD_OK)
	(void)brlobol_draw_proxy_cache_invalidate(dbip, name, kind, status);
    return ret;
}

static int
brlobol_draw_metadata_int_attr(const bu_attribute_value_set *avs,
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
brlobol_draw_metadata_color_attr(const bu_attribute_value_set *avs,
				 BRLObolDrawMetadataRecord *record)
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
brlobol_draw_metadata_string_attr(const bu_attribute_value_set *avs,
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
brlobol_draw_metadata_to_disk(BRLObolDrawMetadataDiskRecord *disk,
			      const BRLObolDrawMetadataRecord *record)
{
    memset(disk, 0, sizeof(*disk));
    disk->magic = BRLOBOL_DRAW_METADATA_DISK_MAGIC;
    disk->version = BRLOBOL_DRAW_CACHE_CURRENT_FORMAT;
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
brlobol_draw_metadata_from_disk(BRLObolDrawMetadataRecord *record,
				const BRLObolDrawMetadataDiskRecord *disk)
{
    brlobol_draw_metadata_record_init(record);
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
brlobol_draw_metadata_disk_valid(const void *data, size_t dsize)
{
    const BRLObolDrawMetadataDiskRecord *disk =
	(const BRLObolDrawMetadataDiskRecord *)data;
    return data && dsize == sizeof(*disk) &&
	   disk->magic == BRLOBOL_DRAW_METADATA_DISK_MAGIC &&
	   disk->version == BRLOBOL_DRAW_CACHE_CURRENT_FORMAT;
}

static int
brlobol_draw_metadata_target_found(db_i *dbip, const char *name, int pathAware)
{
    if (!dbip || !name || !name[0])
	return 0;

    if (!pathAware)
	return db_lookup(dbip, name, LOOKUP_QUIET) != RT_DIR_NULL;

    struct db_full_path path;
    db_full_path_init(&path);
    int found = (db_string_to_path(&path, dbip, name) == 0 &&
		 DB_FULL_PATH_CUR_DIR(&path)) ? 1 : 0;
    db_free_full_path(&path);
    return found;
}

static void
brlobol_draw_metadata_record_apply_directory(
    BRLObolDrawMetadataRecord *record,
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
brlobol_draw_metadata_color_channel(float value)
{
    long color = lround(value * 255.0f);
    if (color < 0)
	color = 0;
    if (color > 255)
	color = 255;
    return (int)color;
}

static void
brlobol_draw_metadata_record_apply_tree_state(
    BRLObolDrawMetadataRecord *record,
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

    if (state->ts_mater.ma_color_valid) {
	record->hasColor = 1;
	record->color[0] =
	    (unsigned char)brlobol_draw_metadata_color_channel(
		state->ts_mater.ma_color[0]);
	record->color[1] =
	    (unsigned char)brlobol_draw_metadata_color_channel(
		state->ts_mater.ma_color[1]);
	record->color[2] =
	    (unsigned char)brlobol_draw_metadata_color_channel(
		state->ts_mater.ma_color[2]);
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

static void
brlobol_draw_metadata_status_current_for_component(
    db_i *dbip,
    const char *name,
    const char *component,
    int pathAware,
    BRLObolDrawCacheStatus *status)
{
    void *data = NULL;
    char key[BU_CACHE_KEY_MAXLEN] = {0};

    if (!status)
	return;
    brlobol_draw_cache_status_init(status);
    if (!dbip || !name)
	return;

    status->directoryFound =
	brlobol_draw_metadata_target_found(dbip, name, pathAware);

    brlobol_draw_cache_key(key, name, component);
    if (!key[0])
	return;

    int sem = brlobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BRLObolDrawCacheHandle handle;
    if (brlobol_draw_cache_open(&handle, dbip)) {
	size_t dsize = bu_cache_get(&data, key, handle.cache, NULL);
	status->hasCachedPayload =
	    brlobol_draw_metadata_disk_valid(data, dsize) ? 1 : 0;
	brlobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    if (data)
	bu_free(data, "brlobol draw metadata status data");
}

static int
brlobol_draw_metadata_cache_invalidate_for_component(
    db_i *dbip,
    const char *name,
    const char *component,
    int pathAware,
    BRLObolDrawCacheStatus *status)
{
    BRLObolDrawCacheStatus current;
    char key[BU_CACHE_KEY_MAXLEN] = {0};
    char pathIndexKey[BU_CACHE_KEY_MAXLEN] = {0};

    if (status)
	brlobol_draw_cache_status_init(status);
    if (!dbip || !name || !component)
	return BRLCAD_ERROR;
    brlobol_draw_metadata_status_current_for_component(dbip, name, component,
	    pathAware, &current);
    brlobol_draw_cache_key(key, name, component);
    if (!key[0])
	return BRLCAD_ERROR;
    if (BU_STR_EQUAL(component, BRLOBOL_DRAW_CACHE_PATH_METADATA))
	brlobol_draw_path_metadata_index_key(pathIndexKey, name);

    int sem = brlobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BRLObolDrawCacheHandle handle;
    if (brlobol_draw_cache_open(&handle, dbip)) {
	bu_cache_clear(key, handle.cache, NULL);
	if (pathIndexKey[0])
	    bu_cache_clear(pathIndexKey, handle.cache, NULL);
	brlobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    current.clearedCacheEntry = 1;
    current.hasCachedPayload = 0;
    if (status)
	*status = current;
    return BRLCAD_OK;
}

static int
brlobol_draw_metadata_cache_store_for_component(
    db_i *dbip,
    const char *name,
    const char *component,
    int pathAware,
    const BRLObolDrawMetadataRecord *record,
    BRLObolDrawCacheStatus *status)
{
    BRLObolDrawCacheStatus current;
    BRLObolDrawMetadataDiskRecord disk;
    char key[BU_CACHE_KEY_MAXLEN] = {0};
    char pathIndexKey[BU_CACHE_KEY_MAXLEN] = {0};

    if (status)
	brlobol_draw_cache_status_init(status);
    if (!dbip || !name || !component || !record)
	return BRLCAD_ERROR;
    brlobol_draw_cache_status_init(&current);

    current.directoryFound =
	brlobol_draw_metadata_target_found(dbip, name, pathAware);
    if (!current.directoryFound) {
	if (status)
	    *status = current;
	return BRLCAD_ERROR;
    }

    brlobol_draw_cache_key(key, name, component);
    if (!key[0])
	return BRLCAD_ERROR;
    if (BU_STR_EQUAL(component, BRLOBOL_DRAW_CACHE_PATH_METADATA))
	brlobol_draw_path_metadata_index_key(pathIndexKey, name);
    brlobol_draw_metadata_to_disk(&disk, record);

    int sem = brlobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BRLObolDrawCacheHandle handle;
    size_t wsize = 0;
    size_t iwsize = 0;
    if (brlobol_draw_cache_open(&handle, dbip)) {
	wsize = bu_cache_write((void *)&disk, sizeof(disk), key,
			       handle.cache, NULL);
	if (pathIndexKey[0] && wsize == sizeof(disk))
	    iwsize = bu_cache_write((void *)name, strlen(name) + 1,
				    pathIndexKey, handle.cache, NULL);
	brlobol_draw_cache_close(&handle);
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
brlobol_draw_metadata_cache_get_for_component(
    db_i *dbip,
    const char *name,
    const char *component,
    BRLObolDrawMetadataRecord *record)
{
    void *data = NULL;
    char key[BU_CACHE_KEY_MAXLEN] = {0};

    if (!dbip || !name || !component || !record)
	return BRLCAD_ERROR;
    brlobol_draw_metadata_record_init(record);

    brlobol_draw_cache_key(key, name, component);
    if (!key[0])
	return BRLCAD_ERROR;

    int sem = brlobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BRLObolDrawCacheHandle handle;
    size_t dsize = 0;
    if (brlobol_draw_cache_open(&handle, dbip)) {
	dsize = bu_cache_get(&data, key, handle.cache, NULL);
	brlobol_draw_cache_close(&handle);
    }
    bu_semaphore_release(sem);

    if (!brlobol_draw_metadata_disk_valid(data, dsize)) {
	if (data)
	    bu_free(data, "brlobol draw metadata get data");
	return BRLCAD_ERROR;
    }

    brlobol_draw_metadata_from_disk(record,
				    (const BRLObolDrawMetadataDiskRecord *)data);
    bu_free(data, "brlobol draw metadata get data");
    return BRLCAD_OK;
}

extern "C" int
brlobol_draw_metadata_cache_status(db_i *dbip,
				   const char *name,
				   BRLObolDrawCacheStatus *status)
{
    if (!status || !dbip || !name)
	return BRLCAD_ERROR;
    brlobol_draw_metadata_status_current_for_component(dbip, name,
	    BRLOBOL_DRAW_CACHE_METADATA, 0, status);
    return BRLCAD_OK;
}

extern "C" int
brlobol_draw_metadata_cache_invalidate(db_i *dbip,
				       const char *name,
				       BRLObolDrawCacheStatus *status)
{
    return brlobol_draw_metadata_cache_invalidate_for_component(
	dbip, name, BRLOBOL_DRAW_CACHE_METADATA, 0, status);
}

extern "C" int
brlobol_draw_metadata_cache_store(db_i *dbip,
				  const char *name,
				  const BRLObolDrawMetadataRecord *record,
				  BRLObolDrawCacheStatus *status)
{
    return brlobol_draw_metadata_cache_store_for_component(
	dbip, name, BRLOBOL_DRAW_CACHE_METADATA, 0, record, status);
}

extern "C" int
brlobol_draw_metadata_cache_get(db_i *dbip,
				const char *name,
				BRLObolDrawMetadataRecord *record)
{
    return brlobol_draw_metadata_cache_get_for_component(
	dbip, name, BRLOBOL_DRAW_CACHE_METADATA, record);
}

extern "C" int
brlobol_draw_metadata_cache_refresh(db_i *dbip,
				    const char *name,
				    BRLObolDrawCacheStatus *status)
{
    BRLObolDrawMetadataRecord record;
    bu_attribute_value_set avs = BU_AVS_INIT_ZERO;
    int haveAttrs = 0;
    int attrInt = 0;

    if (status)
	brlobol_draw_cache_status_init(status);
    if (!dbip || !name)
	return BRLCAD_ERROR;

    directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
	if (status)
	    status->directoryFound = 0;
	return BRLCAD_ERROR;
    }

    brlobol_draw_metadata_record_init(&record);
    brlobol_draw_metadata_record_apply_directory(&record, dp);

    if (db5_get_attributes(dbip, &avs, dp) == 0)
	haveAttrs = 1;

    if (haveAttrs) {
	const char *region = bu_avs_get(&avs, "region");
	if (region && bu_str_true(region))
	    record.isRegion = 1;

	if (brlobol_draw_metadata_int_attr(&avs, "region_id", &attrInt)) {
	    record.hasRegionId = 1;
	    record.regionId = attrInt;
	} else if (record.isRegion) {
	    record.hasRegionId = 1;
	    record.regionId = 0;
	}
	if (brlobol_draw_metadata_int_attr(&avs, "air", &attrInt)) {
	    record.hasAircode = 1;
	    record.aircode = attrInt;
	}
	if (brlobol_draw_metadata_int_attr(&avs, "los", &attrInt)) {
	    record.hasLos = 1;
	    record.los = attrInt;
	}
	if (brlobol_draw_metadata_int_attr(&avs, "material_id", &attrInt)) {
	    record.hasMaterialId = 1;
	    record.materialId = attrInt;
	}

	const char *inherit = bu_avs_get(&avs, "inherit");
	if (inherit && inherit[0]) {
	    record.hasInherit = 1;
	    record.inherit = bu_str_true(inherit) ? 1 : 0;
	}

	(void)brlobol_draw_metadata_color_attr(&avs, &record);
	brlobol_draw_metadata_string_attr(&avs, "material_name",
					  record.materialName, sizeof(record.materialName),
					  &record.hasMaterialName);
	brlobol_draw_metadata_string_attr(&avs, "shader",
					  record.shader, sizeof(record.shader), &record.hasShader);
    } else if (record.isRegion) {
	record.hasRegionId = 1;
	record.regionId = 0;
    }

    if (haveAttrs)
	bu_avs_free(&avs);

    return brlobol_draw_metadata_cache_store(dbip, name, &record, status);
}

extern "C" int
brlobol_draw_path_metadata_cache_status(db_i *dbip,
					const char *path,
					BRLObolDrawCacheStatus *status)
{
    if (!status || !dbip || !path)
	return BRLCAD_ERROR;
    brlobol_draw_metadata_status_current_for_component(dbip, path,
	    BRLOBOL_DRAW_CACHE_PATH_METADATA, 1, status);
    return BRLCAD_OK;
}

extern "C" int
brlobol_draw_path_metadata_cache_invalidate(db_i *dbip,
					    const char *path,
					    BRLObolDrawCacheStatus *status)
{
    return brlobol_draw_metadata_cache_invalidate_for_component(
	dbip, path, BRLOBOL_DRAW_CACHE_PATH_METADATA, 1, status);
}

extern "C" int
brlobol_draw_path_metadata_cache_invalidate_object(
    db_i *dbip,
    const char *name,
    BRLObolDrawCacheStatus *status)
{
    char **keys = NULL;
    int nkeys = 0;
    int cleared = 0;
    int opened = 0;

    if (status)
	brlobol_draw_cache_status_init(status);
    if (!dbip || !name || !name[0])
	return BRLCAD_ERROR;

    int sem = brlobol_draw_cache_semaphore();
    bu_semaphore_acquire(sem);
    BRLObolDrawCacheHandle handle;
    opened = brlobol_draw_cache_open(&handle, dbip);
    if (opened) {
	nkeys = bu_cache_keys(&keys, handle.cache);
	for (int i = 0; i < nkeys; i++) {
	    if (!brlobol_draw_cache_key_component_is(keys[i],
		    BRLOBOL_DRAW_CACHE_PATH_METADATA_INDEX))
		continue;

	    void *data = NULL;
	    size_t dsize = bu_cache_get(&data, keys[i], handle.cache, NULL);
	    const char *path = (const char *)data;
	    if (data && dsize > 0 && ((const char *)data)[dsize - 1] == '\0' &&
		brlobol_draw_path_contains_object_name(path, name)) {
		char pathMetadataKey[BU_CACHE_KEY_MAXLEN] = {0};
		brlobol_draw_cache_key(pathMetadataKey, path,
				       BRLOBOL_DRAW_CACHE_PATH_METADATA);
		if (pathMetadataKey[0])
		    bu_cache_clear(pathMetadataKey, handle.cache, NULL);
		bu_cache_clear(keys[i], handle.cache, NULL);
		cleared = 1;
	    }
	    if (data)
		bu_free(data, "brlobol draw path metadata index data");
	}
	brlobol_draw_cache_close(&handle);
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
brlobol_draw_path_metadata_cache_store(db_i *dbip,
				       const char *path,
				       const BRLObolDrawMetadataRecord *record,
				       BRLObolDrawCacheStatus *status)
{
    return brlobol_draw_metadata_cache_store_for_component(
	dbip, path, BRLOBOL_DRAW_CACHE_PATH_METADATA, 1, record, status);
}

extern "C" int
brlobol_draw_path_metadata_cache_get(db_i *dbip,
				     const char *path,
				     BRLObolDrawMetadataRecord *record)
{
    return brlobol_draw_metadata_cache_get_for_component(
	dbip, path, BRLOBOL_DRAW_CACHE_PATH_METADATA, record);
}

extern "C" int
brlobol_draw_path_metadata_cache_refresh(db_i *dbip,
					 const char *path,
					 BRLObolDrawCacheStatus *status)
{
    BRLObolDrawMetadataRecord record;
    struct db_tree_state state;
    struct db_full_path fullPath;

    if (status)
	brlobol_draw_cache_status_init(status);
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

    brlobol_draw_metadata_record_init(&record);
    brlobol_draw_metadata_record_apply_directory(&record,
	    DB_FULL_PATH_CUR_DIR(&fullPath));
    brlobol_draw_metadata_record_apply_tree_state(&record, &state);

    db_free_full_path(&fullPath);
    db_free_db_tree_state(&state);

    return brlobol_draw_path_metadata_cache_store(dbip, path, &record,
	    status);
}
