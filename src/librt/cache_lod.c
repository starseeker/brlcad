/*                     C A C H E _ L O D . C
 * BRL-CAD
 *
 * Copyright (c) 2016-2026 United States Government as represented by
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
/** @file cache.c
 *
 * Caching of LoD drawing data
 */

#include "common.h"

/* implementation headers */
#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "bu/app.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "bu/path.h"
#include "bu/process.h"
#include "bu/time.h"
#include "bsg/lod.h"
#include "rt/db_instance.h"
#include "rt/view.h"

#include "./librt_private.h"

struct rt_mesh_lod {
    struct bsg_mesh_lod *bsg;
};

struct rt_mesh_lod_context {
    struct bsg_mesh_lod_context *bsg;
};

struct rt_mesh_lod_bot_detail {
    struct db_i *dbip;
    struct directory *dp;
    struct rt_db_internal *intern;
};

struct rt_mesh_lod_detail_callbacks {
    rt_mesh_lod_detail_setup_callback setup_clbk;
    rt_mesh_lod_detail_clear_callback clear_clbk;
    rt_mesh_lod_detail_free_callback free_clbk;
    void *cb_data;
};

static void
rt_mesh_lod_bsg_detail_clear(struct bsg_mesh_lod *lod)
{
    if (!lod)
	return;

    lod->faces = NULL;
    lod->fcnt = 0;
    lod->pcnt = 0;
    lod->points = NULL;
    lod->porig_cnt = 0;
    lod->points_orig = NULL;
    lod->normals = NULL;
}

static int
rt_mesh_lod_size_to_int(size_t count, int *out)
{
    if (!out || count > (size_t)INT_MAX)
	return 0;

    *out = (int)count;
    return 1;
}

static int
rt_mesh_lod_detail_setup_bsg(struct bsg_mesh_lod *lod, void *cb_data)
{
    struct rt_mesh_lod_detail_callbacks *callbacks =
	(struct rt_mesh_lod_detail_callbacks *)cb_data;
    struct rt_mesh_lod_detail detail;
    int face_count = 0;
    int point_count = 0;
    int point_orig_count = 0;

    if (!lod || !callbacks || !callbacks->setup_clbk)
	return -1;

    rt_mesh_lod_detail_init(&detail);
    if (callbacks->setup_clbk(&detail, callbacks->cb_data) != 0)
	return -1;

    if (!rt_mesh_lod_size_to_int(detail.face_count, &face_count) ||
	    !rt_mesh_lod_size_to_int(detail.point_count, &point_count) ||
	    !rt_mesh_lod_size_to_int(detail.point_orig_count,
		&point_orig_count)) {
	if (callbacks->clear_clbk)
	    callbacks->clear_clbk(callbacks->cb_data);
	rt_mesh_lod_bsg_detail_clear(lod);
	return -1;
    }

    lod->faces = detail.faces;
    lod->fcnt = face_count;
    lod->points = detail.points;
    lod->pcnt = point_count;
    lod->points_orig = detail.points_orig;
    lod->porig_cnt = point_orig_count;
    lod->normals = detail.normals;

    return 0;
}

static int
rt_mesh_lod_detail_clear_bsg(struct bsg_mesh_lod *lod, void *cb_data)
{
    struct rt_mesh_lod_detail_callbacks *callbacks =
	(struct rt_mesh_lod_detail_callbacks *)cb_data;
    int ret = 0;

    if (callbacks && callbacks->clear_clbk)
	ret = callbacks->clear_clbk(callbacks->cb_data);

    rt_mesh_lod_bsg_detail_clear(lod);
    return ret;
}

static int
rt_mesh_lod_detail_free_bsg(struct bsg_mesh_lod *lod, void *cb_data)
{
    struct rt_mesh_lod_detail_callbacks *callbacks =
	(struct rt_mesh_lod_detail_callbacks *)cb_data;
    int ret = 0;

    if (callbacks) {
	if (callbacks->free_clbk)
	    ret = callbacks->free_clbk(callbacks->cb_data);
	else if (callbacks->clear_clbk)
	    ret = callbacks->clear_clbk(callbacks->cb_data);
	BU_PUT(callbacks, struct rt_mesh_lod_detail_callbacks);
    }

    rt_mesh_lod_bsg_detail_clear(lod);
    return ret;
}

static struct bsg_mesh_lod_context *
rt_mesh_lod_context_bsg(struct rt_mesh_lod_context *c)
{
    return c ? c->bsg : NULL;
}

static struct bsg_mesh_lod_context *
rt_mesh_lod_context_ensure(struct db_i *dbip)
{
    if (!dbip || !dbip->i)
	return NULL;

    if (!dbip->i->mesh_c) {
	const char *ctx_name = dbip->dbi_filename;
	char inmem_ctx_name[128] = {0};

	if (!ctx_name || !ctx_name[0]) {
	    snprintf(inmem_ctx_name, sizeof(inmem_ctx_name),
		    "rt_inmem_mesh_lod_%p", (void *)dbip);
	    ctx_name = inmem_ctx_name;
	}

	struct bsg_mesh_lod_context *bsg = bsg_mesh_lod_context_create(ctx_name);
	if (!bsg) {
	    bu_log("Error creating mesh LoD context for %s\n", ctx_name);
	    return NULL;
	}

	BU_GET(dbip->i->mesh_c, struct rt_mesh_lod_context);
	dbip->i->mesh_c->bsg = bsg;
    }

    return rt_mesh_lod_context_bsg(dbip->i->mesh_c);
}

static int
rt_mesh_lod_bot_detail_setup(struct bsg_mesh_lod *lod, void *cb_data)
{
    if (!lod || !cb_data)
	return -1;

    struct rt_mesh_lod_bot_detail *cd = (struct rt_mesh_lod_bot_detail *)cb_data;
    if (!cd->dbip || !cd->dp)
	return -1;

    if (cd->intern)
	return 0;

    BU_GET(cd->intern, struct rt_db_internal);
    RT_DB_INTERNAL_INIT(cd->intern);
    if (rt_db_get_internal(cd->intern, cd->dp, cd->dbip, NULL) < 0) {
	BU_PUT(cd->intern, struct rt_db_internal);
	cd->intern = NULL;
	return -1;
    }

    if (cd->intern->idb_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	rt_db_free_internal(cd->intern);
	BU_PUT(cd->intern, struct rt_db_internal);
	cd->intern = NULL;
	return -1;
    }

    struct rt_bot_internal *bot = (struct rt_bot_internal *)cd->intern->idb_ptr;
    RT_BOT_CK_MAGIC(bot);

    lod->faces = bot->faces;
    lod->fcnt = (int)bot->num_faces;
    lod->points = (const point_t *)bot->vertices;
    lod->pcnt = (int)bot->num_vertices;
    lod->points_orig = (const point_t *)bot->vertices;
    lod->porig_cnt = (int)bot->num_vertices;
    lod->normals = NULL;

    return 0;
}

static int
rt_mesh_lod_bot_detail_clear(struct bsg_mesh_lod *lod, void *cb_data)
{
    struct rt_mesh_lod_bot_detail *cd = (struct rt_mesh_lod_bot_detail *)cb_data;

    if (cd && cd->intern) {
	rt_db_free_internal(cd->intern);
	BU_PUT(cd->intern, struct rt_db_internal);
	cd->intern = NULL;
    }

    if (lod) {
	lod->faces = NULL;
	lod->fcnt = 0;
	lod->points = NULL;
	lod->pcnt = 0;
	lod->points_orig = NULL;
	lod->porig_cnt = 0;
	lod->normals = NULL;
    }

    return 0;
}

static int
rt_mesh_lod_bot_detail_free(struct bsg_mesh_lod *lod, void *cb_data)
{
    rt_mesh_lod_bot_detail_clear(lod, cb_data);
    if (cb_data) {
	struct rt_mesh_lod_bot_detail *cd = (struct rt_mesh_lod_bot_detail *)cb_data;
	BU_PUT(cd, struct rt_mesh_lod_bot_detail);
    }

    return 0;
}

void
_rt_mesh_lod_context_destroy(struct rt_mesh_lod_context *c)
{
    if (!c)
	return;
    if (c->bsg) {
	bsg_mesh_lod_context_destroy(c->bsg);
	c->bsg = NULL;
    }
    BU_PUT(c, struct rt_mesh_lod_context);
}

void
db_mesh_lod_init(struct db_i *dbip, int verbose) {

    if (!dbip || !dbip->i)
	return;

    struct bsg_mesh_lod_context *mesh_c = rt_mesh_lod_context_ensure(dbip);
    if (!mesh_c)
	return;

    int64_t start, elapsed, overall_start;
    fastf_t seconds;
    dbip->i->mesh_c_completed = 0;
    dbip->i->mesh_c_target = 0;
    struct directory *dp;
    FOR_ALL_DIRECTORY_START(dp, dbip)
	if (dp->d_addr == RT_DIR_PHONY_ADDR)
	    continue;
	if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT)
	    dbip->i->mesh_c_target++;
    FOR_ALL_DIRECTORY_END;

    // Total target count is known, proceed
    start = bu_gettime();
    overall_start = start;
    FOR_ALL_DIRECTORY_START(dp, dbip)
	if (dp->d_addr == RT_DIR_PHONY_ADDR)
	    continue;
	if (dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
	    continue;

	// If we already have a match, assume it is valid.  Resetting
	// invalid data in the cache is outside the scope of cache init.
	unsigned long long key = bsg_mesh_lod_key_get(mesh_c, dp->d_namep);
	if (key)
	    continue;

	if (verbose > 1)
	    bu_log("Processing(%d):  %s\n", dbip->i->mesh_c_completed+1, dp->d_namep);

	// Process object
	db_mesh_lod_update(dbip, dp->d_namep);

	// Increment completed count
	dbip->i->mesh_c_completed++;

	elapsed = bu_gettime() - start;
	seconds = elapsed / 1000000.0;

	if (verbose > 1)
	    bu_log("Completed. (%g seconds)", seconds);

	if (seconds > 5.0) {
	    if (verbose) {
		elapsed = bu_gettime() - overall_start;
		seconds = elapsed / 1000000.0;
		bu_log("LoD cache processing (%g seconds): completed %d of %d BoTs\n", seconds, dbip->i->mesh_c_completed, dbip->i->mesh_c_target);
	    }
	    start = bu_gettime();
	}
    FOR_ALL_DIRECTORY_END;

    elapsed = bu_gettime() - overall_start;
    int rseconds = elapsed / 1000000;
    int rminutes = rseconds / 60;
    int rhours = rminutes / 60;
    rminutes = rminutes % 60;
    rseconds = rseconds % 60;
    bu_log("Mesh LoD caching complete (Elapsed time: %02d:%02d:%02d)\n", rhours, rminutes, rseconds);
}

void
db_mesh_lod_clear(struct db_i *dbip)
{
    if (!dbip || !dbip->i)
	return;

    struct bsg_mesh_lod_context *mesh_c = rt_mesh_lod_context_ensure(dbip);
    if (!mesh_c)
	return;

    bsg_mesh_lod_clear_cache(mesh_c, 0);
}

void
rt_mesh_lod_cache_clear_all(void)
{
    bsg_mesh_lod_clear_cache(NULL, 0);
}

void
rt_mesh_lod_cache_status_init(struct rt_mesh_lod_cache_status *status)
{
    struct rt_mesh_lod_cache_status defaults = RT_MESH_LOD_CACHE_STATUS_INIT;
    if (status)
	*status = defaults;
}

static int
rt_mesh_lod_payload_available(struct bsg_mesh_lod_context *mesh_c,
	unsigned long long key)
{
    struct bsg_mesh_lod *lod;

    if (!mesh_c || !key)
	return 0;

    lod = bsg_mesh_lod_create(mesh_c, key);
    if (!lod)
	return 0;

    bsg_mesh_lod_destroy(lod);
    return 1;
}

static void
rt_mesh_lod_status_current(struct db_i *dbip,
	struct bsg_mesh_lod_context *mesh_c,
	const char *name,
	struct rt_mesh_lod_cache_status *status)
{
    struct directory *dp;

    if (!status)
	return;

    rt_mesh_lod_cache_status_init(status);
    if (!dbip || !name)
	return;

    dp = db_lookup(dbip, name, LOOKUP_QUIET);
    status->directory_found = (dp != RT_DIR_NULL) ? 1 : 0;
    status->is_bot = (dp != RT_DIR_NULL &&
	    dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) ? 1 : 0;

    if (!mesh_c)
	return;

    status->cache_key = bsg_mesh_lod_key_get(mesh_c, name);
    status->has_cache_key = status->cache_key ? 1 : 0;
    status->has_cached_payload = rt_mesh_lod_payload_available(mesh_c,
	    status->cache_key);
    status->stale_cache_entry = (status->has_cache_key &&
	    !status->has_cached_payload) ? 1 : 0;
}

int
db_mesh_lod_status(struct db_i *dbip,
	const char *name,
	struct rt_mesh_lod_cache_status *status)
{
    if (!dbip || !dbip->i)
	return BRLCAD_ERROR;

    struct bsg_mesh_lod_context *mesh_c = rt_mesh_lod_context_ensure(dbip);
    if (!mesh_c) {
	bu_log("Error processing %s - unable to create LoD context\n",
		name ? name : "<all>");
	return BRLCAD_ERROR;
    }

    if (!name || !status)
	return BRLCAD_ERROR;

    rt_mesh_lod_status_current(dbip, mesh_c, name, status);
    return BRLCAD_OK;
}

int
db_mesh_lod_invalidate(struct db_i *dbip,
	const char *name,
	struct rt_mesh_lod_cache_status *status)
{
    struct rt_mesh_lod_cache_status current = RT_MESH_LOD_CACHE_STATUS_INIT;

    if (status)
	rt_mesh_lod_cache_status_init(status);

    if (!dbip || !dbip->i || !name)
	return BRLCAD_ERROR;

    struct bsg_mesh_lod_context *mesh_c = rt_mesh_lod_context_ensure(dbip);
    if (!mesh_c) {
	bu_log("Error invalidating %s - unable to create LoD context\n",
		name ? name : "<null>");
	return BRLCAD_ERROR;
    }

    rt_mesh_lod_status_current(dbip, mesh_c, name, &current);
    if (current.has_cache_key) {
	current.cleared_cache_entry = 1;
	current.cleared_cache_key = current.cache_key;
	bsg_mesh_lod_clear_cache(mesh_c, current.cache_key);
	current.cache_key = 0;
	current.has_cache_key = 0;
	current.has_cached_payload = 0;
	current.stale_cache_entry = 0;
    }

    if (status)
	*status = current;

    return BRLCAD_OK;
}

int
db_mesh_lod_refresh(struct db_i *dbip,
	const char *name,
	struct rt_mesh_lod_cache_status *status)
{
    struct rt_mesh_lod_cache_status current = RT_MESH_LOD_CACHE_STATUS_INIT;

    if (status)
	rt_mesh_lod_cache_status_init(status);

    if (!dbip || !dbip->i)
	return BRLCAD_ERROR;

    struct bsg_mesh_lod_context *mesh_c = rt_mesh_lod_context_ensure(dbip);
    if (!mesh_c) {
	bu_log("Error processing %s - unable to create LoD context\n",
		name ? name : "<all>");
	return BRLCAD_ERROR;
    }

    if (!name)
	return BRLCAD_OK;

    rt_mesh_lod_status_current(dbip, mesh_c, name, &current);

    if (current.has_cache_key) {
	current.cleared_cache_entry = 1;
	current.cleared_cache_key = current.cache_key;
	bsg_mesh_lod_clear_cache(mesh_c, current.cache_key);
	current.cache_key = 0;
	current.has_cache_key = 0;
	current.has_cached_payload = 0;
	current.stale_cache_entry = 0;
    }

    // If this isn't an active BoT, we're done.
    struct directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    current.directory_found = (dp != RT_DIR_NULL) ? 1 : 0;
    current.is_bot = (dp != RT_DIR_NULL &&
	    dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) ? 1 : 0;
    if (dp == RT_DIR_NULL || dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	if (status)
	    *status = current;
	return BRLCAD_OK;
    }

    // Unpack the data for LoD generation
    struct rt_db_internal dbintern;
    RT_DB_INTERNAL_INIT(&dbintern);
    struct rt_db_internal *ip = &dbintern;
    int ret = rt_db_get_internal(ip, dp, dbip, NULL);
    if (ret < 0) {
	bu_log("Error processing %s - internal get failed\n", dp->d_namep);
	return BRLCAD_ERROR;
    }
    if (ip->idb_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	bu_log("Error processing %s - mismatch between d_minor_type (%c) and idb_minor_type (%c)\n", dp->d_namep, dp->d_minor_type, ip->idb_minor_type);
	rt_db_free_internal(&dbintern);
	return BRLCAD_ERROR;
    }
    struct rt_bot_internal *bot = (struct rt_bot_internal *)ip->idb_ptr;
    RT_BOT_CK_MAGIC(bot);

    // Generate and write new data
    unsigned long long key = bsg_mesh_lod_cache(mesh_c, (const point_t *)bot->vertices, bot->num_vertices, NULL, bot->faces, bot->num_faces, 0, 0.66);
    if (!key) {
	bu_log("Error processing %s - unable to generate LoD data\n", dp->d_namep);
	rt_db_free_internal(&dbintern);
	return BRLCAD_ERROR;
    }
    if (bsg_mesh_lod_key_put(mesh_c, dp->d_namep, key) != 0) {
	bu_log("Error processing %s - unable to store LoD cache key\n", dp->d_namep);
	rt_db_free_internal(&dbintern);
	return BRLCAD_ERROR;
    }

    // Done with BoT
    rt_db_free_internal(&dbintern);

    // Make sure we can retrieve the cached data
    // TODO - may not really be necessary to verify this here once we're
    // working - including during early stages for testing.
    struct bsg_mesh_lod *lod = bsg_mesh_lod_create(mesh_c, key);
    if (!lod) {
	bu_log("Error processing %s - unable to retrieve LoD data\n", dp->d_namep);
	return BRLCAD_ERROR;
    }

    bsg_mesh_lod_destroy(lod);

    current.cache_key = key;
    current.has_cache_key = 1;
    current.has_cached_payload = 1;
    current.stale_cache_entry = 0;
    current.generated_cache_entry = 1;
    if (status)
	*status = current;

    return BRLCAD_OK;
}

int
db_mesh_lod_update(struct db_i *dbip, const char *name)
{
    return db_mesh_lod_refresh(dbip, name, NULL);
}

int
db_mesh_lod_store_mesh(struct db_i *dbip,
	const char *name,
	const point_t *vertices,
	size_t vertex_count,
	const vect_t *normals,
	const int *faces,
	size_t face_count,
	unsigned long long user_key,
	fastf_t fidelity_ratio,
	struct rt_mesh_lod_cache_status *status)
{
    struct rt_mesh_lod_cache_status current = RT_MESH_LOD_CACHE_STATUS_INIT;

    if (status)
	rt_mesh_lod_cache_status_init(status);

    if (!dbip || !dbip->i || !name || !vertices || !vertex_count ||
	    !faces || !face_count)
	return BRLCAD_ERROR;

    struct bsg_mesh_lod_context *mesh_c = rt_mesh_lod_context_ensure(dbip);
    if (!mesh_c) {
	bu_log("Error storing mesh LoD for %s - unable to create LoD context\n",
		name);
	return BRLCAD_ERROR;
    }

    rt_mesh_lod_status_current(dbip, mesh_c, name, &current);
    if (!current.directory_found) {
	if (status)
	    *status = current;
	return BRLCAD_ERROR;
    }

    if (current.has_cache_key &&
	    (current.cache_key != user_key || !current.has_cached_payload)) {
	current.cleared_cache_entry = 1;
	current.cleared_cache_key = current.cache_key;
	bsg_mesh_lod_clear_cache(mesh_c, current.cache_key);
	current.cache_key = 0;
	current.has_cache_key = 0;
	current.has_cached_payload = 0;
	current.stale_cache_entry = 0;
    }

    if (current.has_cache_key && current.has_cached_payload) {
	if (status)
	    *status = current;
	return BRLCAD_OK;
    }

    unsigned long long key = bsg_mesh_lod_cache(mesh_c, vertices, vertex_count,
	    normals, (int *)faces, face_count, user_key, fidelity_ratio);
    if (!key) {
	bu_log("Error storing mesh LoD for %s - unable to generate cache data\n",
		name);
	if (status)
	    *status = current;
	return BRLCAD_ERROR;
    }
    if (bsg_mesh_lod_key_put(mesh_c, name, key) != 0) {
	bu_log("Error storing mesh LoD for %s - unable to store cache key\n",
		name);
	if (status)
	    *status = current;
	return BRLCAD_ERROR;
    }

    current.cache_key = key;
    current.has_cache_key = 1;
    current.has_cached_payload = rt_mesh_lod_payload_available(mesh_c, key);
    current.stale_cache_entry = current.has_cached_payload ? 0 : 1;
    current.generated_cache_entry = current.has_cached_payload ? 1 : 0;
    if (status)
	*status = current;

    return current.has_cached_payload ? BRLCAD_OK : BRLCAD_ERROR;
}

struct rt_mesh_lod *
db_mesh_lod_get(struct db_i *dbip, const char *name)
{
    if (!dbip || !dbip->i || !name)
	return NULL;

    struct bsg_mesh_lod_context *mesh_c = rt_mesh_lod_context_ensure(dbip);
    if (!mesh_c)
	return NULL;

    struct directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
	return NULL;

    unsigned long long key = bsg_mesh_lod_key_get(mesh_c, name);
    if (!key)
	return NULL;

    struct bsg_mesh_lod *bsg_lod = bsg_mesh_lod_create(mesh_c, key);
    if (!bsg_lod)
	return NULL;

    struct rt_mesh_lod *lod;
    BU_GET(lod, struct rt_mesh_lod);
    lod->bsg = bsg_lod;

    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
	struct rt_mesh_lod_bot_detail *cbd;
	BU_GET(cbd, struct rt_mesh_lod_bot_detail);
	memset(cbd, 0, sizeof(*cbd));
	cbd->dbip = dbip;
	cbd->dp = dp;
	bsg_mesh_lod_detail_setup_clbk(bsg_lod, &rt_mesh_lod_bot_detail_setup, (void *)cbd);
	bsg_mesh_lod_detail_clear_clbk(bsg_lod, &rt_mesh_lod_bot_detail_clear);
	bsg_mesh_lod_detail_free_clbk(bsg_lod, &rt_mesh_lod_bot_detail_free);
    }

    return lod;
}

int
rt_mesh_lod_load_level(struct rt_mesh_lod *lod, int level, int reset)
{
    if (!lod || !lod->bsg)
	return -1;

    return bsg_mesh_lod_load_level(lod->bsg, level, reset);
}

int
rt_mesh_lod_load_view(struct rt_mesh_lod *lod, const struct rt_view_info *info, int reset)
{
    struct rt_view_info sanitized = RT_VIEW_INFO_INIT;

    if (!lod || !lod->bsg)
	return -1;

    if (info)
	sanitized = *info;
    rt_view_info_sanitize(&sanitized);

    return bsg_mesh_lod_load_view_size(lod->bsg, sanitized.size,
	    sanitized.lod.scale, reset);
}

int
rt_mesh_lod_current_level(const struct rt_mesh_lod *lod)
{
    if (!lod || !lod->bsg)
	return -1;

    return bsg_mesh_lod_load_level(lod->bsg, -1, 0);
}

int
rt_mesh_lod_has_active_data(const struct rt_mesh_lod *lod)
{
    if (!lod || !lod->bsg)
	return 0;

    const struct bsg_mesh_lod *bsg = lod->bsg;
    return (bsg->faces && bsg->fcnt > 0 && bsg->points && bsg->pcnt > 0) ?
	1 : 0;
}

int
rt_mesh_lod_data_get(const struct rt_mesh_lod *lod, struct rt_mesh_lod_data *data)
{
    if (!data)
	return 0;

    memset(data, 0, sizeof(*data));
    if (!lod || !lod->bsg)
	return 0;

    struct bsg_mesh_lod *bsg = lod->bsg;
    data->faces = bsg->faces;
    data->face_count = (bsg->fcnt > 0) ? (size_t)bsg->fcnt : 0;
    data->points = bsg->points;
    data->point_count = (bsg->pcnt > 0) ? (size_t)bsg->pcnt : 0;
    data->points_orig = bsg->points_orig;
    data->point_orig_count = (bsg->porig_cnt > 0) ? (size_t)bsg->porig_cnt : 0;
    data->normals = bsg->normals;
    data->normal_count = (bsg->normals && data->face_count) ?
	data->face_count * 3 : 0;
    VMOVE(data->bmin, bsg->bmin);
    VMOVE(data->bmax, bsg->bmax);

    return (data->faces && data->face_count &&
	    data->points && data->point_count);
}

void
rt_mesh_lod_info_init(struct rt_mesh_lod_info *info)
{
    struct rt_mesh_lod_info defaults = RT_MESH_LOD_INFO_INIT;
    if (info)
	*info = defaults;
}

int
rt_mesh_lod_info_get(const struct rt_mesh_lod *lod, struct rt_mesh_lod_info *info)
{
    if (!info)
	return 0;

    rt_mesh_lod_info_init(info);
    if (!lod || !lod->bsg)
	return 0;

    struct bsg_mesh_lod *bsg = lod->bsg;
    info->active_level = bsg_mesh_lod_load_level(bsg, -1, 0);
    info->face_count = (bsg->fcnt > 0) ? (size_t)bsg->fcnt : 0;
    info->point_count = (bsg->pcnt > 0) ? (size_t)bsg->pcnt : 0;
    info->point_orig_count = (bsg->porig_cnt > 0) ? (size_t)bsg->porig_cnt : 0;
    info->normal_count = (bsg->normals && info->face_count) ?
	info->face_count * 3 : 0;
    info->has_faces = (bsg->faces && info->face_count) ? 1 : 0;
    info->has_points = (bsg->points && info->point_count) ? 1 : 0;
    info->has_original_points = (bsg->points_orig && info->point_orig_count) ? 1 : 0;
    info->has_snapped_points = (info->has_points && info->has_original_points &&
	    bsg->points != bsg->points_orig) ? 1 : 0;
    info->has_normals = (bsg->normals && info->normal_count) ? 1 : 0;
    VMOVE(info->bmin, bsg->bmin);
    VMOVE(info->bmax, bsg->bmax);

    return info->has_faces && info->has_points;
}

void
rt_mesh_lod_detail_init(struct rt_mesh_lod_detail *detail)
{
    if (detail)
	memset(detail, 0, sizeof(*detail));
}

int
rt_mesh_lod_detail_callbacks_set(struct rt_mesh_lod *lod,
	rt_mesh_lod_detail_setup_callback setup_clbk,
	rt_mesh_lod_detail_clear_callback clear_clbk,
	rt_mesh_lod_detail_free_callback free_clbk,
	void *cb_data)
{
    if (!lod || !lod->bsg || !setup_clbk)
	return 0;

    struct rt_mesh_lod_detail_callbacks *callbacks;
    BU_GET(callbacks, struct rt_mesh_lod_detail_callbacks);
    callbacks->setup_clbk = setup_clbk;
    callbacks->clear_clbk = clear_clbk;
    callbacks->free_clbk = free_clbk;
    callbacks->cb_data = cb_data;

    bsg_mesh_lod_detail_setup_clbk(lod->bsg,
	    &rt_mesh_lod_detail_setup_bsg, (void *)callbacks);
    bsg_mesh_lod_detail_clear_clbk(lod->bsg,
	    &rt_mesh_lod_detail_clear_bsg);
    bsg_mesh_lod_detail_free_clbk(lod->bsg,
	    &rt_mesh_lod_detail_free_bsg);

    return 1;
}

void
rt_mesh_lod_memshrink(struct rt_mesh_lod *lod)
{
    if (!lod || !lod->bsg)
	return;

    bsg_mesh_lod_shrink_memory(lod->bsg);
}

void
rt_mesh_lod_destroy(struct rt_mesh_lod *lod)
{
    if (!lod)
	return;
    if (lod->bsg) {
	bsg_mesh_lod_destroy(lod->bsg);
	lod->bsg = NULL;
    }
    BU_PUT(lod, struct rt_mesh_lod);
}

struct bsg_mesh_lod *
_rt_mesh_lod_bsg(struct rt_mesh_lod *lod)
{
    return lod ? lod->bsg : NULL;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
