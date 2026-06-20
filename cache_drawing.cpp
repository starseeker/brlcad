/*               C A C H E _ D R A W I N G . C P P
 * BRL-CAD
 *
 * Copyright (c) 2016-2025 United States Government as represented by
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
/** @file cache_drawing.cpp
 *
 * 5-stage concurrent pipeline for pre-computing drawing data:
 *
 *   [main thread] → q_init (object names)
 *       ↓ attr_worker  — cracks rt_db_internal, writes region attrs to dcache
 *       → q_aabb
 *       ↓ aabb_worker  — ft_bbox, writes AABB to dcache, posts DRAWRESULT_AABB
 *       → q_obb
 *       ↓ obb_worker   — ft_oriented_bbox, writes OBB to dcache,
 *                        posts DRAWRESULT_OBB
 *       → q_lod
 *       ↓ lod_worker   — bsg_mesh_lod_cache for BoTs, posts DRAWRESULT_LOD
 *
 *   [write_worker]  — drains q_write and serialises all LMDB writes
 *
 * All DrawResult notifications are posted to results_q, which the main
 * thread drains via DrawPipeline::drain() (called from
 * DbiState::drain_geom_results).
 *
 * Cache key format: "<hash>:<component>", e.g. "12345678:aabb"
 * where <hash> = bu_data_hash(dp->d_namep, ...).
 */

#include "common.h"

#include <atomic>
#include <chrono>
#include <cstring>
#include <string>
#include <thread>

#include "bu/app.h"
#include "bu/cache.h"
#include "bu/color.h"
#include "bu/file.h"
#include "bu/opt.h"
#include "bu/time.h"
#include "bsg/lod.h"
#include "raytrace.h"

#include "./librt_private.h"

/* Cache key component identifiers.  Must match those used in dbi_state.cpp. */
#define CACHE_OBJ_BOUNDS   "bb"
#define CACHE_OBB          "obb"
#define CACHE_REGION_ID    "rid"
#define CACHE_REGION_FLAG  "rf"
#define CACHE_INHERIT_FLAG "if"
#define CACHE_COLOR        "c"
#define CACHE_TIMESTAMP    "timestmp"

/* ------------------------------------------------------------------
 * Debug delay helpers — read once at pipeline startup.
 *
 * Set environment variables to slow down individual pipeline stages so
 * intermediate draw states (AABB placeholder, OBB placeholder, …) are
 * visible long enough to screenshot for debugging:
 *
 *   BRLCAD_CACHE_ATTR_DELAY_MS  — sleep after each attr_worker object
 *   BRLCAD_CACHE_AABB_DELAY_MS  — sleep after each aabb_worker object
 *   BRLCAD_CACHE_OBB_DELAY_MS   — sleep after each obb_worker  object
 *   BRLCAD_CACHE_LOD_DELAY_MS   — sleep after each lod_worker  object
 *
 * Example: BRLCAD_CACHE_AABB_DELAY_MS=500 adds a 500 ms pause after
 * each AABB computation, making the AABB wireframe placeholder state
 * linger for many seconds on a model with dozens of BoTs.
 * ------------------------------------------------------------------ */
static int g_delay_attr_ms = 0;
static int g_delay_aabb_ms = 0;
static int g_delay_obb_ms  = 0;
static int g_delay_lod_ms  = 0;

static void
init_debug_delays(void)
{
    static bool done = false;
    if (done) return;
    done = true;
    auto getenv_int = [](const char *name) -> int {
	const char *v = getenv(name);
	if (!v) return 0;
	int n = atoi(v);
	return (n > 0) ? n : 0;
    };
    g_delay_attr_ms = getenv_int("BRLCAD_CACHE_ATTR_DELAY_MS");
    g_delay_aabb_ms = getenv_int("BRLCAD_CACHE_AABB_DELAY_MS");
    g_delay_obb_ms  = getenv_int("BRLCAD_CACHE_OBB_DELAY_MS");
    g_delay_lod_ms  = getenv_int("BRLCAD_CACHE_LOD_DELAY_MS");
    if (g_delay_attr_ms || g_delay_aabb_ms ||
	g_delay_obb_ms  || g_delay_lod_ms)
	bu_log("DrawPipeline: debug delays enabled "
	       "(attr=%d aabb=%d obb=%d lod=%d ms)\n",
	       g_delay_attr_ms, g_delay_aabb_ms,
	       g_delay_obb_ms,  g_delay_lod_ms);
}

/* ------------------------------------------------------------------
 * CacheWriteItem implementation
 * ------------------------------------------------------------------ */

CacheWriteItem::CacheWriteItem()
    : erase_op(false), data_len(0), data(nullptr)
{
    key[0] = '\0';
}

CacheWriteItem::CacheWriteItem(const char *k, const void *d, size_t len)
{
    snprintf(key, BU_CACHE_KEY_MAXLEN, "%s", k);
    erase_op = (!d || !len);
    data_len  = len;
    if (d && len) {
	data = bu_malloc(len, "CacheWriteItem data");
	memcpy(data, d, len);
    } else {
	data = nullptr;
    }
}

CacheWriteItem::CacheWriteItem(const CacheWriteItem &o)
    : erase_op(o.erase_op), data_len(o.data_len), data(nullptr)
{
    snprintf(key, BU_CACHE_KEY_MAXLEN, "%s", o.key);
    if (o.data && o.data_len) {
	data = bu_malloc(o.data_len, "CacheWriteItem copy");
	memcpy(data, o.data, o.data_len);
    }
}

CacheWriteItem &
CacheWriteItem::operator=(const CacheWriteItem &o)
{
    if (this == &o)
	return *this;
    if (data)
	bu_free(data, "CacheWriteItem free");
    data = nullptr;
    erase_op = o.erase_op;
    data_len  = o.data_len;
    snprintf(key, BU_CACHE_KEY_MAXLEN, "%s", o.key);
    if (o.data && o.data_len) {
	data = bu_malloc(o.data_len, "CacheWriteItem assign");
	memcpy(data, o.data, o.data_len);
    }
    return *this;
}

CacheWriteItem::~CacheWriteItem()
{
    if (data)
	bu_free(data, "CacheWriteItem dtor");
    data = nullptr;
}

/* ------------------------------------------------------------------
 * Helper: build a cache key string
 * ------------------------------------------------------------------ */
static inline void
make_key(char *buf, unsigned long long hash, const char *component)
{
    snprintf(buf, BU_CACHE_KEY_MAXLEN, "%llu:%s", hash, component);
}

/* ------------------------------------------------------------------
 * Pipeline stage 1: attr_worker
 *
 * Consumes object names from q_init.
 * - Cracks rt_db_internal using a per-thread resource.
 * - Writes region_flag, region_id, inherit_flag, color to q_write.
 * - Passes rt_db_internal to q_aabb.
 * ------------------------------------------------------------------ */
static void
attr_worker(std::shared_ptr<DrawPipelineState> p)
{
    char ckey[BU_CACHE_KEY_MAXLEN];
    struct resource bres;
    memset(&bres, 0, sizeof(bres));
    rt_init_resource(&bres, 1, NULL);
    init_debug_delays();

    while (!p->shutdown) {
	if (p->q_init.size_approx() == 0) {
	    std::this_thread::sleep_for(std::chrono::milliseconds(10));
	    continue;
	}

	std::string name;
	if (!p->q_init.try_dequeue(name))
	    continue;

	struct directory *dp = db_lookup(p->dbip, name.c_str(), LOOKUP_QUIET);
	if (dp == RT_DIR_NULL)
	    continue;

	unsigned long long hash =
	    bu_data_hash(dp->d_namep, strlen(dp->d_namep) * sizeof(char));

	/* Read attributes (non-blocking, doesn't need the rt_db_internal) */
	struct bu_attribute_value_set avs = BU_AVS_INIT_ZERO;
	if (db5_get_attributes(p->dbip, &avs, dp) < 0)
	    bu_avs_free(&avs);

	/* --- Region flag --- */
	int rflag = 0;
	{
	    const char *s = bu_avs_get(&avs, "region");
	    if (s && (BU_STR_EQUAL(s, "R") || BU_STR_EQUAL(s, "1")))
		rflag = 1;
	}
	make_key(ckey, hash, CACHE_REGION_FLAG);
	p->q_write.enqueue(CacheWriteItem(ckey, &rflag, sizeof(int)));

	/* --- Region id --- */
	int region_id = -1;
	{
	    const char *s = bu_avs_get(&avs, "region_id");
	    if (s)
		bu_opt_int(NULL, 1, &s, (void *)&region_id);
	}
	make_key(ckey, hash, CACHE_REGION_ID);
	p->q_write.enqueue(CacheWriteItem(ckey, &region_id, sizeof(int)));

	/* --- Inherit flag --- */
	int inherit = 0;
	{
	    const char *s = bu_avs_get(&avs, "inherit");
	    if (s && BU_STR_EQUAL(s, "1"))
		inherit = 1;
	}
	make_key(ckey, hash, CACHE_INHERIT_FLAG);
	p->q_write.enqueue(CacheWriteItem(ckey, &inherit, sizeof(int)));

	/* --- Color --- */
	unsigned int colors = UINT_MAX;
	{
	    const char *s = bu_avs_get(&avs, "color");
	    if (!s) s = bu_avs_get(&avs, "rgb");
	    if (s) {
		struct bu_color col;
		bu_opt_color(NULL, 1, &s, (void *)&col);
		int r, g, b;
		bu_color_to_rgb_ints(&col, &r, &g, &b);
		colors = (unsigned int)(r + (g << 8) + (b << 16));
	    }
	}
	make_key(ckey, hash, CACHE_COLOR);
	p->q_write.enqueue(CacheWriteItem(ckey, &colors, sizeof(unsigned int)));

	bu_avs_free(&avs);

	/* Crack the geometry and pass to q_aabb */
	struct rt_db_internal *ip;
	BU_GET(ip, struct rt_db_internal);
	RT_DB_INTERNAL_INIT(ip);
	if (rt_db_get_internal(ip, dp, p->dbip, NULL, &bres) < 0) {
	    BU_PUT(ip, struct rt_db_internal);
	    continue;
	}
	/* Register the name before enqueuing so downstream stages can look it up. */
	{
	    std::lock_guard<std::mutex> lk(p->name_mu);
	    p->ip_names[ip] = std::string(dp->d_namep);
	}
	if (g_delay_attr_ms > 0)
	    std::this_thread::sleep_for(
		std::chrono::milliseconds(g_delay_attr_ms));
	p->q_aabb.enqueue(ip);
    }

    rt_clean_resource_basic(NULL, &bres);
    p->thread_cnt--;
}

/* ------------------------------------------------------------------
 * Pipeline stage 2: aabb_worker
 *
 * Consumes rt_db_internal from q_aabb.
 * - Calls ft_bbox to get AABB.
 * - Writes AABB data to q_write.
 * - Posts DRAWRESULT_AABB to results_q.
 * - Passes rt_db_internal to q_obb.
 * ------------------------------------------------------------------ */
static void
aabb_worker(std::shared_ptr<DrawPipelineState> p)
{
    char ckey[BU_CACHE_KEY_MAXLEN];
    const struct bn_tol btol = BN_TOL_INIT_TOL;

    while (!p->shutdown) {
	if (p->q_aabb.size_approx() == 0) {
	    std::this_thread::sleep_for(std::chrono::milliseconds(10));
	    continue;
	}

	struct rt_db_internal *ip = nullptr;
	if (!p->q_aabb.try_dequeue(ip))
	    continue;

	std::string ip_name;
	{
	    std::lock_guard<std::mutex> lk(p->name_mu);
	    auto it = p->ip_names.find(ip);
	    if (it == p->ip_names.end()) { p->q_obb.enqueue(ip); continue; }
	    ip_name = it->second;
	}
	const char *name = ip_name.c_str();
	unsigned long long hash =
	    bu_data_hash(name, strlen(name) * sizeof(char));

	make_key(ckey, hash, CACHE_OBJ_BOUNDS);

	DrawResult dr; memset(&dr, 0, sizeof(dr));
	dr.type = DRAWRESULT_AABB;
	dr.hash = hash;
	snprintf(dr.dp_name, sizeof(dr.dp_name), "%s", name);

	if (ip->idb_meth && ip->idb_meth->ft_bbox) {
	    point_t bmin, bmax;
	    VSETALL(bmin,  INFINITY);
	    VSETALL(bmax, -INFINITY);
	    if (ip->idb_meth->ft_bbox(ip, &bmin, &bmax, &btol) == 0) {
		/* Store as two consecutive point_t values */
		point_t bb[2];
		VMOVE(bb[0], bmin);
		VMOVE(bb[1], bmax);
		p->q_write.enqueue(
		    CacheWriteItem(ckey, &bb, 2 * sizeof(point_t)));

		VMOVE(dr.bmin, bmin);
		VMOVE(dr.bmax, bmax);
		p->results_q.enqueue(dr);
	    } else {
		/* Clear any stale entry */
		p->q_write.enqueue(CacheWriteItem(ckey, nullptr, 0));
	    }
	} else {
	    p->q_write.enqueue(CacheWriteItem(ckey, nullptr, 0));
	}

	if (g_delay_aabb_ms > 0)
	    std::this_thread::sleep_for(
		std::chrono::milliseconds(g_delay_aabb_ms));
	p->q_obb.enqueue(ip);
    }

    p->thread_cnt--;
}

/* ------------------------------------------------------------------
 * Pipeline stage 3: obb_worker
 *
 * Consumes rt_db_internal from q_obb.
 * - Calls ft_oriented_bbox to get an arb8 OBB.
 * - Writes OBB (8 corners) to q_write.
 * - Posts DRAWRESULT_OBB to results_q.
 * - Passes rt_db_internal to q_lod.
 * ------------------------------------------------------------------ */
static void
obb_worker(std::shared_ptr<DrawPipelineState> p)
{
    char ckey[BU_CACHE_KEY_MAXLEN];

    while (!p->shutdown) {
	if (p->q_obb.size_approx() == 0) {
	    std::this_thread::sleep_for(std::chrono::milliseconds(10));
	    continue;
	}

	struct rt_db_internal *ip = nullptr;
	if (!p->q_obb.try_dequeue(ip))
	    continue;

	std::string ip_name;
	{
	    std::lock_guard<std::mutex> lk(p->name_mu);
	    auto it = p->ip_names.find(ip);
	    if (it == p->ip_names.end()) { p->q_lod.enqueue(ip); continue; }
	    ip_name = it->second;
	}
	const char *name = ip_name.c_str();
	unsigned long long hash =
	    bu_data_hash(name, strlen(name) * sizeof(char));

	make_key(ckey, hash, CACHE_OBB);

	DrawResult dr; memset(&dr, 0, sizeof(dr));
	dr.type      = DRAWRESULT_OBB;
	dr.hash      = hash;
	dr.obb_valid = 0;
	snprintf(dr.dp_name, sizeof(dr.dp_name), "%s", name);

	if (ip->idb_meth && ip->idb_meth->ft_oriented_bbox) {
	    struct rt_arb_internal arb;
	    arb.magic = RT_ARB_INTERNAL_MAGIC;
	    double tol_dist = BN_TOL_DIST;
	    if (ip->idb_meth->ft_oriented_bbox(&arb, ip, tol_dist) == 0) {
		/* The arb8 directly provides the 8 OBB corner points */
		for (int k = 0; k < 8; k++)
		    VMOVE(dr.obb_pts[k], arb.pt[k]);
		dr.obb_valid = 1;
		p->q_write.enqueue(
		    CacheWriteItem(ckey, arb.pt, 8 * sizeof(point_t)));
		p->results_q.enqueue(dr);
	    } else {
		p->q_write.enqueue(CacheWriteItem(ckey, nullptr, 0));
	    }
	}
	/* No ft_oriented_bbox → skip OBB silently, still forward to lod */

	if (g_delay_obb_ms > 0)
	    std::this_thread::sleep_for(
		std::chrono::milliseconds(g_delay_obb_ms));
	p->q_lod.enqueue(ip);
    }

    p->thread_cnt--;
}

/* ------------------------------------------------------------------
 * Pipeline stage 4: lod_worker
 *
 * Consumes rt_db_internal from q_lod.
 * - For BoTs: calls bsg_mesh_lod_cache + bsg_mesh_lod_key_put.
 * - Posts DRAWRESULT_LOD to results_q.
 * - Frees the rt_db_internal (last consumer).
 * ------------------------------------------------------------------ */
static void
lod_worker(std::shared_ptr<DrawPipelineState> p)
{
    while (!p->shutdown) {
	if (p->q_lod.size_approx() == 0) {
	    std::this_thread::sleep_for(std::chrono::milliseconds(10));
	    continue;
	}

	struct rt_db_internal *ip = nullptr;
	if (!p->q_lod.try_dequeue(ip))
	    continue;

	std::string ip_name;
	{
	    std::lock_guard<std::mutex> lk(p->name_mu);
	    auto it = p->ip_names.find(ip);
	    if (it != p->ip_names.end())
		ip_name = it->second;
	    p->ip_names.erase(ip); /* last consumer — remove from map */
	}
	const char *name = ip_name.c_str();

	if (p->lod_ctx && !ip_name.empty() &&
	    ip->idb_minor_type == DB5_MINORTYPE_BRLCAD_BOT)
	{
	    struct rt_bot_internal *bot =
		(struct rt_bot_internal *)ip->idb_ptr;
	    if (bot && bot->magic == RT_BOT_INTERNAL_MAGIC
		&& bot->num_faces > 0 && bot->num_vertices > 0)
	    {
		unsigned long long hash =
		    bu_data_hash(name, strlen(name) * sizeof(char));

		/* Fast path: if the namekey db already has a key for this
		 * object the LoD data is already cached.  bsg_mesh_lod_key_get
		 * now uses a per-call MDB_RDONLY transaction (no mutex) so this
		 * check is cheap even under concurrent read traffic. */
		unsigned long long key =
		    bsg_mesh_lod_key_get(p->lod_ctx, name);
		if (!key) {
		    /* Cache miss — compute and store POP data.  This is the
		     * expensive path: it hashes all vertex/face data, runs
		     * tri_process(), and writes the results to LMDB. */
		    key = bsg_mesh_lod_cache(
			p->lod_ctx,
			(const point_t *)bot->vertices, bot->num_vertices,
			NULL, bot->faces, bot->num_faces, 0, 0.66);
		    if (key)
			bsg_mesh_lod_key_put(p->lod_ctx, name, key);
		}
		if (key) {
		    DrawResult dr; memset(&dr, 0, sizeof(dr));
		    dr.type    = DRAWRESULT_LOD;
		    dr.hash    = hash;
		    dr.lod_key = key;
		    snprintf(dr.dp_name, sizeof(dr.dp_name), "%s", name);
		    p->results_q.enqueue(dr);
		}
	    }
	}

	/* LoD worker is the last stage — free the internal. */
	if (g_delay_lod_ms > 0)
	    std::this_thread::sleep_for(
		std::chrono::milliseconds(g_delay_lod_ms));
	rt_db_free_internal(ip);
	BU_PUT(ip, struct rt_db_internal);
    }

    p->thread_cnt--;
}

/* ------------------------------------------------------------------
 * Pipeline stage 5: write_worker
 *
 * Serialises all LMDB writes from q_write.  A single writer is required
 * because lmdb opens in MDB_NOLOCK mode by default in bu_cache.
 * ------------------------------------------------------------------ */
static void
write_worker(std::shared_ptr<DrawPipelineState> p)
{
    while (!p->shutdown) {
	if (p->q_write.size_approx() == 0) {
	    std::this_thread::sleep_for(std::chrono::milliseconds(10));
	    continue;
	}

	CacheWriteItem item;
	if (!p->q_write.try_dequeue(item))
	    continue;

	if (!p->dcache)
	    continue;

	if (item.erase_op || !item.data || !item.data_len) {
	    bu_cache_clear(item.key, p->dcache, NULL);
	} else {
	    int tries = 0;
	    while (tries < 5 &&
		   !bu_cache_write(item.data, item.data_len,
				   item.key, p->dcache, NULL))
	    {
		tries++;
		std::this_thread::sleep_for(std::chrono::milliseconds(50));
	    }
	}
    }

    p->thread_cnt--;
}

/* ------------------------------------------------------------------
 * Public extern "C" API
 * ------------------------------------------------------------------ */

extern "C" int
db_cache_start(struct db_i *dbip)
{
    if (!dbip || !dbip->i)
	return BRLCAD_ERROR;

    /* Already running */
    if (dbip->i->draw_pipeline)
	return BRLCAD_OK;

    /* In-memory databases have no filename → skip LMDB cache */
    if (dbip->dbi_filename && strlen(dbip->dbi_filename)) {
	/* Size the LMDB map to 2× the .g file size.  The drawing cache stores
	 * preprocessed data (bboxes, OBBs, attr copies) that is typically
	 * smaller than the raw geometry, but a 2× multiplier gives headroom
	 * for growth without requiring frequent remap operations. */
	long long fsize = 2 * bu_file_size(dbip->dbi_filename);
	if (fsize <= 0)
	    fsize = 0;
	dbip->i->dcache = bu_cache_open(dbip->dbi_filename, 1,
					(size_t)fsize);
	/* dcache being NULL is non-fatal — we just skip persistence */
    }

    auto state = std::make_shared<DrawPipelineState>();
    state->dbip   = dbip;
    state->dcache = dbip->i->dcache;
    /* lod_ctx is NULL at this point; set later via db_cache_set_lod_ctx */
    state->lod_ctx = nullptr;

    /* Launch the 5 pipeline threads */
    state->thread_cnt = 5;
    state->threads.emplace_back(attr_worker,  state);
    state->threads.emplace_back(aabb_worker,  state);
    state->threads.emplace_back(obb_worker,   state);
    state->threads.emplace_back(lod_worker,   state);
    state->threads.emplace_back(write_worker, state);

    /* Detach all threads — they run until shutdown is signalled */
    for (auto &t : state->threads)
	t.detach();
    state->threads.clear();

    /* Store the shared_ptr as a raw pointer in the opaque void* slot.
     * db_cache_stop() will recover and reset it. */
    dbip->i->draw_pipeline =
	new std::shared_ptr<DrawPipelineState>(state);

    /* Enqueue all existing directory entries for initial processing */
    struct directory *dp;
    for (int i = 0; i < RT_DBNHASH; i++) {
	for (dp = dbip->dbi_Head[i]; dp != RT_DIR_NULL; dp = dp->d_forw) {
	    if (dp->d_major_type != DB5_MAJORTYPE_BRLCAD)
		continue;
	    if (dp->d_addr == RT_DIR_PHONY_ADDR)
		continue;
	    if (!dp->d_namep || !strlen(dp->d_namep))
		continue;
	    /* Skip combs — they don't have geometry to bbox/lod */
	    if (dp->d_flags & RT_DIR_COMB)
		continue;
	    state->q_init.enqueue(std::string(dp->d_namep));
	}
    }

    return BRLCAD_OK;
}

extern "C" void
db_cache_stop(struct db_i *dbip)
{
    if (!dbip || !dbip->i || !dbip->i->draw_pipeline)
	return;

    auto *pp = static_cast<std::shared_ptr<DrawPipelineState> *>(
	dbip->i->draw_pipeline);
    if (pp && *pp) {
	(*pp)->shutdown = true;
	/* Give threads a moment to notice the shutdown flag */
	int waited = 0;
	while ((*pp)->thread_cnt > 0 && waited < 200) {
	    std::this_thread::sleep_for(std::chrono::milliseconds(10));
	    waited++;
	}
    }

    delete pp;
    dbip->i->draw_pipeline = nullptr;

    if (dbip->i->dcache) {
	bu_cache_close(dbip->i->dcache);
	dbip->i->dcache = nullptr;
    }
}

extern "C" void
db_cache_queue_obj(struct db_i *dbip, const char *name)
{
    if (!dbip || !dbip->i || !dbip->i->draw_pipeline || !name)
	return;
    auto *pp = static_cast<std::shared_ptr<DrawPipelineState> *>(
	dbip->i->draw_pipeline);
    if (pp && *pp)
	(*pp)->q_init.enqueue(std::string(name));
}

extern "C" int
db_cache_settled(struct db_i *dbip)
{
    if (!dbip || !dbip->i || !dbip->i->draw_pipeline)
	return 1;
    auto *pp = static_cast<std::shared_ptr<DrawPipelineState> *>(
	dbip->i->draw_pipeline);
    if (!pp || !*pp)
	return 1;
    auto &s = **pp;
    return (s.q_init.size_approx()  == 0 &&
	    s.q_aabb.size_approx()  == 0 &&
	    s.q_obb.size_approx()   == 0 &&
	    s.q_lod.size_approx()   == 0 &&
	    s.q_write.size_approx() == 0) ? 1 : 0;
}

extern "C" void
db_cache_set_lod_ctx(struct db_i *dbip, bsg_mesh_lod_context *lod_ctx)
{
    if (!dbip || !dbip->i || !dbip->i->draw_pipeline)
	return;
    auto *pp = static_cast<std::shared_ptr<DrawPipelineState> *>(
	dbip->i->draw_pipeline);
    if (pp && *pp)
	(*pp)->lod_ctx = lod_ctx;
}

/* Drain all pending DrawResult notifications into out[].
 * Returns the number of results drained.
 * MAIN THREAD ONLY. */
extern "C" size_t
db_cache_drain_results(struct db_i *dbip, void *out_vec)
{
    if (!dbip || !dbip->i || !dbip->i->draw_pipeline || !out_vec)
	return 0;
    auto *pp = static_cast<std::shared_ptr<DrawPipelineState> *>(
	dbip->i->draw_pipeline);
    if (!pp || !*pp)
	return 0;

    auto *results =
	static_cast<std::vector<DrawResult> *>(out_vec);
    size_t n = 0;
    DrawResult dr; memset(&dr, 0, sizeof(dr));
    while ((*pp)->results_q.try_dequeue(dr)) {
	results->push_back(dr);
	n++;
    }
    return n;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
