/*                         L O D . C P P
 * BRL-CAD
 *
 * Copyright (c) 2008-2026 United States Government as represented by
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
/** @file libged/view/lod.cpp
 *
 * The view lod subcommand.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "BObol/BDrawCache.h"
#include "BObol/BMeshLodCache.h"
#include "bu/cmd.h"
#include "bu/hash.h"
#include "bu/snooze.h"
#include "bu/str.h"
#include "bu/time.h"
#include "bu/vls.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "rt/view.h"

#include "../ged_private.h"
#include "./ged_view.h"

static int
lod_service_status_has_work(
    const ged_draw_obol_lod_service_status_t *status)
{
    if (!status)
	return 0;

    return status->in_flight > 0 ||
	   status->pending_tasks > 0 ||
	   status->queued_results > 0 ||
	   status->queued_cache_writes > 0 ||
	   status->delayed_tasks > 0;
}

int
_view_cmd_lod(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] lod [subcommand] [vals]";
    const char *purpose_string = "manage Level of Detail drawing settings";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct ged *gedp = gd->gedp;
    void *view_ctx;
    int print_help = 0;
    static const char *usage = "view lod [csg|mesh] [0|1]\n"
			       "view lod cache [clear [all_files] | exists] \n"
			       "view lod service [status|start [workers]|stop|poll [max_results]|wait [timeout_ms] [max_results]|prewarm [all|bot ...]]\n"
			       "view lod scale [factor]\n"
			       "view lod point_scale [factor]\n"
			       "view lod curve_scale [factor]\n"
			       "view lod bot_threshold [face_cnt]\n";

    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bu_opt_desc d[2];
    BU_OPT(d[0], "h", "help",        "",  NULL, &print_help,      "Print help");
    BU_OPT_NULL(d[1]);

    // We know we're the lod command - start processing args
    argc--;
    argv++;

    int ac = bu_opt_parse(NULL, argc, argv, d);
    argc = ac;

    if (print_help) {
	bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
	return GED_HELP;
    }

    if (argc > 4 ||
	(argc > 3 && !(argc >= 2 && BU_STR_EQUAL(argv[0], "service") &&
		       (BU_STR_EQUAL(argv[1], "prewarm") ||
			BU_STR_EQUAL(argv[1], "wait"))))) {
	bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
	return BRLCAD_ERROR;
    }

    view_ctx = gd->cv;
    if (view_ctx == NULL) {
	bu_vls_printf(gedp->ged_result_str, "no current view defined\n");
	return BRLCAD_ERROR;
    }

    ged_draw_view_lod_policy lod_policy;
    if (!ged_draw_view_context_lod_policy_get(&lod_policy, view_ctx)) {
	bu_vls_printf(gedp->ged_result_str, "unable to read LoD policy\n");
	return BRLCAD_ERROR;
    }
    auto commit_lod_policy = [&]() {
	ged_draw_view_context_lod_policy_apply(view_ctx, &lod_policy);
	if (!ged_draw_obol_view_lod_policy_changed(gedp, view_ctx)) {
	    int rac = 1;
	    const char *rav[1] = {"redraw"};
	    ged_exec_redraw(gedp, rac, (const char **)rav);
	}
    };
    auto redraw_view = [&]() {
	int rac = 1;
	const char *rav[1] = {"redraw"};
	ged_exec_redraw(gedp, rac, (const char **)rav);
    };

    /* Print current state if no args are supplied */
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "enabled(mesh/csg): %d/%d\n", lod_policy.mesh_enabled, lod_policy.csg_enabled);
	bu_vls_printf(gedp->ged_result_str, "scale: %g\n", lod_policy.scale);
	bu_vls_printf(gedp->ged_result_str, "point_scale: %g\n", lod_policy.point_scale);
	bu_vls_printf(gedp->ged_result_str, "curve_scale: %g\n", lod_policy.curve_scale);
	bu_vls_printf(gedp->ged_result_str, "bot_threshold: %zu\n", lod_policy.bot_threshold);
	return BRLCAD_OK;
    }

    if (BU_STR_EQUIV(argv[0], "1")) {
	if (!lod_policy.mesh_enabled || !lod_policy.csg_enabled) {
	    lod_policy.mesh_enabled = 1;
	    lod_policy.csg_enabled = 1;
	    lod_policy.zoom_refresh = 1;
	    commit_lod_policy();
	}
	return BRLCAD_OK;
    }
    if (BU_STR_EQUIV(argv[0], "0")) {
	if (lod_policy.mesh_enabled || lod_policy.csg_enabled) {
	    lod_policy.mesh_enabled = 0;
	    lod_policy.csg_enabled = 0;
	    lod_policy.zoom_refresh = 0;
	    commit_lod_policy();
	}
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "csg")) {
	if (argc == 1) {
	    bu_vls_printf(gedp->ged_result_str, "%d\n", lod_policy.csg_enabled);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL(argv[1], "1")) {
	    if (!lod_policy.csg_enabled) {
		lod_policy.csg_enabled = 1;
		lod_policy.zoom_refresh = 1;
		commit_lod_policy();
	    }
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL(argv[1], "0")) {
	    if (lod_policy.csg_enabled) {
		lod_policy.csg_enabled = 0;
		if (!lod_policy.mesh_enabled)
		    lod_policy.zoom_refresh = 0;
		commit_lod_policy();
	    }
	    return BRLCAD_OK;
	}
	bu_vls_printf(gedp->ged_result_str, "Error - invalid arg: \"%s\".  Valid args are 0 or 1", argv[1]);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(argv[0], "mesh")) {
	if (argc == 1) {
	    bu_vls_printf(gedp->ged_result_str, "%d\n", lod_policy.mesh_enabled);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL(argv[1], "1")) {
	    if (!lod_policy.mesh_enabled) {
		lod_policy.mesh_enabled = 1;
		lod_policy.zoom_refresh = 1;
		commit_lod_policy();
	    }
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL(argv[1], "0")) {
	    if (lod_policy.mesh_enabled) {
		lod_policy.mesh_enabled = 0;
		if (!lod_policy.csg_enabled)
		    lod_policy.zoom_refresh = 0;
		commit_lod_policy();
	    }
	    return BRLCAD_OK;
	}
	bu_vls_printf(gedp->ged_result_str, "Error - invalid arg: \"%s\".  Valid args are 0 or 1", argv[1]);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(argv[0], "cache")) {
	if (argc == 1) {
	    int64_t elapsedtime = bu_gettime();

	    if (!gedp || !gedp->dbip)
		return BRLCAD_ERROR;

	    struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);

	    // Clear any old cache in memory
	    bobol_mesh_lod_cache_clear_database(gedp->dbip);

	    int done = 0;
	    int total = 0;
	    struct directory *dp;
	    FOR_ALL_DIRECTORY_START(dp, gedp->dbip)
	    if (dp->d_addr == RT_DIR_PHONY_ADDR)
		continue;
	    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT)
		total++;
	    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP)
		total++;
	    FOR_ALL_DIRECTORY_END;

	    FOR_ALL_DIRECTORY_START(dp, gedp->dbip)
	    if (dp->d_addr == RT_DIR_PHONY_ADDR)
		continue;

	    // No need to open up the internal unless it's a BoT or a BRep
	    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
		done++;
		bu_log("Caching BoT %s (%d of %d)\n", dp->d_namep, done, total);
		struct BObolMeshLodCacheStatus status =
			BOBOL_MESH_LOD_CACHE_STATUS_INIT;
		if (bobol_mesh_lod_cache_refresh(gedp->dbip, dp->d_namep, &status) != BRLCAD_OK ||
		    !status.has_cache_key)
		    continue;
	    }

	    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP) {
		unsigned long long key = 0;
		struct bu_external ext = BU_EXTERNAL_INIT_ZERO;
		if (db_get_external(&ext, dp, gedp->dbip))
		    continue;
		key = bu_data_hash((void *)ext.ext_buf,  ext.ext_nbytes);
		bu_free_external(&ext);
		if (!key)
		    continue;

		struct rt_db_internal dbintern;
		RT_DB_INTERNAL_INIT(&dbintern);
		struct rt_db_internal *ip = &dbintern;
		int ret = rt_db_get_internal(ip, dp, gedp->dbip, NULL);
		if (ret < 0)
		    continue;

		if (ip->idb_minor_type != DB5_MINORTYPE_BRLCAD_BREP) {
		    rt_db_free_internal(&dbintern);
		    continue;
		}
		done++;
		struct bu_vls pname = BU_VLS_INIT_ZERO;
		bu_log("Caching BRep %s (%d of %d)\n", dp->d_namep, done, total);
		bu_vls_free(&pname);
		struct rt_brep_internal *bi = (struct rt_brep_internal *)ip->idb_ptr;
		RT_BREP_CK_MAGIC(bi);

		// Unlike a BoT, which has the mesh data already, we need to generate the
		// mesh from the brep
		int *faces = NULL;
		int face_cnt = 0;
		vect_t *normals = NULL;
		point_t *pnts = NULL;
		int pnt_cnt = 0;
		struct bn_tol *tol = &wdbp->wdb_tol;
		struct bg_tess_tol *ttol = &wdbp->wdb_ttol;

		int bret = brep_cdt_fast(&faces, &face_cnt, &normals, &pnts, &pnt_cnt, bi->brep, -1, ttol, tol);
		if (bret != BRLCAD_OK) {
		    bu_free(faces, "faces");
		    bu_free(normals, "normals");
		    bu_free(pnts, "pnts");
		    rt_db_free_internal(&dbintern);
		    continue;
		}

		// Because BRep LoD uses generated mesh data rather than a database
		// full-detail mesh payload, store it with a fidelity ratio of 1.
		struct BObolMeshLodCacheStatus status =
			BOBOL_MESH_LOD_CACHE_STATUS_INIT;
		(void)bobol_mesh_lod_cache_store_mesh(gedp->dbip, dp->d_namep,
							(const point_t *)pnts, (size_t)pnt_cnt, normals,
							faces, (size_t)face_cnt, key, 1.0, &status);

		rt_db_free_internal(&dbintern);
		bu_free(faces, "faces");
		bu_free(normals, "normals");
		bu_free(pnts, "pnts");
	    }
	    FOR_ALL_DIRECTORY_END;

	    elapsedtime = bu_gettime() - elapsedtime;
	    {
		int seconds = elapsedtime / 1000000;
		int minutes = seconds / 60;
		int hours = minutes / 60;

		minutes = minutes % 60;
		seconds = seconds %60;
		bu_vls_printf(gedp->ged_result_str, "Caching complete (Elapsed time: %02d:%02d:%02d)\n", hours, minutes, seconds);
	    }
	    return BRLCAD_OK;
	}
	if (argc == 2) {
	    if (BU_STR_EQUAL(argv[1], "clear")) {
		bobol_mesh_lod_cache_clear_database(gedp->dbip);
		(void)bobol_draw_cache_clear_database(gedp->dbip);
		return BRLCAD_OK;
	    } else if (BU_STR_EQUAL(argv[1], "exists")) {
		struct directory *dp;
		FOR_ALL_DIRECTORY_START(dp, gedp->dbip)
		if (dp->d_addr == RT_DIR_PHONY_ADDR)
		    continue;
		// checking both BoTs and BREPs
		if ((dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) ||
		    (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP)) {
		    struct BObolMeshLodCacheStatus status =
			    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
		    if (bobol_mesh_lod_cache_status(gedp->dbip, dp->d_namep, &status) != BRLCAD_OK ||
			!status.has_cache_key) {
			return BRLCAD_ERROR;
		    }
		}
		FOR_ALL_DIRECTORY_END;
		return BRLCAD_OK;
	    }
	}
	if (argc == 3) {
	    if (BU_STR_EQUAL(argv[1], "clear") && BU_STR_EQUAL(argv[2], "all_files")) {
		bobol_mesh_lod_cache_clear_all();
		bobol_draw_cache_clear_all();
		return BRLCAD_OK;
	    }
	}
	bu_vls_printf(gedp->ged_result_str, "unknown argument to cache: %s\n", argv[1]);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(argv[0], "service")) {
	if (argc == 1 || BU_STR_EQUAL(argv[1], "status")) {
	    ged_draw_obol_lod_service_status_t status = {};
	    if (!ged_draw_obol_lod_service_status(gedp, view_ctx, &status)) {
		bu_vls_printf(gedp->ged_result_str,
			      "no Obol LoD service is attached to the current view\n");
		return BRLCAD_ERROR;
	    }
	    bu_vls_printf(gedp->ged_result_str, "attached: %d\n", status.attached);
	    bu_vls_printf(gedp->ged_result_str, "running: %d\n", status.running);
	    bu_vls_printf(gedp->ged_result_str, "auto_submit: %d\n", status.auto_submit);
	    bu_vls_printf(gedp->ged_result_str, "workers: %zu\n", status.worker_count);
	    bu_vls_printf(gedp->ged_result_str, "in_flight: %zu\n", status.in_flight);
	    bu_vls_printf(gedp->ged_result_str, "pending_tasks: %zu\n", status.pending_tasks);
	    bu_vls_printf(gedp->ged_result_str, "queued_results: %zu\n", status.queued_results);
	    bu_vls_printf(gedp->ged_result_str, "queued_cache_writes: %zu\n", status.queued_cache_writes);
	    bu_vls_printf(gedp->ged_result_str, "delayed_tasks: %zu\n", status.delayed_tasks);
	    bu_vls_printf(gedp->ged_result_str, "last_visited_meshes: %u\n", status.last_visited_mesh_count);
	    bu_vls_printf(gedp->ged_result_str, "last_submitted_tasks: %u\n", status.last_submitted_task_count);
	    bu_vls_printf(gedp->ged_result_str, "last_skipped_meshes: %u\n", status.last_skipped_mesh_count);
	    bu_vls_printf(gedp->ged_result_str, "last_results: %zu\n", status.last_result_count);
	    bu_vls_printf(gedp->ged_result_str, "last_matched_results: %u\n", status.last_matched_result_count);
	    bu_vls_printf(gedp->ged_result_str, "last_applied_results: %u\n", status.last_applied_result_count);
	    bu_vls_printf(gedp->ged_result_str, "last_rejected_results: %u\n", status.last_rejected_result_count);
	    bu_vls_printf(gedp->ged_result_str, "last_unmatched_results: %u\n", status.last_unmatched_result_count);
	    bu_vls_printf(gedp->ged_result_str, "active_lod_mesh_payloads: %zu\n", status.active_mesh_payloads);
	    bu_vls_printf(gedp->ged_result_str, "active_lod_aabb_proxies: %zu\n", status.active_aabb_proxy_payloads);
	    bu_vls_printf(gedp->ged_result_str, "active_lod_obb_proxies: %zu\n", status.active_obb_proxy_payloads);
	    if (status.last_diagnostics[0])
		bu_vls_printf(gedp->ged_result_str, "last_diagnostics:\n%s\n",
			      status.last_diagnostics);
	    return BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "start")) {
	    int workers = 0;
	    if (argc == 3 &&
		(bu_opt_int(NULL, 1, (const char **)&argv[2],
			    (void *)&workers) != 1 || workers < 0)) {
		bu_vls_printf(gedp->ged_result_str,
			      "invalid worker count: %s\n", argv[2]);
		return BRLCAD_ERROR;
	    }
	    if (argc > 3) {
		bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
		return BRLCAD_ERROR;
	    }
	    if (!ged_draw_obol_lod_service_start(gedp, view_ctx,
						 (size_t)workers)) {
		bu_vls_printf(gedp->ged_result_str,
			      "unable to start Obol LoD service for the current view\n");
		return BRLCAD_ERROR;
	    }
	    redraw_view();
	    return BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "stop")) {
	    if (argc != 2) {
		bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
		return BRLCAD_ERROR;
	    }
	    if (!ged_draw_obol_lod_service_stop(gedp, view_ctx)) {
		bu_vls_printf(gedp->ged_result_str,
			      "unable to stop Obol LoD service for the current view\n");
		return BRLCAD_ERROR;
	    }
	    redraw_view();
	    return BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "poll")) {
	    int max_results = 0;
	    if (argc == 3 &&
		(bu_opt_int(NULL, 1, (const char **)&argv[2],
			    (void *)&max_results) != 1 || max_results < 0)) {
		bu_vls_printf(gedp->ged_result_str,
			      "invalid max result count: %s\n", argv[2]);
		return BRLCAD_ERROR;
	    }
	    if (argc > 3) {
		bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
		return BRLCAD_ERROR;
	    }
	    ged_draw_obol_lod_service_status_t status = {};
	    if (!ged_draw_obol_lod_service_poll(gedp, view_ctx,
						(size_t)max_results, &status)) {
		bu_vls_printf(gedp->ged_result_str,
			      "unable to poll Obol LoD service for the current view\n");
		return BRLCAD_ERROR;
	    }
	    redraw_view();
	    bu_vls_printf(gedp->ged_result_str,
			  "submitted=%u applied=%u queued=%zu in_flight=%zu pending=%zu\n",
			  status.last_submitted_task_count,
			  status.last_applied_result_count,
			  status.queued_results,
			  status.in_flight,
			  status.pending_tasks);
	    if (status.last_diagnostics[0])
		bu_vls_printf(gedp->ged_result_str, "%s\n",
			      status.last_diagnostics);
	    return BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "wait")) {
	    int timeout_ms = 5000;
	    int max_results = 0;
	    if (argc >= 3 &&
		(bu_opt_int(NULL, 1, (const char **)&argv[2],
			    (void *)&timeout_ms) != 1 || timeout_ms < 0)) {
		bu_vls_printf(gedp->ged_result_str,
			      "invalid timeout in milliseconds: %s\n", argv[2]);
		return BRLCAD_ERROR;
	    }
	    if (argc == 4 &&
		(bu_opt_int(NULL, 1, (const char **)&argv[3],
			    (void *)&max_results) != 1 || max_results < 0)) {
		bu_vls_printf(gedp->ged_result_str,
			      "invalid max result count: %s\n", argv[3]);
		return BRLCAD_ERROR;
	    }
	    if (argc > 4) {
		bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
		return BRLCAD_ERROR;
	    }

	    ged_draw_obol_lod_service_status_t status = {};
	    int timed_out = 0;
	    const int64_t start = bu_gettime();
	    const int64_t timeout_us = (int64_t)timeout_ms * 1000;

	    while (1) {
		if (!ged_draw_obol_lod_service_poll(gedp, view_ctx,
						    (size_t)max_results,
						    &status)) {
		    bu_vls_printf(gedp->ged_result_str,
				  "unable to wait on Obol LoD service for the current view\n");
		    return BRLCAD_ERROR;
		}
		if (!lod_service_status_has_work(&status))
		    break;
		if (timeout_us >= 0 && bu_gettime() - start >= timeout_us) {
		    timed_out = 1;
		    break;
		}
		bu_snooze(25000);
	    }

	    redraw_view();
	    bu_vls_printf(gedp->ged_result_str,
			  "wait_timed_out=%d submitted=%u applied=%u queued=%zu in_flight=%zu pending=%zu delayed=%zu\n",
			  timed_out,
			  status.last_submitted_task_count,
			  status.last_applied_result_count,
			  status.queued_results,
			  status.in_flight,
			  status.pending_tasks,
			  status.delayed_tasks);
	    if (status.last_diagnostics[0])
		bu_vls_printf(gedp->ged_result_str, "%s\n",
			      status.last_diagnostics);
	    return timed_out ? BRLCAD_ERROR : BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "prewarm")) {
	    ged_draw_obol_lod_service_status_t status = {};
	    size_t submitted = ged_draw_obol_lod_service_prewarm(gedp,
			       view_ctx, argc - 2, argc > 2 ? argv + 2 : NULL, &status);
	    if (submitted == 0 && !status.running) {
		bu_vls_printf(gedp->ged_result_str,
			      "unable to start Obol LoD service for cache prewarm\n");
		return BRLCAD_ERROR;
	    }
	    redraw_view();
	    bu_vls_printf(gedp->ged_result_str,
			  "prewarm_submitted=%zu queued=%zu in_flight=%zu pending=%zu\n",
			  submitted,
			  status.queued_results,
			  status.in_flight,
			  status.pending_tasks);
	    if (status.last_diagnostics[0])
		bu_vls_printf(gedp->ged_result_str, "%s\n",
			      status.last_diagnostics);
	    return BRLCAD_OK;
	}

	bu_vls_printf(gedp->ged_result_str,
		      "unknown argument to service: %s\n", argv[1]);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(argv[0], "scale")) {
	if (argc == 1) {
	    bu_vls_printf(gedp->ged_result_str, "%g\n", lod_policy.scale);
	    return BRLCAD_OK;
	}
	fastf_t scale = 1.0;
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1], (void *)&scale) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "unknown argument to point_scale: %s\n", argv[1]);
	    return BRLCAD_ERROR;
	}
	lod_policy.scale = scale;
	commit_lod_policy();
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "point_scale")) {
	if (argc == 1) {
	    bu_vls_printf(gedp->ged_result_str, "%g\n", lod_policy.point_scale);
	    return BRLCAD_OK;
	}
	fastf_t scale = 1.0;
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1], (void *)&scale) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "unknown argument to point_scale: %s\n", argv[1]);
	    return BRLCAD_ERROR;
	}
	lod_policy.point_scale = scale;
	commit_lod_policy();
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "curve_scale")) {
	if (argc == 1) {
	    bu_vls_printf(gedp->ged_result_str, "%g\n", lod_policy.curve_scale);
	    return BRLCAD_OK;
	}
	fastf_t scale = 1.0;
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1], (void *)&scale) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "unknown argument to curve_scale: %s\n", argv[1]);
	    return BRLCAD_ERROR;
	}
	lod_policy.curve_scale = scale;
	commit_lod_policy();
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "bot_threshold")) {
	if (argc == 1) {
	    bu_vls_printf(gedp->ged_result_str, "%zu\n", lod_policy.bot_threshold);
	    return BRLCAD_OK;
	}
	int bcnt = 0;
	if (bu_opt_int(NULL, 1, (const char **)&argv[1], (void *)&bcnt) != 1 || bcnt < 0) {
	    bu_vls_printf(gedp->ged_result_str, "unknown argument to bot_threshold: %s\n", argv[1]);
	    return BRLCAD_ERROR;
	}
	lod_policy.bot_threshold = (size_t)bcnt;
	commit_lod_policy();
	return BRLCAD_OK;

    }

    bu_vls_printf(gedp->ged_result_str, "unknown subcommand: %s\n", argv[0]);
    return BRLCAD_ERROR;
}



// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
