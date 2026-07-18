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
#include "BObol/BLodService.h"
#include "BObol/BMeshLodCache.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "bu/cmd.h"
#include "bu/hash.h"
#include "bu/parallel.h"
#include "bu/snooze.h"
#include "bu/str.h"
#include "bu/time.h"
#include "bu/vls.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "rt/view.h"

#include "../ged_bobol_private.hpp"
#include "../ged_draw_private.h"
#include "../ged_draw_view_private.h"
#include "../ged_private.h"
#include "./ged_view.h"

static int
lod_service_has_work(BObolViewController *controller)
{
    BObolLodService *service = controller ? controller->getLodService() : NULL;
    if (!service)
	return 0;

    return service->inFlightCount() > 0 ||
	service->pendingTaskCountForDiagnostics() > 0 ||
	service->queuedResultCountForDiagnostics() > 0 ||
	service->queuedCacheWriteCountForDiagnostics() > 0 ||
	service->delayedTaskCountForDiagnostics() > 0;
}

static int
lod_service_poll(BObolViewController *controller, size_t max_results)
{
    if (!controller)
	return 0;
    if (controller->hasPendingLodResults())
	(void)controller->processPendingLodResults(max_results);
    if (controller->isLodAutoSubmitEnabled())
	(void)controller->submitLodRequestsIfNeeded();
    return 1;
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
    struct ged_view_context *view_ctx;
    int print_help = 0;
    static const char *usage = "view lod [csg|mesh] [0|1]\n"
			       "view lod cache [clear [all_files] | exists] \n"
			       "view lod service [status|start [workers]|stop|poll [max_results]|wait [timeout_ms] [max_results]|prewarm [all|bot ...]]\n"
			       "view lod frontier [prewarm|expand] path draw_mode max_sources max_children\n"
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

    if (argc > 6 ||
	(argc > 4 && !(argc == 6 && BU_STR_EQUAL(argv[0], "frontier"))) ||
	(argc > 3 && !(argc >= 2 && BU_STR_EQUAL(argv[0], "service") &&
		       (BU_STR_EQUAL(argv[1], "prewarm") ||
			BU_STR_EQUAL(argv[1], "wait"))) &&
	 !(argc == 6 && BU_STR_EQUAL(argv[0], "frontier")))) {
	bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
	return BRLCAD_ERROR;
    }

    view_ctx = gd->cv;
    if (view_ctx == NULL) {
	bu_vls_printf(gedp->ged_result_str, "no current view defined\n");
	return BRLCAD_ERROR;
    }

    BObolViewAttachment *attachment =
	ged_view_context_obol_attachment(view_ctx);
    if (!attachment) {
	bu_vls_printf(gedp->ged_result_str, "unable to read LoD policy\n");
	return BRLCAD_ERROR;
    }
    ged_draw_source_lod_policy lod_policy;
    attachment->getLodPolicy(&lod_policy);
    BObolViewController *view_controller =
	ged_bobol_view_controller(view_ctx);
    auto commit_lod_policy = [&]() {
	attachment->setLodPolicy(&lod_policy);
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
	    if (!view_controller) {
		bu_vls_printf(gedp->ged_result_str,
			      "no Obol LoD service is attached to the current view\n");
		return BRLCAD_ERROR;
	    }
	    BObolLodService *service = view_controller->getLodService();
	    const SbString &diagnostics = view_controller->getLastLodDiagnostics();
	    bu_vls_printf(gedp->ged_result_str, "attached: 1\n");
	    bu_vls_printf(gedp->ged_result_str, "running: %d\n",
		service && service->isRunning() ? 1 : 0);
	    bu_vls_printf(gedp->ged_result_str, "auto_submit: %d\n",
		view_controller->isLodAutoSubmitEnabled() ? 1 : 0);
	    bu_vls_printf(gedp->ged_result_str, "workers: %zu\n",
		view_controller->getManagedLodWorkerCount());
	    bu_vls_printf(gedp->ged_result_str, "in_flight: %zu\n",
		service ? service->inFlightCount() : 0);
	    bu_vls_printf(gedp->ged_result_str, "pending_tasks: %zu\n",
		service ? service->pendingTaskCountForDiagnostics() : 0);
	    bu_vls_printf(gedp->ged_result_str, "queued_results: %zu\n",
		service ? service->queuedResultCountForDiagnostics() : 0);
	    bu_vls_printf(gedp->ged_result_str, "queued_cache_writes: %zu\n",
		service ? service->queuedCacheWriteCountForDiagnostics() : 0);
	    bu_vls_printf(gedp->ged_result_str, "delayed_tasks: %zu\n",
		service ? service->delayedTaskCountForDiagnostics() : 0);
	    bu_vls_printf(gedp->ged_result_str, "last_visited_meshes: %u\n",
		view_controller->getLastLodVisitedMeshCount());
	    bu_vls_printf(gedp->ged_result_str, "last_submitted_tasks: %u\n",
		view_controller->getLastLodSubmittedTaskCount());
	    bu_vls_printf(gedp->ged_result_str, "last_skipped_meshes: %u\n",
		view_controller->getLastLodSkippedMeshCount());
	    bu_vls_printf(gedp->ged_result_str, "last_results: %zu\n",
		view_controller->getLastLodResultCount());
	    bu_vls_printf(gedp->ged_result_str, "last_matched_results: %u\n",
		view_controller->getLastLodMatchedResultCount());
	    bu_vls_printf(gedp->ged_result_str, "last_applied_results: %u\n",
		view_controller->getLastLodAppliedResultCount());
	    bu_vls_printf(gedp->ged_result_str, "last_rejected_results: %u\n",
		view_controller->getLastLodRejectedResultCount());
	    bu_vls_printf(gedp->ged_result_str, "last_unmatched_results: %u\n",
		view_controller->getLastLodUnmatchedResultCount());
	    bu_vls_printf(gedp->ged_result_str, "active_lod_mesh_payloads: %zu\n",
		view_controller->getActiveLodMeshPayloadCount());
	    bu_vls_printf(gedp->ged_result_str, "active_lod_aabb_proxies: %zu\n",
		view_controller->getActiveLodProxyPayloadCount(BOBOL_LOD_PROXY_AABB));
	    bu_vls_printf(gedp->ged_result_str, "active_lod_obb_proxies: %zu\n",
		view_controller->getActiveLodProxyPayloadCount(BOBOL_LOD_PROXY_OBB));
	    if (diagnostics.getLength() > 0)
		bu_vls_printf(gedp->ged_result_str, "last_diagnostics:\n%s\n",
			      diagnostics.getString());
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
	    size_t worker_count = workers > 0 ? (size_t)workers : bu_avail_cpus();
	    if (!view_controller ||
		!view_controller->ensureManagedLodService(
		    worker_count ? worker_count : 1)) {
		bu_vls_printf(gedp->ged_result_str,
			      "unable to start Obol LoD service for the current view\n");
		return BRLCAD_ERROR;
	    }
	    view_controller->setLodAutoSubmit(TRUE);
	    view_controller->requestRender("lod-service-start");
	    redraw_view();
	    return BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "stop")) {
	    if (argc != 2) {
		bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
		return BRLCAD_ERROR;
	    }
	    if (!view_controller) {
		bu_vls_printf(gedp->ged_result_str,
			      "unable to stop Obol LoD service for the current view\n");
		return BRLCAD_ERROR;
	    }
	    view_controller->stopManagedLodService();
	    view_controller->requestRender("lod-service-stop");
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
	    if (!lod_service_poll(view_controller, (size_t)max_results)) {
		bu_vls_printf(gedp->ged_result_str,
			      "unable to poll Obol LoD service for the current view\n");
		return BRLCAD_ERROR;
	    }
	    BObolLodService *service = view_controller->getLodService();
	    const SbString &diagnostics = view_controller->getLastLodDiagnostics();
	    redraw_view();
	    bu_vls_printf(gedp->ged_result_str,
			  "submitted=%u applied=%u queued=%zu in_flight=%zu pending=%zu\n",
			  view_controller->getLastLodSubmittedTaskCount(),
			  view_controller->getLastLodAppliedResultCount(),
			  service ? service->queuedResultCountForDiagnostics() : 0,
			  service ? service->inFlightCount() : 0,
			  service ? service->pendingTaskCountForDiagnostics() : 0);
	    if (diagnostics.getLength() > 0)
		bu_vls_printf(gedp->ged_result_str, "%s\n",
			      diagnostics.getString());
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

	    int timed_out = 0;
	    const int64_t start = bu_gettime();
	    const int64_t timeout_us = (int64_t)timeout_ms * 1000;

	    while (1) {
		if (!lod_service_poll(view_controller, (size_t)max_results)) {
		    bu_vls_printf(gedp->ged_result_str,
				  "unable to wait on Obol LoD service for the current view\n");
		    return BRLCAD_ERROR;
		}
		if (!lod_service_has_work(view_controller))
		    break;
		if (timeout_us >= 0 && bu_gettime() - start >= timeout_us) {
		    timed_out = 1;
		    break;
		}
		bu_snooze(25000);
	    }

	    BObolLodService *service = view_controller ?
		view_controller->getLodService() : NULL;
	    const SbString diagnostics = view_controller ?
		view_controller->getLastLodDiagnostics() : SbString();
	    redraw_view();
	    bu_vls_printf(gedp->ged_result_str,
			  "wait_timed_out=%d submitted=%u applied=%u queued=%zu in_flight=%zu pending=%zu delayed=%zu\n",
			  timed_out,
			  view_controller ? view_controller->getLastLodSubmittedTaskCount() : 0,
			  view_controller ? view_controller->getLastLodAppliedResultCount() : 0,
			  service ? service->queuedResultCountForDiagnostics() : 0,
			  service ? service->inFlightCount() : 0,
			  service ? service->pendingTaskCountForDiagnostics() : 0,
			  service ? service->delayedTaskCountForDiagnostics() : 0);
	    if (diagnostics.getLength() > 0)
		bu_vls_printf(gedp->ged_result_str, "%s\n",
			      diagnostics.getString());
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

    /* Frontier realization is an adapter detail.  Keep it behind the command
     * boundary so diagnostics consume semantic paths and stable key/value
     * results rather than installed Obol records or scene-node operations. */
    if (BU_STR_EQUAL(argv[0], "frontier")) {
	if (argc != 6 ||
	    (!BU_STR_EQUAL(argv[1], "prewarm") &&
	     !BU_STR_EQUAL(argv[1], "expand"))) {
	    bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
	    return BRLCAD_ERROR;
	}

	int draw_mode = 0;
	int max_sources = 0;
	int max_children = 0;
	if (bu_opt_int(NULL, 1, (const char **)&argv[3], &draw_mode) != 1 ||
	    bu_opt_int(NULL, 1, (const char **)&argv[4], &max_sources) != 1 ||
	    bu_opt_int(NULL, 1, (const char **)&argv[5], &max_children) != 1 ||
	    max_sources < 0 || max_children < 0) {
	    bu_vls_printf(gedp->ged_result_str,
		"invalid frontier limits or draw mode\n");
	    return BRLCAD_ERROR;
	}

	if (BU_STR_EQUAL(argv[1], "prewarm")) {
	    struct ged_draw_obol_source_prewarm_status status = {};
	    const size_t submitted =
		ged_draw_obol_database_source_prewarm_visible_child_aabb_proxies(
		    gedp, view_ctx, argv[2], draw_mode,
		    (size_t)max_sources, (size_t)max_children, &status);
	    bu_vls_printf(gedp->ged_result_str,
		"result=%zu child_count=%zu considered=%zu submitted=%zu "
		"already_cached=%zu skipped_non_union=%zu "
		"skipped_duplicate_instance=%zu shared_request=%zu "
		"non_union_children=%zu duplicate_instances=%zu "
		"skipped_invalid=%zu remaining=%zu comb_sources=%zu "
		"leaf_sources=%zu\n",
		submitted, status.child_count, status.considered,
		status.submitted, status.already_cached,
		status.skipped_non_union, status.skipped_duplicate_instance,
		status.shared_request, status.non_union_children,
		status.duplicate_instances, status.skipped_invalid,
		status.remaining, status.comb_sources, status.leaf_sources);
	    return BRLCAD_OK;
	}

	struct ged_draw_obol_source_expansion_status status = {};
	const int changed =
	    ged_draw_obol_database_source_expand_visible_children(
		gedp, view_ctx, argv[2], draw_mode, (size_t)max_sources,
		(size_t)max_children, &status);
	bu_vls_printf(gedp->ged_result_str,
	    "result=%d child_count=%zu considered=%zu expanded=%zu "
	    "existing=%zu skipped_non_union=%zu "
	    "skipped_duplicate_instance=%zu expanded_non_union=%zu "
	    "expanded_duplicate_instance=%zu skipped_invalid=%zu "
	    "remaining=%zu proxy_published=%zu metadata_applied=%zu "
	    "comb_sources=%zu leaf_sources=%zu\n",
	    changed, status.child_count, status.considered, status.expanded,
	    status.existing, status.skipped_non_union,
	    status.skipped_duplicate_instance, status.expanded_non_union,
	    status.expanded_duplicate_instance, status.skipped_invalid,
	    status.remaining, status.proxy_published, status.metadata_applied,
	    status.comb_sources, status.leaf_sources);
	return BRLCAD_OK;
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
