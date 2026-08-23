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
#include <cmath>
#include <limits>
#include <string.h>

#include "BObol/BDrawCache.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshLodCache.h"
#include "BObol/BViewController.h"
#include "bu/cmd.h"
#include "bu/parallel.h"
#include "bu/snooze.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "bu/vls.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "rt/view.h"

#include "../ged_bobol_private.hpp"
#include "../ged_draw_private.h"
#include "../ged_private.h"
#include "./ged_view.h"

static int
lod_valid_fps_target(fastf_t target)
{
    return std::isfinite(target) && target > 0.0 &&
	target <= static_cast<fastf_t>(std::numeric_limits<float>::max());
}

static int
lod_service_has_work(BObolViewController *controller)
{
    if (!controller)
	return 0;

    BObolLodService *service = controller->getLodService();
    return controller->hasPendingLodSubmissions() ||
	controller->hasPendingLodResults() ||
	controller->hasProgressiveWorkPending() ||
	(service && (service->inFlightCount() > 0 ||
	service->pendingTaskCountForDiagnostics() > 0 ||
	service->queuedResultCountForDiagnostics() > 0 ||
	service->queuedCacheWriteCountForDiagnostics() > 0 ||
	service->delayedTaskCountForDiagnostics() > 0));
}

static int
lod_service_poll(BObolViewController *controller, size_t max_results)
{
    if (!controller)
	return 0;
    BObolProgressiveOptions options =
	controller->getDefaultProgressiveOptions();
    options.maxLodResults = max_results;
    (void)controller->advanceProgressiveWork(&options, NULL);
    return 1;
}

static int
lod_cache_summary(struct ged *gedp, int require_complete)
{
    if (!gedp || !gedp->dbip)
	return BRLCAD_ERROR;

    struct BObolMeshLodCacheSummary summary =
	BOBOL_MESH_LOD_CACHE_SUMMARY_INIT;
    if (bobol_mesh_lod_cache_summary(gedp->dbip, &summary) != BRLCAD_OK)
	return BRLCAD_ERROR;
    bu_vls_printf(gedp->ged_result_str,
	"mapped %llu of %llu BoTs (%llu missing)\n",
	(unsigned long long)summary.mapped_bot_count,
	(unsigned long long)summary.database_bot_count,
	(unsigned long long)summary.missing_bot_count);
    return require_complete && !summary.all_bots_mapped ?
	BRLCAD_ERROR : BRLCAD_OK;
}

static int
lod_cache_command(struct ged *gedp, int argc, const char **argv)
{
    if (!gedp || !gedp->dbip || argc < 1 || !argv)
	return BRLCAD_ERROR;

    /* Cache creation belongs to the bounded LoD service.  The historical bare
     * command first erased the database cache and then regenerated every BoT
     * synchronously, freezing the UI and defeating view-aware admission. */
    if (argc == 1)
	return lod_cache_summary(gedp, 0);
    if (argc == 2) {
	if (BU_STR_EQUAL(argv[1], "clear")) {
	    bobol_mesh_lod_cache_clear_database(gedp->dbip);
	    (void)bobol_draw_cache_clear_database(gedp->dbip);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL(argv[1], "status") ||
	    BU_STR_EQUAL(argv[1], "exists"))
	    return lod_cache_summary(gedp,
		BU_STR_EQUAL(argv[1], "exists") ? 1 : 0);
    }
    if (argc == 3) {
	if (BU_STR_EQUAL(argv[1], "exists") &&
	    BU_STR_EQUAL(argv[2], "deep")) {
	    struct directory *dp;
	    FOR_ALL_DIRECTORY_START(dp, gedp->dbip)
	    if (dp->d_addr == RT_DIR_PHONY_ADDR)
		continue;
	    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
		struct BObolMeshLodCacheStatus status =
			BOBOL_MESH_LOD_CACHE_STATUS_INIT;
		if (bobol_mesh_lod_cache_status(gedp->dbip, dp->d_namep,
		    &status) != BRLCAD_OK || !status.has_cache_key ||
		    !status.has_cached_payload || status.stale_cache_entry)
		    return BRLCAD_ERROR;
	    }
	    FOR_ALL_DIRECTORY_END;
	    return BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "clear") &&
	    BU_STR_EQUAL(argv[2], "all_files")) {
	    bobol_mesh_lod_cache_clear_all();
	    bobol_draw_cache_clear_all();
	    return BRLCAD_OK;
	}
    }
    bu_vls_printf(gedp->ged_result_str, "unknown cache argument\n");
    return BRLCAD_ERROR;
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
			       "view lod cache [clear [all_files] | status | exists [deep]] \n"
			       "view lod service [status|start [workers]|stop|poll [max_results]|wait [timeout_ms] [max_results]|prewarm [all|bot ...]]\n"
			       "view lod frontier [prewarm|expand] path draw_mode max_sources max_children\n"
			       "view lod scale [factor]\n"
			       "view lod point_scale [factor]\n"
			       "view lod curve_scale [factor]\n"
			       "view lod fps [interactive [stable]]\n"
			       "view lod memory [auto|available_percent]\n"
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

    /* Cache inspection and maintenance are database operations, not view
     * operations.  Keeping them ahead of view lookup makes the cheap status
     * path usable by headless gsh and avoids constructing a graphical view
     * merely to inspect an LMDB namespace. */
    if (argc > 0 && BU_STR_EQUAL(argv[0], "cache"))
	return lod_cache_command(gedp, argc, argv);

    view_ctx = gd->cv;
    if (view_ctx == NULL) {
	bu_vls_printf(gedp->ged_result_str, "no current view defined\n");
	return BRLCAD_ERROR;
    }

    ged_view_lod_policy lod_policy;
    if (!ged_view_lod_policy_get(&lod_policy, view_ctx)) {
	bu_vls_printf(gedp->ged_result_str, "unable to read LoD policy\n");
	return BRLCAD_ERROR;
    }
    BObolViewController *view_controller =
	ged_bobol_view_controller(view_ctx);
    auto commit_lod_policy = [&]() {
	(void)ged_view_lod_policy_apply(view_ctx, &lod_policy);
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
	if (view_controller) {
	    bu_vls_printf(gedp->ged_result_str,
		"fps(interactive/stable): %g/%g\n",
		view_controller->getLodInteractiveTargetFps(),
		view_controller->getLodStableTargetFps());
	    BObolLodService *service = view_controller->getLodService();
	    if (service) {
		const double percent =
		    service->getResidentMeshAvailableMemoryPercent();
		if (percent > 0.0)
		    bu_vls_printf(gedp->ged_result_str,
			"memory(available_percent/bytes): %g/%zu\n",
			percent, service->getResidentMeshLimit());
		else
		    bu_vls_printf(gedp->ged_result_str,
			"memory(available_percent/bytes): auto/%zu\n",
			service->getResidentMeshLimit());
	    }
	    bu_vls_printf(gedp->ged_result_str,
		"scene_faces(active): %zu\n",
		view_controller->getActiveLodFaceCount());
	    bu_vls_printf(gedp->ged_result_str,
		"scene_render_cost(active/budget): %zu/%zu\n",
		view_controller->getActiveLodRenderCost(),
		view_controller->getCurrentLodRenderCostBudget());
	}
	return BRLCAD_OK;
    }

    if (BU_STR_EQUIV(argv[0], "1")) {
	if (lod_policy.policy != BV_LOD_AUTO ||
	    !lod_policy.mesh_enabled || !lod_policy.csg_enabled) {
	    lod_policy.policy = BV_LOD_AUTO;
	    lod_policy.mesh_enabled = 1;
	    lod_policy.csg_enabled = 1;
	    lod_policy.zoom_refresh = 1;
	    commit_lod_policy();
	}
	return BRLCAD_OK;
    }
    if (BU_STR_EQUIV(argv[0], "0")) {
	if (lod_policy.policy != BV_LOD_OFF ||
	    lod_policy.mesh_enabled || lod_policy.csg_enabled) {
	    lod_policy.policy = BV_LOD_OFF;
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
		lod_policy.policy = BV_LOD_AUTO;
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
		if (!lod_policy.mesh_enabled)
		    lod_policy.policy = BV_LOD_OFF;
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
		lod_policy.policy = BV_LOD_AUTO;
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
		if (!lod_policy.csg_enabled)
		    lod_policy.policy = BV_LOD_OFF;
		commit_lod_policy();
	    }
	    return BRLCAD_OK;
	}
	bu_vls_printf(gedp->ged_result_str, "Error - invalid arg: \"%s\".  Valid args are 0 or 1", argv[1]);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(argv[0], "fps")) {
	if (!view_controller) {
	    bu_vls_printf(gedp->ged_result_str,
		"no Obol LoD controller is attached to the current view\n");
	    return BRLCAD_ERROR;
	}
	if (argc == 1) {
	    bu_vls_printf(gedp->ged_result_str, "%g %g\n",
		view_controller->getLodInteractiveTargetFps(),
		view_controller->getLodStableTargetFps());
	    return BRLCAD_OK;
	}
	if (argc > 3) {
	    bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
	    return BRLCAD_ERROR;
	}
	fastf_t interactive = 0.0;
	fastf_t stable = 0.0;
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1],
		(void *)&interactive) != 1 ||
	    !lod_valid_fps_target(interactive)) {
	    bu_vls_printf(gedp->ged_result_str,
		"invalid interactive FPS target: %s\n", argv[1]);
	    return BRLCAD_ERROR;
	}
	if (argc == 3) {
	    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[2],
		    (void *)&stable) != 1 || !lod_valid_fps_target(stable)) {
		bu_vls_printf(gedp->ged_result_str,
		    "invalid stable FPS target: %s\n", argv[2]);
		return BRLCAD_ERROR;
	    }
	} else {
	    stable = interactive;
	}
	view_controller->setLodFrameRateTargets(
	    static_cast<float>(interactive), static_cast<float>(stable));
	redraw_view();
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "memory")) {
	if (!view_controller || !view_controller->getLodService()) {
	    bu_vls_printf(gedp->ged_result_str,
		"no Obol LoD service is attached to the current view\n");
	    return BRLCAD_ERROR;
	}
	if (argc > 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage:\n%s", usage);
	    return BRLCAD_ERROR;
	}
	BObolLodService *service = view_controller->getLodService();
	if (argc == 2) {
	    if (BU_STR_EQUAL(argv[1], "auto")) {
		service->setResidentMeshLimit(0);
	    } else {
		fastf_t percent = 0.0;
		if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1],
			(void *)&percent) != 1 || !std::isfinite(percent) ||
		    percent <= 0.0 || percent >
			service->getMaximumResidentMeshAvailableMemoryPercent()) {
		    bu_vls_printf(gedp->ged_result_str,
			"invalid available-memory percentage: %s "
			"(expected > 0 and <= %g)\n", argv[1],
			service->getMaximumResidentMeshAvailableMemoryPercent());
		    return BRLCAD_ERROR;
		}
		if (!service->setResidentMeshAvailableMemoryPercent(
			static_cast<double>(percent))) {
		    bu_vls_printf(gedp->ged_result_str,
			"unable to query available system memory\n");
		    return BRLCAD_ERROR;
		}
	    }
	    /* The service admission revision unlocks richer memory-limited
	     * prefixes after an increase.  A progressive/render edge is also
	     * needed for a lower ceiling to schedule stable compaction. */
	    view_controller->markProgressiveWorkPending();
	    view_controller->requestRender("lod-memory-limit");
	    redraw_view();
	}
	const double percent =
	    service->getResidentMeshAvailableMemoryPercent();
	if (percent > 0.0) {
	    bu_vls_printf(gedp->ged_result_str,
		"available_percent: %g\n"
		"available_basis_bytes: %zu\n",
		percent,
		service->getResidentMeshAvailableMemoryBasisBytes());
	} else {
	    bu_vls_printf(gedp->ged_result_str,
		"available_percent: auto\n");
	}
	bu_vls_printf(gedp->ged_result_str,
	    "resident_limit_bytes: %zu\n"
	    "maximum_available_percent: %g\n",
	    service->getResidentMeshLimit(),
	    service->getMaximumResidentMeshAvailableMemoryPercent());
	return BRLCAD_OK;
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
	    bu_vls_printf(gedp->ged_result_str, "pending_submissions: %d\n",
		view_controller->hasPendingLodSubmissions() ? 1 : 0);
	    bu_vls_printf(gedp->ged_result_str, "pending_results: %d\n",
		view_controller->hasPendingLodResults() ? 1 : 0);
	    bu_vls_printf(gedp->ged_result_str, "progressive_pending: %d\n",
		view_controller->hasProgressiveWorkPending() ? 1 : 0);
	    bu_vls_printf(gedp->ged_result_str, "workers: %zu\n",
		view_controller->getManagedLodWorkerCount());
	    bu_vls_printf(gedp->ged_result_str,
		"target_fps_interactive: %g\n",
		view_controller->getLodInteractiveTargetFps());
	    bu_vls_printf(gedp->ged_result_str,
		"target_fps_stable: %g\n",
		view_controller->getLodStableTargetFps());
	    if (service) {
		const double percent =
		    service->getResidentMeshAvailableMemoryPercent();
		if (percent > 0.0) {
		    bu_vls_printf(gedp->ged_result_str,
			"resident_memory_available_percent: %g\n"
			"resident_memory_available_basis_bytes: %zu\n",
			percent,
			service->getResidentMeshAvailableMemoryBasisBytes());
		} else {
		    bu_vls_printf(gedp->ged_result_str,
			"resident_memory_available_percent: auto\n");
		}
		bu_vls_printf(gedp->ged_result_str,
		    "resident_memory_limit_bytes: %zu\n",
		    service->getResidentMeshLimit());
	    }
	    bu_vls_printf(gedp->ged_result_str,
		"active_scene_faces: %zu\n",
		view_controller->getActiveLodFaceCount());
	    bu_vls_printf(gedp->ged_result_str,
		"active_scene_render_cost: %zu\n",
		view_controller->getActiveLodRenderCost());
	    bu_vls_printf(gedp->ged_result_str,
		"scene_render_cost_budget: %zu\n",
		view_controller->getCurrentLodRenderCostBudget());
	    bu_vls_printf(gedp->ged_result_str,
		"calibrated_render_cost_per_second: %.0f\n",
		view_controller->getCalibratedLodRenderCostPerSecond());
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
	    bu_vls_printf(gedp->ged_result_str, "last_updated_cuts: %u\n",
		view_controller->getLastLodUpdatedCutCount());
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
			  "wait_timed_out=%d submitted=%u applied=%u "
			  "pending_submissions=%d pending_results=%d "
			  "progressive_pending=%d queued=%zu in_flight=%zu "
			  "pending=%zu delayed=%zu\n",
			  timed_out,
			  view_controller ? view_controller->getLastLodSubmittedTaskCount() : 0,
			  view_controller ? view_controller->getLastLodAppliedResultCount() : 0,
			  view_controller && view_controller->hasPendingLodSubmissions() ? 1 : 0,
			  view_controller && view_controller->hasPendingLodResults() ? 1 : 0,
			  view_controller && view_controller->hasProgressiveWorkPending() ? 1 : 0,
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
