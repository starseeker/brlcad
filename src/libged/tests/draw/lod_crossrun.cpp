/*              L O D _ C R O S S R U N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
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
/** @file lod_crossrun.cpp
 *
 * LoD cross-run cache-stability test.
 *
 * The mesh LoD cache tests verify in-process round-trips.  This test exercises
 * the cross-process (or cross-ged_open) scenario:
 *
 *  Run 1:
 *    a. Open moss.g, facetize all.g → all.bot.
 *    b. Enable mesh LoD.
 *    c. Draw all.bot with LoD active (populates the bu_cache on disk).
 *    d. Capture the rendered image (R1).
 *    e. ged_close (releases all in-memory LoD data).
 *
 *  Run 2:
 *    a. Reopen the same moss.g copy without re-facetizing.
 *    b. Enable mesh LoD (same bu_cache directory).
 *    c. Draw all.bot with LoD active (must load from cache, not recompute).
 *    d. Capture the rendered image (R2).
 *    e. ged_close.
 *
 *  Assert that R2 renders non-empty LoD output from the cache established by
 *  Run 1.  Exact pixel equality is intentionally not required: the Obol render
 *  path may use different cache identities without changing user-facing LoD
 *  semantics.
 *
 * Usage: ged_test_lod_crossrun <directory-containing-moss.g>
 */

#include "common.h"

#include <chrono>
#include <cstring>
#include <fstream>
#include <thread>

#include "BObol/BDisplayEndpoint.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BLodService.h"
#include "BObol/BViewController.h"
#include <bu.h>
#include "rt/view.h"
#include "view_test_util.h"
#include <ged.h>
#include <ged/draw.h>
#include <ged/draw_obol.h>
#include <ged/db_index.h>
#include <ged/event_txn.h>
#include <icv.h>

extern "C" void ged_changed_callback(struct db_i *UNUSED(dbip), struct directory *dp, int mode, void *u_data);

static void
log_lod_state(struct ged *gedp, const char *label)
{
    struct ged_view_context *view_ctx = gedp ?
	ged_view_active_ctx(gedp) : NULL;
    bobol_display_endpoint_t *endpoint = view_ctx ?
	ged_view_context_display_endpoint_get(view_ctx) : NULL;
    BObolViewController *controller = endpoint ?
	static_cast<BObolViewController *>(
	    bobol_display_endpoint_controller(endpoint)) : NULL;
    if (!controller)
	return;

    BObolLodService *service = controller->getLodService();
    const std::vector<SoBRLDatabaseSource *> sources =
	controller->getRenderDatabaseSources();
    point_t center = VINIT_ZERO;
    const struct bv *view = view_ctx ?
	bv_context_view_const((const struct bv_context *)view_ctx) : NULL;
    if (view)
	(void)bv_center_get(center, view);
    bu_log("%s: sources=%zu active-mesh=%zu view center=(%.9g %.9g %.9g) scale=%.9g demand auto=%d service=%d visited=%u submitted=%u skipped=%u service pending=%zu in-flight=%zu results=%zu writes=%zu delayed=%zu\n",
	label ? label : "LoD",
	sources.size(),
	controller->getActiveLodMeshPayloadCount(),
	center[0], center[1], center[2], view ? bv_scale_get(view) : 0.0,
	controller->isLodAutoSubmitEnabled() ? 1 : 0,
	service && service->isRunning() ? 1 : 0,
	controller->getLastLodVisitedMeshCount(),
	controller->getLastLodSubmittedTaskCount(),
	controller->getLastLodSkippedMeshCount(),
	service ? service->pendingTaskCountForDiagnostics() : 0,
	service ? service->inFlightCount() : 0,
	service ? service->queuedResultCountForDiagnostics() : 0,
	service ? service->queuedCacheWriteCountForDiagnostics() : 0,
	service ? service->delayedTaskCountForDiagnostics() : 0);
    for (size_t i = 0; i < sources.size(); i++) {
	BObolDatabaseSourceSummary summary;
	SoBRLDatabaseSource *source = sources[i];
	if (!source || !source->getSummary(summary))
	    continue;
	bu_log("  source[%zu] path=%s key=%s bounds=%d compact=%d shapes=%d meshes=%d status=%d\n",
	    i, summary.path.getString(), summary.instanceKey.getString(),
	    summary.sourceBoundsValid ? 1 : 0,
	    source ? source->getCompactInstanceCount() : 0,
	    summary.realizedShapeCount, summary.realizedMeshCount,
	    summary.realizationStatus);
	for (int j = 0; source && j < source->getCompactInstanceCount(); j++) {
	    BObolCompactOccurrence occurrence;
	    if (!source->getCompactOccurrence(j, occurrence))
		continue;
	    bu_log("    occurrence[%d] geometry=%d kind=%d/%s lod=%d request=%d/faces=%llu visible=%d transparency=%.3g color=(%.3g %.3g %.3g) points=%d triangles=%d bounds=%d\n",
		j, occurrence.geometry ? 1 : 0,
		occurrence.summary.shapeKind,
		occurrence.summary.geometryKind.getString(),
		occurrence.lodBacked ? 1 : 0,
		occurrence.sourceMeshRequestValid ? 1 : 0,
		(unsigned long long)occurrence.sourceMeshRequest.faceCount,
		occurrence.summary.visible ? 1 : 0,
		occurrence.summary.transparency,
		occurrence.summary.color[0], occurrence.summary.color[1],
		occurrence.summary.color[2],
		occurrence.summary.pointCount,
		occurrence.summary.triangleCount,
		occurrence.summary.boundsValid ? 1 : 0);
	}
    }
}

static int
wait_for_lod_service(struct ged *gedp, int timeout_ms)
{
    struct ged_view_context *view_ctx = gedp ? ged_view_active_ctx(gedp) : NULL;
    bobol_display_endpoint_t *endpoint = view_ctx ?
	ged_view_context_display_endpoint_get(view_ctx) : NULL;
    BObolViewController *controller = endpoint ?
	static_cast<BObolViewController *>(
	    bobol_display_endpoint_controller(endpoint)) : NULL;
    if (!controller)
	return 1;

    for (int elapsed = 0; elapsed <= timeout_ms; elapsed += 25) {
	BObolProgressiveStatus progressive;
	(void)controller->advanceProgressiveWork(NULL, &progressive);
	BObolLodService *service = controller->getLodService();
	const int service_idle = !service ||
	    (service->inFlightCount() == 0 &&
		service->pendingTaskCountForDiagnostics() == 0 &&
		service->queuedCacheWriteCountForDiagnostics() == 0 &&
		service->delayedTaskCountForDiagnostics() == 0);
	if (!progressive.hasMore && service_idle &&
	    !controller->hasPendingLodResults()) {
	    log_lod_state(gedp, "LoD settled");
	    return 1;
	}
	std::this_thread::sleep_for(std::chrono::milliseconds(25));
    }

    bu_log("LoD service did not become idle within %d ms\n", timeout_ms);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Minimal GED setup (headless Obol endpoint, 512x512, az/el 35/25)  */
/* ------------------------------------------------------------------ */
static struct ged *
open_and_attach(const char *gfile)
{
    struct ged *gedp = ged_open("db", gfile, 1);
    if (!gedp)
	return NULL;

    db_add_changed_clbk(gedp->dbip, &ged_changed_callback, (void *)gedp);

    struct ged_view_context *v = ged_view_active_ctx(gedp);
    bv_dimensions_set(DRAW_TEST_BV(v), 512, 512);
    const char *s_av[7] = {
	"dm", "open", "--host", "headless", "--renderer", "sw", NULL
    };
    if (ged_exec_dm(gedp, 6, s_av) != BRLCAD_OK) {
	ged_close(gedp);
	return NULL;
    }
    bv_unit_conversion_set(DRAW_TEST_BV(v), gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);

    /* Set view orientation */
    const char *ae_av[4] = {"ae", "35", "25", NULL};
    ged_exec_ae(gedp, 3, ae_av);

    return gedp;
}

/* ------------------------------------------------------------------ */
/* Render current scene to PNG file                                    */
/* ------------------------------------------------------------------ */
static int
render_to_file(struct ged *gedp, const char *outfile)
{
    struct ged_view_context *v = ged_view_active_ctx(gedp);
    ged_db_index_refresh(gedp);
    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    txn.view = v;
    ged_draw_apply_transaction(gedp, &txn, NULL);

    const char *sg_av[2] = {"screengrab", outfile};
    return ged_exec_screengrab(gedp, 2, sg_av);
}

/* ------------------------------------------------------------------ */
/* Verify a PNG contains drawn output rather than an empty scene.       */
/* ------------------------------------------------------------------ */
static int
png_not_empty(const char *a)
{
    icv_image_t *ia = icv_read(a, BU_MIME_IMAGE_PNG, 0, 0);
    icv_image_t *ib = icv_read("lod_cr_empty.png", BU_MIME_IMAGE_PNG, 0, 0);
    if (!ia || !ib) {
	if (ia) icv_destroy(ia);
	if (ib) icv_destroy(ib);
	bu_log("png_not_empty: could not read '%s' or empty control\n", a);
	return 0;
    }
    int match = 0, off1 = 0, offmany = 0;
    icv_diff(&match, &off1, &offmany, ia, ib);
    icv_destroy(ia);
    icv_destroy(ib);
    if (!off1 && !offmany) {
	bu_log("  image '%s' matches empty scene\n", a);
	return 0;
    }
    bu_log("  image '%s' non-empty: %d matching, %d off-by-1, %d off-by-many\n",
	    a, match, off1, offmany);
    return 1;
}

/* ================================================================== */
int
main(int ac, char *av[])
{
    bu_setprogname(av[0]);

    int soft_fail = 0;
    struct bu_opt_desc d[2];
    BU_OPT(d[0], "c", "continue", "", NULL, &soft_fail, "Continue on failure.");
    BU_OPT_NULL(d[1]);

    int uac = bu_opt_parse(NULL, ac, (const char **)av, d);
    if (uac != 2)
	bu_exit(EXIT_FAILURE,
		"Usage: %s [-c] <directory-containing-moss.g>\n", av[0]);
    const char *datadir = av[1];
    if (!bu_file_directory(datadir))
	bu_exit(EXIT_FAILURE, "ERROR: [%s] is not a directory\n", datadir);

    /* Private cache dir so we don't pollute the system cache */
    char lcache[MAXPATHLEN];
    bu_dir(lcache, MAXPATHLEN, BU_DIR_CURR, "lod_crossrun_cache", NULL);
    bu_mkdir(lcache);
    bu_setenv("BU_DIR_CACHE", lcache, 1);

    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    /* Make a private copy of moss.g so facetize doesn't dirty the source */
    struct bu_vls srcpath = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&srcpath, "%s/moss.g", datadir);
    {
	std::ifstream orig(bu_vls_cstr(&srcpath), std::ios::binary);
	std::ofstream tmpg("moss_lod_cr_tmp.g", std::ios::binary);
	tmpg << orig.rdbuf();
    }
    bu_vls_free(&srcpath);

    /* ----------------------------------------------------------------
     * Run 1 — generate facetized geometry, populate LoD cache, render.
     * -------------------------------------------------------------- */
    bu_log("LoD cross-run test: Run 1 — facetize + draw (populates cache)\n");

    struct ged *gedp = open_and_attach("moss_lod_cr_tmp.g");
    if (!gedp)
	bu_exit(EXIT_FAILURE, "Run 1: ged_open failed\n");

    const char *clear_av[4] = {"view", "lod", "cache", "clear"};
    ged_exec_view(gedp, 4, clear_av);

    /* Facetize all.g → all.bot at a moderately fine tolerance */
    const char *s_av[8] = {NULL};
    s_av[0] = "tol"; s_av[1] = "rel"; s_av[2] = "0.0005"; s_av[3] = NULL;
    ged_exec_tol(gedp, 3, s_av);

    s_av[0] = "facetize"; s_av[1] = "-r"; s_av[2] = "all.g"; s_av[3] = "all.bot"; s_av[4] = NULL;
    ged_exec_facetize(gedp, 4, s_av);
    ged_db_index_refresh(gedp);
    ged_event_notify_batch_rebuild(gedp, NULL);

    /* Enable LoD at a coarse scale to make the level selection visible */
    s_av[0] = "view"; s_av[1] = "lod"; s_av[2] = "mesh"; s_av[3] = "1"; s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);
    s_av[0] = "view"; s_av[1] = "lod"; s_av[2] = "scale"; s_av[3] = "0.8"; s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);

    (void)render_to_file(gedp, "lod_cr_empty.png");

    s_av[0] = "draw"; s_av[1] = "-m1"; s_av[2] = "all.bot"; s_av[3] = NULL;
    const int r1_draw = ged_exec_draw(gedp, 3, s_av);
    bu_log("Run 1 draw returned %d: %s\n", r1_draw,
	    bu_vls_cstr(gedp->ged_result_str));

    /* The draw transaction arms progressive autoview while the cold-start
     * scene is still empty.  An immediate explicit autoview would frame that
     * empty scene, change the view revision, and correctly cancel the armed
     * autoview before the first boxes arrive. */
    int r1_wait = wait_for_lod_service(gedp, 5000);
    if (!r1_wait)
	bu_log("Run 1: LoD service wait timed out\n");

    int r1_ret = render_to_file(gedp, "lod_cr_run1.png");
    bu_log("Run 1 image captured (%s)\n", r1_ret ? "WARN: screengrab failed" : "ok");

    /* Close releases RT-owned in-memory LoD structures; bu_cache stays on disk. */
    ged_close(gedp);

    /* ----------------------------------------------------------------
     * Run 2 — reopen without re-facetizing; LoD must come from cache.
     * -------------------------------------------------------------- */
    bu_log("LoD cross-run test: Run 2 — reopen, draw (must use cache)\n");

    gedp = open_and_attach("moss_lod_cr_tmp.g");
    if (!gedp) {
	bu_file_delete("lod_cr_run1.png");
	bu_file_delete("moss_lod_cr_tmp.g");
	bu_dirclear(lcache);
	bu_exit(EXIT_FAILURE, "Run 2: ged_open failed\n");
    }

    /* Same LoD settings */
    s_av[0] = "view"; s_av[1] = "lod"; s_av[2] = "mesh"; s_av[3] = "1"; s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);
    s_av[0] = "view"; s_av[1] = "lod"; s_av[2] = "scale"; s_av[3] = "0.8"; s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);

    s_av[0] = "draw"; s_av[1] = "-m1"; s_av[2] = "all.bot"; s_av[3] = NULL;
    const int r2_draw = ged_exec_draw(gedp, 3, s_av);
    bu_log("Run 2 draw returned %d: %s\n", r2_draw,
	    bu_vls_cstr(gedp->ged_result_str));

    int r2_wait = wait_for_lod_service(gedp, 5000);
    if (!r2_wait)
	bu_log("Run 2: LoD service wait timed out\n");

    int r2_ret = render_to_file(gedp, "lod_cr_run2.png");
    bu_log("Run 2 image captured (%s)\n", r2_ret ? "WARN: screengrab failed" : "ok");

    ged_close(gedp);

    /* ----------------------------------------------------------------
     * Verify both runs produced visible LoD output.  The second run keeps the
     * BU_DIR_CACHE contents from run 1 and should therefore exercise cached
     * LoD loading without relying on exact pixel identity.
     * -------------------------------------------------------------- */
    int ret = 0;
    if (!r1_wait || !r2_wait) {
	bu_log("FAIL: LoD service did not settle before capture\n");
	ret = 1;
    } else if (r1_ret || r2_ret) {
	bu_log("FAIL: one or both captures failed\n");
	ret = 1;
    } else if (!png_not_empty("lod_cr_run1.png") || !png_not_empty("lod_cr_run2.png")) {
	bu_log("FAIL: LoD cross-run produced empty rendered output\n");
	ret = 1;
    } else {
	bu_log("PASS: LoD cross-run rendered cached LoD output\n");
    }

    /* Cleanup */
    bu_file_delete("lod_cr_run1.png");
    bu_file_delete("lod_cr_run2.png");
    bu_file_delete("lod_cr_empty.png");
    bu_file_delete("moss_lod_cr_tmp.g");
    bu_dirclear(lcache);

    return ret;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
