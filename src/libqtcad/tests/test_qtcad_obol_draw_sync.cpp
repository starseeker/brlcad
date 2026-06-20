/*          T E S T _ Q T C A D _ O B O L _ D R A W _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/view_controller.h"
#include "bsg/util.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "qtcad/QgObolDatabaseSync.h"
#include "qtcad/QgObolDrawSync.h"
#include "qtcad/QgView.h"
#include "raytrace.h"
#include "rt/view_legacy_bsg.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoCamera.h>

#include <QApplication>
#include <QImage>

#include <math.h>
#include <stdio.h>
#include <string.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static int
lit_pixel_count(const QImage &image)
{
    QImage rgba = image.convertToFormat(QImage::Format_RGBA8888);
    int count = 0;
    for (int y = 0; y < rgba.height(); y++) {
	const unsigned char *line = rgba.constScanLine(y);
	for (int x = 0; x < rgba.width(); x++) {
	    const unsigned char *p = line + x * 4;
	    if (p[0] > 32 || p[1] > 32 || p[2] > 32)
		count++;
	}
    }
    return count;
}

static int
camera_positions_differ(const SbVec3f &a, const SbVec3f &b, float min_delta)
{
    return fabsf(a[0] - b[0]) > min_delta ||
	fabsf(a[1] - b[1]) > min_delta ||
	fabsf(a[2] - b[2]) > min_delta;
}

static int
make_draw_sync_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    point_t center = {5.0, 0.0, 0.0};

    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0 &&
	mk_sph(wdbp, "ball.s", center, 1.0) == 0;
    wdb_close(wdbp);
    return ret;
}

static SoBRLDatabaseSource *
source_for_path(BRLObolViewController *controller, const char *path)
{
    if (!controller || !path)
	return NULL;
    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (source && BU_STR_EQUAL(source->path.getValue().getString(), path))
	    return source;
    }
    return NULL;
}

struct draw_observer_sync_context {
    QgView *view;
    int calls;
    int changed;
};

static void
test_draw_observer(struct ged *gedp,
	const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result,
	void *client_data)
{
    struct draw_observer_sync_context *ctx =
	static_cast<struct draw_observer_sync_context *>(client_data);
    if (!ctx || !ctx->view)
	return;
    ctx->calls++;
    ctx->changed += qg_obol_sync_draw_transaction(gedp, txn, result,
	    ctx->view);
}

static int
apply_and_sync(struct ged *gedp,
	QgView *view,
	struct ged_draw_transaction *txn,
	int expect_changed)
{
    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int draw_ret = ged_draw_apply_transaction(gedp, txn, &result);
    int changed = qg_obol_sync_draw_transaction(gedp, txn, &result, view);
    ged_draw_transaction_result_free(&result);

    if (draw_ret < 0)
	return 0;
    return expect_changed ? changed != 0 : changed == 0;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_draw_sync_tmp.g";
    if (!make_draw_sync_db(dbpath))
	FAIL("failed to create qtcad Obol draw-sync test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol draw-sync test database");

    QgView view(NULL, QgView_SW, NULL);
    view.resize(180, 140);
    gedp->ged_gvp = view.view();

    BRLObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();

    struct draw_observer_sync_context obs = {&view, 0, 0};
    ged_draw_observer_token observerToken =
	ged_draw_observer_add(gedp, test_draw_observer, &obs);
    if (!observerToken)
	FAIL("GED draw observer should register for qtcad Obol draw sync");

    const char *draw_cmd[2] = {"draw", "box.s"};
    if (ged_exec_draw(gedp, 2, draw_cmd) != BRLCAD_OK)
	FAIL("real GED draw command should succeed for observer sync");
    if (obs.calls <= 0 || obs.changed <= 0)
	FAIL("real GED draw command should notify and sync qtcad Obol");
    if (controller->getDatabaseSourceCount() != 1)
	FAIL("observer-synced GED draw should create one Obol database source");
    SoBRLDatabaseSource *observerSource = source_for_path(controller, "box.s");
    if (!observerSource ||
	    observerSource->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    observerSource->getRealizedShapeCount() <= 0)
	FAIL("observer-synced GED draw should realize Obol wire geometry");

    controller->getViewport()->viewAll();
    controller->requestRender("observer-draw-visible");
    QImage observerImage;
    view.get_viewport_image(observerImage);
    if (observerImage.isNull() || lit_pixel_count(observerImage) < 10)
	FAIL("observer-synced GED draw should be visible through qtcad capture");

    SoCamera *camera = controller->getCamera();
    if (!camera)
	FAIL("qtcad Obol controller should expose a camera for view sync");
    point_t offcenter = {100.0, 0.0, 0.0};
    rt_view_center_vec_set_bsg(view.view(), offcenter);
    rt_view_scale_set_bsg(view.view(), 250.0);
    view.need_update(QG_VIEW_REFRESH);
    SbVec3f offTargetCamera = camera->position.getValue();
    if (offTargetCamera[0] < 50.0f)
	FAIL("qtcad refresh should sync Obol camera from BSG view state");

    const char *autoview_cmd[1] = {"autoview"};
    if (ged_exec_autoview(gedp, 1, autoview_cmd) != BRLCAD_OK)
	FAIL("real GED autoview command should succeed for qtcad Obol view sync");
    view.need_update(QG_VIEW_REFRESH);
    SbVec3f autoviewCamera = camera->position.getValue();
    if (!camera_positions_differ(offTargetCamera, autoviewCamera, 10.0f))
	FAIL("GED autoview should update qtcad Obol camera through view refresh");
    if (fabsf(autoviewCamera[0]) > 25.0f)
	FAIL("GED autoview should recenter qtcad Obol camera near drawn geometry");
    controller->requestRender("observer-autoview-visible");
    QImage autoviewImage;
    view.get_viewport_image(autoviewImage);
    if (autoviewImage.isNull() || lit_pixel_count(autoviewImage) < 10)
	FAIL("GED-autoviewed Obol scene should stay visible through qtcad capture");

    obs.calls = 0;
    obs.changed = 0;
    const char *erase_cmd[2] = {"erase", "box.s"};
    if (ged_exec_erase(gedp, 2, erase_cmd) != BRLCAD_OK)
	FAIL("real GED erase command should succeed for observer sync");
    if (obs.calls <= 0 || obs.changed <= 0)
	FAIL("real GED erase command should notify and sync qtcad Obol");
    if (controller->getDatabaseSourceCount() != 0)
	FAIL("observer-synced GED erase should remove Obol database sources");
    if (ged_draw_observer_remove(gedp, observerToken) != 1)
	FAIL("GED draw observer should unregister after qtcad Obol sync test");
    controller->clearDatabaseSources();

    struct ged_draw_transaction draw_box =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "box.s");
    draw_box.view = view.view();
    if (!apply_and_sync(gedp, &view, &draw_box, 1))
	FAIL("GED draw should sync a wire Obol database source");
    if (controller->getDatabaseSourceCount() != 1)
	FAIL("Obol draw sync should create one database source");
    SoBRLDatabaseSource *source = controller->getDatabaseSource(0);
    if (!source ||
	    !BU_STR_EQUAL(source->path.getValue().getString(), "box.s") ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    source->getRealizedShapeCount() <= 0)
	FAIL("wire Obol database source should preserve path and realized geometry");

    controller->getViewport()->viewAll();
    controller->requestRender("draw-sync-visible");
    QImage visibleImage;
    view.get_viewport_image(visibleImage);
    if (visibleImage.isNull() || lit_pixel_count(visibleImage) < 10)
	FAIL("Obol-synced GED draw should be visible through qtcad capture");

    struct ged_draw_transaction erase_box =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "box.s");
    erase_box.view = view.view();
    if (!apply_and_sync(gedp, &view, &erase_box, 1))
	FAIL("GED erase should remove an Obol database source");
    if (controller->getDatabaseSourceCount() != 0)
	FAIL("Obol draw sync should remove erased database sources");

    struct bsg_appearance_settings shaded_settings = BSG_APPEARANCE_SETTINGS_INIT;
    shaded_settings.draw_mode = BSG_DRAW_MODE_SHADED;
    struct ged_draw_transaction draw_ball =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "ball.s");
    draw_ball.view = view.view();
    draw_ball.appearance = &shaded_settings;
    if (!apply_and_sync(gedp, &view, &draw_ball, 1))
	FAIL("GED shaded draw should sync a shaded Obol database source");
    source = controller->getDatabaseSource(0);
    if (!source ||
	    !BU_STR_EQUAL(source->path.getValue().getString(), "ball.s") ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::SHADED ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    source->getRealizedMeshCount() <= 0)
	FAIL("shaded Obol database source should preserve draw mode and mesh geometry");

    const char *paths[2] = {"box.s", "ball.s"};
    struct bsg_appearance_settings wire_settings = BSG_APPEARANCE_SETTINGS_INIT;
    struct ged_draw_transaction draw_both =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, NULL);
    draw_both.view = view.view();
    draw_both.paths = paths;
    draw_both.path_count = 2;
    draw_both.appearance = &wire_settings;
    if (!apply_and_sync(gedp, &view, &draw_both, 1))
	FAIL("multi-path GED draw should sync multiple Obol database sources");
    if (controller->getDatabaseSourceCount() != 2)
	FAIL("multi-path Obol draw sync should retain two database sources");

    const char *direct_path = "box.s";
    if (!qg_obol_sync_database_sources(gedp->dbip, &direct_path, 1,
	    QG_OBOL_DATABASE_SHADED, 123, &view))
	FAIL("direct Obol database sync should replace a source without BSG transaction input");
    source = source_for_path(controller, "box.s");
    if (!source ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::SHADED ||
	    source->sourceRevision.getValue() != 123 ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    source->getRealizedMeshCount() <= 0)
	FAIL("direct Obol database sync should preserve draw mode, revision, and mesh geometry");
    if (!qg_obol_remove_database_sources(&direct_path, 1, &view))
	FAIL("direct Obol database remove should remove one source without BSG transaction input");
    if (controller->getDatabaseSourceCount() != 1 ||
	    source_for_path(controller, "box.s"))
	FAIL("direct Obol database remove should leave only unrelated sources");

    controller->clearDatabaseSources();
    if (controller->getDatabaseSourceCount() != 0)
	FAIL("test setup should clear Obol sources before full redraw sync");
    struct ged_draw_transaction redraw_all =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    redraw_all.view = view.view();
    if (!apply_and_sync(gedp, &view, &redraw_all, 1))
	FAIL("GED redraw should rebuild Obol sources from retained draw state");
    if (controller->getDatabaseSourceCount() != 2)
	FAIL("Obol full redraw sync should rebuild retained GED draw paths");

    struct ged_draw_transaction clear_all =
	ged_draw_transaction_make(GED_DRAW_TXN_CLEAR, NULL);
    clear_all.view = view.view();
    if (!apply_and_sync(gedp, &view, &clear_all, 1))
	FAIL("GED clear should clear Obol database sources");
    if (controller->getDatabaseSourceCount() != 0)
	FAIL("Obol draw sync should clear all database sources");

    ged_close(gedp);
    bu_file_delete(dbpath);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
