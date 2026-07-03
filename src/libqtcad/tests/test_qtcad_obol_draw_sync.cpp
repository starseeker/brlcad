/*          T E S T _ Q T C A D _ O B O L _ D R A W _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/lod_realization.h"
#include "brlobol/view_controller.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "ged.h"
#include "ged/draw.h"
#include "QgLegacyViewContext.h"
#include "QgObolDatabaseSyncPrivate.h"
#include "QgObolDrawSyncPrivate.h"
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgView.h"
#include "raytrace.h"
#include "rt/db_fullpath.h"
#include "rt/view.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoCamera.h>

#include <QApplication>
#include <QImage>

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <string>

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
    struct wmember pair;
    BU_LIST_INIT(&pair.l);
    ret = ret &&
	mk_addmember("box.s", &pair.l, NULL, WMOP_UNION) != NULL &&
	mk_lcomb(wdbp, "pair.c", &pair, 0, NULL, NULL, NULL, 0) == 0;
    wdb_close(wdbp);
    return ret;
}

static uint32_t
test_fold_revision(uint64_t revision)
{
    if (!revision)
	return 0;

    uint32_t folded = (uint32_t)(revision ^ (revision >> 32));
    return folded ? folded : 1;
}

static const char *
test_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static int
test_path_equal(const char *a, const char *b)
{
    if (!a || !b)
	return 0;
    if (BU_STR_EQUAL(a, b))
	return 1;
    return BU_STR_EQUAL(test_skip_leading_slash(a),
	    test_skip_leading_slash(b));
}

static int
test_shape_record_matches_path(const struct ged_draw_shape_record *record,
	const char *path)
{
    if (!record || !path || !path[0])
	return 0;
    if (test_path_equal(record->display_name, path) ||
	    test_path_equal(record->leaf_name, path))
	return 1;

    if (!record->fullpath || record->fullpath->fp_len <= 0)
	return 0;

    char *recordPath = db_path_to_string(record->fullpath);
    if (!recordPath)
	return 0;
    int matched = test_path_equal(recordPath, path);
    bu_free(recordPath, "qtcad Obol draw-sync test record path");
    return matched;
}

struct test_shape_record_path_context {
    const char *path;
    struct ged_draw_shape_record *out;
    int found;
};

static int
test_shape_record_by_path_cb(const struct ged_draw_shape_record *record,
	void *userdata)
{
    struct test_shape_record_path_context *ctx =
	(struct test_shape_record_path_context *)userdata;
    if (!ctx || !record || !test_shape_record_matches_path(record,
	    ctx->path))
	return 1;

    if (ctx->out)
	*ctx->out = *record;
    ctx->found = 1;
    return 0;
}

static int
test_shape_record_by_path(struct ged *gedp,
	const char *path,
	struct ged_draw_shape_record *out)
{
    if (!gedp || !path || !path[0])
	return 0;

    struct test_shape_record_path_context ctx;
    ctx.path = path;
    ctx.out = out;
    ctx.found = 0;
    ged_draw_foreach_shape_record(gedp, test_shape_record_by_path_cb, &ctx);
    return ctx.found;
}

static SoBRLDatabaseSource *
source_for_path(BRLObolViewController *controller, const char *path)
{
    if (!controller || !path)
	return NULL;
    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (source && test_path_equal(source->path.getValue().getString(),
		path))
	    return source;
    }
    return NULL;
}

static std::string
source_instance_key_for_view(const char *view_name, const char *path)
{
    std::string key("ged-view:");
    if (view_name)
	key += view_name;
    key += ":";
    key += test_skip_leading_slash(path);
    return key;
}

static std::string
source_mode_instance_key(const char *path, int mode)
{
    std::string key(test_skip_leading_slash(path));
    if (mode != GED_DRAW_MODE_WIRE) {
	char mode_buf[64] = {0};
	snprintf(mode_buf, sizeof(mode_buf), ":ged-draw-mode:%d", mode);
	key += mode_buf;
    }
    return key;
}

static SoBRLDatabaseSource *
source_for_instance(BRLObolViewController *controller, const char *instanceKey)
{
    if (!controller || !instanceKey)
	return NULL;
    return controller->findDatabaseSourceInstance(instanceKey);
}

static int
scene_display_summary_by_path(BRLObolViewController *controller,
	const char *path,
	int nodeKind,
	BRLObolSceneDisplaySummary *out)
{
    if (!controller || !controller->getSceneController() || !path)
	return 0;

    const SoBRLSceneController *scene = controller->getSceneController();
    for (int i = 0; i < scene->getSceneDisplaySummaryCount(); i++) {
	BRLObolSceneDisplaySummary summary;
	if (!scene->getSceneDisplaySummary(i, summary) || !summary.valid)
	    continue;
	if (summary.nodeKind == nodeKind &&
		test_path_equal(summary.path.getString(), path)) {
	    if (out)
		*out = summary;
	    return 1;
	}
    }
    return 0;
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
    ctx->changed += qg_obol_sync_ged_draw_transaction(gedp, txn, result,
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
    int changed = qg_obol_sync_ged_draw_transaction(gedp, txn, &result, view);
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

    QgView view(NULL, QgView_SW);
    view.resize(180, 140);
    qg_legacy_view_ged_active_set(gedp, view.view());

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
    void *view_ctx = qg_legacy_view_to_context(view.view());
    rt_view_context_center_set(view_ctx, offcenter);
    rt_view_context_scale_set(view_ctx, 250.0);
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

    struct ged_draw_appearance_settings box_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    box_appearance.color_override = 1;
    box_appearance.color[0] = 10;
    box_appearance.color[1] = 20;
    box_appearance.color[2] = 30;
    box_appearance.s_line_width = 5;
    struct ged_draw_transaction draw_box =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "box.s");
    draw_box.view = qg_legacy_view_to_context(view.view());
    draw_box.appearance = &box_appearance;
    if (!apply_and_sync(gedp, &view, &draw_box, 1))
	FAIL("GED draw should sync a wire Obol database source");
    if (controller->getDatabaseSourceCount() != 1)
	FAIL("Obol draw sync should create one database source");
    SoBRLDatabaseSource *source = controller->getDatabaseSource(0);
    if (!source ||
	    !test_path_equal(source->path.getValue().getString(), "box.s") ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    source->getRealizedShapeCount() <= 0)
	FAIL("wire Obol database source should preserve path and realized geometry");
    struct ged_draw_shape_record box_record;
    if (!test_shape_record_by_path(gedp, "box.s", &box_record))
	FAIL("GED draw should expose neutral source state for box.s");
    if (box_record.source_revision &&
	    source->sourceRevision.getValue() !=
	    test_fold_revision(box_record.source_revision))
	FAIL("Obol draw sync should copy nonzero GED source revision");
    if (source->inputsRevision.getValue() !=
	    test_fold_revision(box_record.inputs_revision) ||
	    (source->visible.getValue() ? 1 : 0) != box_record.visible ||
	    (source->highlighted.getValue() ? 1 : 0) != box_record.highlighted ||
	    source->lineWidth.getValue() != box_record.line_width ||
	    fabs(source->transparency.getValue() - box_record.transparency) >
	    1.0e-6)
	FAIL("Obol draw sync should seed source display state from GED records");
    void *box_scene_ctx = ged_draw_shape_ref_context(gedp, box_record.ref);
    struct ged_draw_scene_display_summary box_display;
    if (!box_scene_ctx ||
	    !ged_draw_scene_context_display_summary(box_scene_ctx,
		&box_display) ||
	    !box_display.valid ||
	    !box_display.material_valid)
	FAIL("GED draw should expose neutral display/material state for box.s");
    if (source->lineStyle.getValue() != box_display.line_style ||
	    (source->materialColorValid.getValue() ? 1 : 0) !=
	    box_display.material_valid)
	FAIL("Obol draw sync should copy GED display line/material validity");
    SbColor sourceMaterial = source->materialColor.getValue();
    if (fabsf(sourceMaterial[0] -
		(float)box_display.material_color[0] / 255.0f) > 1.0e-6f ||
	    fabsf(sourceMaterial[1] -
		(float)box_display.material_color[1] / 255.0f) > 1.0e-6f ||
	    fabsf(sourceMaterial[2] -
		(float)box_display.material_color[2] / 255.0f) > 1.0e-6f)
	FAIL("Obol draw sync should copy GED material color");
    struct ged_draw_shape_material_summary box_material;
    if (ged_draw_shape_ref_material_summary(gedp, box_record.ref,
	    &box_material) && box_material.valid &&
	    source->materialRevision.getValue() !=
	    test_fold_revision(box_material.material_revision))
	FAIL("Obol draw sync should copy GED material revision");
    BRLObolSceneDisplaySummary box_group_display;
    if (!scene_display_summary_by_path(controller, "box.s",
	    BRLObolSceneTreeSummary::NODE_GROUP, &box_group_display) ||
	    !box_group_display.hasDrawIntent ||
	    box_group_display.intentDrawMode != BRLOBOL_LOD_DRAW_WIRE ||
	    !box_group_display.visible ||
	    box_group_display.lineWidth != box_record.line_width ||
	    controller->getSceneController()->getGroupChildCount("box.s") != 1)
	FAIL("Obol draw sync should retain the GED draw group around the source");

    controller->getViewport()->viewAll();
    controller->requestRender("draw-sync-visible");
    QImage visibleImage;
    view.get_viewport_image(visibleImage);
    if (visibleImage.isNull() || lit_pixel_count(visibleImage) < 10)
	FAIL("Obol-synced GED draw should be visible through qtcad capture");

    struct ged_draw_transaction stale_box =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "box.s");
    stale_box.view = qg_legacy_view_to_context(view.view());
    if (!apply_and_sync(gedp, &view, &stale_box, 1))
	FAIL("GED stale-source event should refresh Obol source state");
    struct ged_draw_shape_record stale_record;
    if (!test_shape_record_by_path(gedp, "box.s", &stale_record) ||
	    !stale_record.source_revision)
	FAIL("GED stale-source event should expose a nonzero source revision");
    source = source_for_path(controller, "box.s");
    if (!source)
	FAIL("Obol stale-source sync should retain the box source");
    uint32_t expectedSourceRevision =
	test_fold_revision(stale_record.source_revision);
    uint32_t expectedInputsRevision =
	test_fold_revision(stale_record.inputs_revision);
    if (source->sourceRevision.getValue() != expectedSourceRevision ||
	    source->inputsRevision.getValue() != expectedInputsRevision ||
	    source->realizedSourceRevision.getValue() != expectedSourceRevision ||
	    source->realizedInputsRevision.getValue() != expectedInputsRevision ||
	    source->stale.getValue() ||
	    source->staleReason.getValue() != SoBRLDatabaseSource::STALE_NONE)
	FAIL("Obol stale-source sync should realize current GED source revisions");
    BRLObolRealizedShapeSummary shape_summary;
    if (source->getRealizedShapeSummaryCount() <= 0 ||
	    !source->getRealizedShapeSummary(0, shape_summary) ||
	    shape_summary.ownerSourceRevision != expectedSourceRevision ||
	    shape_summary.ownerInputsRevision != expectedInputsRevision ||
	    shape_summary.ownerSourceStale)
	FAIL("Obol realized summary should carry current GED source lineage");

    struct ged_draw_transaction erase_box =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "box.s");
    erase_box.view = qg_legacy_view_to_context(view.view());
    if (!apply_and_sync(gedp, &view, &erase_box, 1))
	FAIL("GED erase should remove an Obol database source");
    if (controller->getDatabaseSourceCount() != 0)
	FAIL("Obol draw sync should remove erased database sources");
    if (scene_display_summary_by_path(controller, "box.s",
	    BRLObolSceneTreeSummary::NODE_GROUP, NULL))
	FAIL("Obol draw sync should prune empty GED draw groups after erase");

    struct ged_draw_appearance_settings shaded_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    shaded_appearance.draw_mode = GED_DRAW_MODE_SHADED;
    struct ged_draw_transaction draw_ball =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "ball.s");
    draw_ball.view = qg_legacy_view_to_context(view.view());
    draw_ball.appearance = &shaded_appearance;
    int drew_ball = apply_and_sync(gedp, &view, &draw_ball, 1);
    if (!drew_ball)
	FAIL("GED shaded draw should sync a shaded Obol database source");
    source = controller->getDatabaseSource(0);
    if (!source ||
	    !test_path_equal(source->path.getValue().getString(), "ball.s") ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::SHADED ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    source->getRealizedMeshCount() <= 0)
	FAIL("shaded Obol database source should preserve draw mode and mesh geometry");

    const char *paths[2] = {"box.s", "ball.s"};
    struct ged_draw_appearance_settings wire_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    wire_appearance.draw_mode = GED_DRAW_MODE_WIRE;
    struct ged_draw_transaction draw_both =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, NULL);
    draw_both.view = qg_legacy_view_to_context(view.view());
    draw_both.paths = paths;
    draw_both.path_count = 2;
    draw_both.appearance = &wire_appearance;
    int drew_both = apply_and_sync(gedp, &view, &draw_both, 1);
    if (!drew_both)
	FAIL("multi-path GED draw should sync multiple Obol database sources");
    const std::string shaded_ball_instance =
	source_mode_instance_key("ball.s", GED_DRAW_MODE_SHADED);
    if (controller->getDatabaseSourceCount() != 3 ||
	    !source_for_instance(controller, "box.s") ||
	    !source_for_instance(controller, "ball.s") ||
	    !source_for_instance(controller, shaded_ball_instance.c_str()))
	FAIL("multi-path Obol draw sync should retain shared and representation-specific database sources");

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
    if (controller->getDatabaseSourceCount() != 2 ||
	    source_for_path(controller, "box.s"))
	FAIL("direct Obol database remove should leave only unrelated ball sources");

    if (!source_for_instance(controller, "ball.s") ||
	    !source_for_instance(controller, shaded_ball_instance.c_str()))
	FAIL("direct Obol database remove should retain shared and representation-specific ball source owners");
    if (!rt_view_context_name_set(view_ctx, "QV0") ||
	    rt_view_context_independent_scope_is_null(view_ctx, 1) ||
	    !rt_view_context_is_independent(view_ctx))
	FAIL("qtcad direct source-owner parity test should create an independent view scope");
    const std::string scoped_box =
	source_instance_key_for_view("QV0", "box.s");
    const std::string scoped_ball =
	source_instance_key_for_view("QV0", "ball.s");
    if (!qg_obol_sync_database_sources(gedp->dbip, &direct_path, 1,
	    QG_OBOL_DATABASE_WIREFRAME, 321, &view))
	FAIL("direct Obol database sync should create scoped independent-view owners");
    SoBRLDatabaseSource *scopedBoxSource =
	source_for_instance(controller, scoped_box.c_str());
    if (!scopedBoxSource ||
	    !BU_STR_EQUAL(scopedBoxSource->path.getValue().getString(),
		"box.s") ||
	    !source_for_instance(controller, "ball.s") ||
	    !source_for_instance(controller, shaded_ball_instance.c_str()) ||
	    controller->getDatabaseSourceCount() != 3)
	FAIL("scoped direct Obol sync should coexist with shared source owners and mode-specific representations");
    if (!qg_obol_remove_database_sources(&direct_path, 1, &view))
	FAIL("direct Obol database remove should target scoped independent-view owners");
    if (source_for_instance(controller, scoped_box.c_str()) ||
	    !source_for_instance(controller, "ball.s") ||
	    !source_for_instance(controller, shaded_ball_instance.c_str()) ||
	    controller->getDatabaseSourceCount() != 2)
	FAIL("scoped direct Obol remove should leave shared source owners and mode-specific representations intact");
    if (!qg_obol_sync_database_sources(gedp->dbip, paths, 2,
	    QG_OBOL_DATABASE_WIREFRAME, 654, &view))
	FAIL("direct Obol database sync should create multiple scoped owners");
    if (!source_for_instance(controller, scoped_box.c_str()) ||
	    !source_for_instance(controller, scoped_ball.c_str()) ||
	    !source_for_instance(controller, "ball.s") ||
	    !source_for_instance(controller, shaded_ball_instance.c_str()) ||
	    controller->getDatabaseSourceCount() != 4)
	FAIL("scoped direct Obol sync should retain shared, scoped, and mode-specific owners separately");
    if (!qg_obol_clear_database_sources(&view))
	FAIL("direct Obol database clear should target the active source owner scope");
    if (source_for_instance(controller, scoped_box.c_str()) ||
	    source_for_instance(controller, scoped_ball.c_str()) ||
	    !source_for_instance(controller, "ball.s") ||
	    !source_for_instance(controller, shaded_ball_instance.c_str()) ||
	    controller->getDatabaseSourceCount() != 2)
	FAIL("scoped direct Obol clear should leave shared source owners and mode-specific representations intact");
    rt_view_context_independent_scope_destroy(view_ctx);
    if (rt_view_context_is_independent(view_ctx))
	FAIL("qtcad direct source-owner parity test should restore shared view scope");

    controller->clearDatabaseSources();
    if (controller->getDatabaseSourceCount() != 0)
	FAIL("test setup should clear Obol sources before full redraw sync");
    struct ged_draw_transaction redraw_all =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    redraw_all.view = qg_legacy_view_to_context(view.view());
    if (!apply_and_sync(gedp, &view, &redraw_all, 1))
	FAIL("GED redraw should rebuild Obol sources from retained draw state");
    if (controller->getDatabaseSourceCount() != 2)
	FAIL("Obol full redraw sync should rebuild retained GED draw paths");

    struct ged_draw_transaction clear_all =
	ged_draw_transaction_make(GED_DRAW_TXN_CLEAR, NULL);
    clear_all.view = qg_legacy_view_to_context(view.view());
    if (!apply_and_sync(gedp, &view, &clear_all, 1))
	FAIL("GED clear should clear Obol database sources");
    if (controller->getDatabaseSourceCount() != 0)
	FAIL("Obol draw sync should clear all database sources");

    struct ged_draw_transaction draw_nested =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "pair.c/box.s");
    draw_nested.view = qg_legacy_view_to_context(view.view());
    if (!apply_and_sync(gedp, &view, &draw_nested, 1))
	FAIL("nested GED draw should sync an Obol database source");
    if (controller->getDatabaseSourceCount() != 1 ||
	    !source_for_path(controller, "pair.c/box.s"))
	FAIL("nested Obol draw sync should retain one full-path database source");
    BRLObolSceneDisplaySummary nested_group_display;
    if (!scene_display_summary_by_path(controller, "pair.c/box.s",
	    BRLObolSceneTreeSummary::NODE_GROUP, &nested_group_display) ||
	    !nested_group_display.hasDrawIntent ||
	    nested_group_display.intentDrawMode != BRLOBOL_LOD_DRAW_WIRE ||
	    controller->getSceneController()->getGroupChildCount(
		"pair.c/box.s") != 1)
	FAIL("full-path Obol draw sync should retain the GED draw group around the source");

    struct ged_draw_transaction erase_nested =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "pair.c/box.s");
    erase_nested.view = qg_legacy_view_to_context(view.view());
    if (!apply_and_sync(gedp, &view, &erase_nested, 1))
	FAIL("nested GED erase should remove an Obol database source");
    if (controller->getDatabaseSourceCount() != 0 ||
	    scene_display_summary_by_path(controller, "pair.c/box.s",
		BRLObolSceneTreeSummary::NODE_GROUP, NULL) ||
	    scene_display_summary_by_path(controller, "pair.c",
		BRLObolSceneTreeSummary::NODE_GROUP, NULL))
	FAIL("nested Obol erase should leave no synthetic GED group ancestors");

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
