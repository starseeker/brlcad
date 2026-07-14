/*          T E S T _ Q T C A D _ O B O L _ D R A W _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bv.h"

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
#include "ged/view.h"
#include "QgObolDrawSyncPrivate.h"
#include "qtcad/QgView.h"
#include "raytrace.h"
#include "rt/db_fullpath.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>

#include <QApplication>
#include <QCoreApplication>
#include <QImage>

#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <string>
#include <vector>

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
nonblack_pixel_count(const QImage &image)
{
    QImage rgba = image.convertToFormat(QImage::Format_RGBA8888);
    int count = 0;
    for (int y = 0; y < rgba.height(); y++) {
	const unsigned char *line = rgba.constScanLine(y);
	for (int x = 0; x < rgba.width(); x++) {
	    const unsigned char *p = line + x * 4;
	    if (p[0] > 2 || p[1] > 2 || p[2] > 2)
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
render_source(BRLObolViewController *controller, int index);

static void
collect_render_sources(SoNode *node,
	std::vector<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source =
	    static_cast<SoBRLDatabaseSource *>(node);
	if (std::find(sources.begin(), sources.end(), source) == sources.end())
	    sources.push_back(source);
	return;
    }
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	collect_render_sources(group->getChild(i), sources);
}

static std::vector<SoBRLDatabaseSource *>
render_sources(BRLObolViewController *controller)
{
    std::vector<SoBRLDatabaseSource *> sources;
    if (controller)
	collect_render_sources(controller->getRenderSceneRoot(), sources);
    if (sources.empty() && controller)
	collect_render_sources(controller->getSceneRoot(), sources);
    return sources;
}

static int
render_source_count(BRLObolViewController *controller)
{
    return static_cast<int>(render_sources(controller).size());
}

static SoBRLDatabaseSource *
render_source(BRLObolViewController *controller, int index)
{
    std::vector<SoBRLDatabaseSource *> sources = render_sources(controller);
    return index >= 0 && static_cast<size_t>(index) < sources.size() ?
	sources[static_cast<size_t>(index)] : NULL;
}

static SoBRLDatabaseSource *
source_for_path(BRLObolViewController *controller, const char *path)
{
    if (!controller || !path)
	return NULL;
    for (SoBRLDatabaseSource *source : render_sources(controller)) {
	if (source && test_path_equal(source->path.getValue().getString(),
		path))
	    return source;
    }
    return NULL;
}

static SoBRLDatabaseSource *
source_for_path_mode(BRLObolViewController *controller,
	const char *path,
	int draw_mode)
{
    if (!controller || !path)
	return NULL;
    for (SoBRLDatabaseSource *source : render_sources(controller)) {
	if (source && test_path_equal(source->path.getValue().getString(),
		path) && source->drawMode.getValue() == draw_mode)
	    return source;
    }
    return NULL;
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

    QgView view(NULL, QgViewType::SW);
    view.resize(180, 140);
	ged_view_active_ctx_set(gedp, view.viewContext());
	(void)ged_view_context_host_attach(gedp, view.viewContext());

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
    if (render_source_count(controller) != 1)
	FAIL("observer-synced GED draw should create one Obol database source");
    SoBRLDatabaseSource *observerSource = source_for_path(controller, "box.s");
    if (!observerSource ||
	    observerSource->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    !observerSource->hasRealizedWireGeometry())
	FAIL("observer-synced GED draw should realize Obol wire geometry");

    controller->getViewport()->viewAll();
    controller->requestRender("observer-draw-visible");
    QCoreApplication::processEvents();
    QImage observerImage;
    view.get_viewport_image(observerImage);
    int observerLit = observerImage.isNull() ? 0 : lit_pixel_count(observerImage);
    if (observerImage.isNull() || observerLit < 10) {
	fprintf(stderr,
		"FAIL: observer-synced GED draw should be visible through qtcad capture (null=%d lit=%d)\n",
		observerImage.isNull() ? 1 : 0, observerLit);
	return 1;
    }

    SoCamera *camera = controller->getCamera();
    if (!camera)
	FAIL("qtcad Obol controller should expose a camera for view sync");
    point_t offcenter = {100.0, 0.0, 0.0};
    bv_center_set(bv_context_view(static_cast<struct bv_context *>(view.viewContext())), offcenter);
    bv_scale_set(bv_context_view(static_cast<struct bv_context *>(view.viewContext())), 250.0);
    view.need_update(QG_VIEW_REFRESH);
    SbVec3f offTargetCamera = camera->position.getValue();
    if (offTargetCamera[0] < 50.0f)
	FAIL("qtcad refresh should sync Obol camera from GED view state");

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
    QCoreApplication::processEvents();
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
    if (render_source_count(controller) != 0)
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
    draw_box.view = view.viewContext();
    draw_box.appearance = &box_appearance;
    if (!apply_and_sync(gedp, &view, &draw_box, 1))
	FAIL("GED draw should sync a wire Obol database source");
    if (render_source_count(controller) != 1)
	FAIL("Obol draw sync should create one database source");
    SoBRLDatabaseSource *source = render_source(controller, 0);
    if (!source ||
	    !test_path_equal(source->path.getValue().getString(), "box.s") ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    !source->hasRealizedWireGeometry())
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
    if (source_for_path(controller, "box.s") != source)
	FAIL("Obol draw sync should publish box.s as the authoritative source owner");

    controller->getViewport()->viewAll();
    controller->requestRender("draw-sync-visible");
    QCoreApplication::processEvents();
    QImage visibleImage;
    view.get_viewport_image(visibleImage);
    if (visibleImage.isNull() || nonblack_pixel_count(visibleImage) < 10)
	FAIL("Obol-synced GED draw should be visible through qtcad capture");

    struct ged_draw_transaction stale_box =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "box.s");
    stale_box.view = view.viewContext();
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
    erase_box.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &erase_box, 1))
	FAIL("GED erase should remove an Obol database source");
    if (render_source_count(controller) != 0)
	FAIL("Obol draw sync should remove erased database sources");
    if (scene_display_summary_by_path(controller, "box.s",
	    BRLObolSceneTreeSummary::NODE_GROUP, NULL))
	FAIL("Obol draw sync should prune empty GED draw groups after erase");

    struct ged_draw_appearance_settings shaded_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    shaded_appearance.draw_mode = GED_DRAW_MODE_SHADED;
    struct ged_draw_transaction draw_ball =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "ball.s");
    draw_ball.view = view.viewContext();
    draw_ball.appearance = &shaded_appearance;
    int drew_ball = apply_and_sync(gedp, &view, &draw_ball, 1);
    if (!drew_ball)
	FAIL("GED shaded draw should sync a shaded Obol database source");
    source = render_source(controller, 0);
    if (!source ||
	    !test_path_equal(source->path.getValue().getString(), "ball.s") ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::SHADED ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    !source->hasRealizedMeshGeometry())
	FAIL("shaded Obol database source should preserve draw mode and mesh geometry");

    const char *paths[2] = {"box.s", "ball.s"};
    struct ged_draw_appearance_settings wire_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    wire_appearance.draw_mode = GED_DRAW_MODE_WIRE;
    wire_appearance.mixed_modes = 1;
    struct ged_draw_transaction draw_both =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, NULL);
    draw_both.view = view.viewContext();
    draw_both.paths = paths;
    draw_both.path_count = 2;
    draw_both.appearance = &wire_appearance;
    int drew_both = apply_and_sync(gedp, &view, &draw_both, 1);
    if (!drew_both)
	FAIL("multi-path GED draw should sync multiple Obol database sources");
    if (render_source_count(controller) != 3 ||
	    !source_for_path_mode(controller, "box.s",
		SoBRLDatabaseSource::WIREFRAME) ||
	    !source_for_path_mode(controller, "ball.s",
		SoBRLDatabaseSource::WIREFRAME) ||
	    !source_for_path_mode(controller, "ball.s",
		SoBRLDatabaseSource::SHADED))
	FAIL("multi-path Obol draw sync should retain shared and representation-specific database sources");

    SoBRLDatabaseSource *boxWireBeforeRedraw =
	source_for_path_mode(controller, "box.s", SoBRLDatabaseSource::WIREFRAME);
    SoBRLDatabaseSource *ballWireBeforeRedraw =
	source_for_path_mode(controller, "ball.s", SoBRLDatabaseSource::WIREFRAME);
    SoBRLDatabaseSource *ballShadedBeforeRedraw =
	source_for_path_mode(controller, "ball.s", SoBRLDatabaseSource::SHADED);
    struct ged_draw_transaction redraw_all =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    redraw_all.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &redraw_all, 1))
	FAIL("GED redraw should invalidate the retained Obol render");
    if (render_source_count(controller) != 3 ||
	    source_for_path_mode(controller, "box.s",
		SoBRLDatabaseSource::WIREFRAME) != boxWireBeforeRedraw ||
	    source_for_path_mode(controller, "ball.s",
		SoBRLDatabaseSource::WIREFRAME) != ballWireBeforeRedraw ||
	    source_for_path_mode(controller, "ball.s",
		SoBRLDatabaseSource::SHADED) != ballShadedBeforeRedraw)
	FAIL("Obol redraw should preserve retained database sources");

    struct ged_draw_transaction clear_all =
	ged_draw_transaction_make(GED_DRAW_TXN_CLEAR, NULL);
    clear_all.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &clear_all, 1))
	FAIL("GED clear should clear Obol database sources");
    if (render_source_count(controller) != 0)
	FAIL("Obol draw sync should clear all database sources");

    struct ged_draw_transaction draw_nested =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "pair.c/box.s");
    draw_nested.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &draw_nested, 1))
	FAIL("nested GED draw should sync an Obol database source");
    if (render_source_count(controller) != 1 ||
	    !source_for_path(controller, "pair.c/box.s"))
	FAIL("nested Obol draw sync should retain one full-path database source");
    source = source_for_path(controller, "pair.c/box.s");
    if (!source ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME)
	FAIL("full-path Obol draw sync should publish the nested source owner");

    struct ged_draw_transaction erase_nested =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "pair.c/box.s");
    erase_nested.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &erase_nested, 1))
	FAIL("nested GED erase should remove an Obol database source");
    if (render_source_count(controller) != 0 ||
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
