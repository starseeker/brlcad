/*          T E S T _ Q T C A D _ O B O L _ D R A W _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bv.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BLodRealization.h"
#include "BObol/BViewController.h"
#include "BObol/BSceneController.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "ged.h"
#include "ged/scene.h"
#include "ged/display_obol_private.h"
#include "ged/draw.h"
#include "ged/event.h"
#include "ged/view.h"
#include "QgSceneSyncPrivate.h"
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
	mk_addmember("ball.s", &pair.l, NULL, WMOP_UNION) != NULL &&
	mk_lcomb(wdbp, "pair.c", &pair, 0, NULL, NULL, NULL, 0) == 0;
    wdb_close(wdbp);
    return ret;
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
test_shape_record_matches_path(const struct ged_scene_occurrence_info *record,
	const char *path)
{
    if (!record || !path || !path[0])
	return 0;
    if (test_path_equal(record->path, path) ||
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
    struct ged_scene_occurrence_info *out;
    int found;
};

static int
test_shape_record_by_path_cb(const struct ged_scene_occurrence_info *record,
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
	struct ged_scene_occurrence_info *out)
{
    if (!gedp || !path || !path[0])
	return 0;

    struct test_shape_record_path_context ctx;
    ctx.path = path;
    ctx.out = out;
    ctx.found = 0;
    ged_scene_occurrences_visit(gedp, test_shape_record_by_path_cb, &ctx);
    return ctx.found;
}

static SoBRLDatabaseSource *
render_source(BObolViewController *controller, int index);

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
render_sources(BObolViewController *controller)
{
    std::vector<SoBRLDatabaseSource *> sources;
    if (controller)
	collect_render_sources(controller->getRenderSceneRoot(), sources);
    if (sources.empty() && controller)
	collect_render_sources(controller->getSceneRoot(), sources);
    return sources;
}

static int
render_source_count(BObolViewController *controller)
{
    return static_cast<int>(render_sources(controller).size());
}

static int
visible_compact_occurrence_count(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;
    int count = 0;
    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactOccurrence occurrence;
	if (source->getCompactOccurrence(i, occurrence) &&
		occurrence.summary.visible)
	    count++;
    }
    return count;
}

static int
compact_occurrence_for_path(SoBRLDatabaseSource *source, const char *path,
	BObolCompactOccurrence *out)
{
    if (!source || !path)
	return 0;
    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactOccurrence occurrence;
	if (!source->getCompactOccurrence(i, occurrence) ||
	    !test_path_equal(occurrence.summary.path.getString(), path))
	    continue;
	if (out)
	    *out = occurrence;
	return 1;
    }
    return 0;
}

static SoBRLDatabaseSource *
render_source(BObolViewController *controller, int index)
{
    std::vector<SoBRLDatabaseSource *> sources = render_sources(controller);
    return index >= 0 && static_cast<size_t>(index) < sources.size() ?
	sources[static_cast<size_t>(index)] : NULL;
}

static SoBRLDatabaseSource *
source_for_path(BObolViewController *controller, const char *path)
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
source_for_path_mode(BObolViewController *controller,
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
scene_display_summary_by_path(BObolViewController *controller,
	const char *path,
	int nodeKind,
	BObolSceneDisplaySummary *out)
{
    if (!controller || !controller->getSceneController() || !path)
	return 0;

    const BObolSceneController *scene = controller->getSceneController();
    for (int i = 0; i < scene->getSceneDisplaySummaryCount(); i++) {
	BObolSceneDisplaySummary summary;
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
	const struct ged_scene_delta *delta,
	void *client_data)
{
    struct draw_observer_sync_context *ctx =
	static_cast<struct draw_observer_sync_context *>(client_data);
    if (!ctx || !ctx->view)
	return;
    ctx->calls++;
    ctx->changed += qg_scene_delta_notify(gedp, delta, ctx->view);
}

static int
scene_result_matches(enum ged_scene_status status,
	struct ged_scene_result *result,
	int expect_changed)
{
    const int changed = ged_scene_result_changed(result);
    const int matches = status == GED_SCENE_OK &&
	(expect_changed ? changed != 0 : changed == 0);
    ged_scene_result_destroy(result);
    return matches;
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

    /* A renderer attachment must consume the complete current semantic
     * snapshot, including visibility-frontier changes committed while the
     * scene was headless. */
    struct ged_view_context *headless_view = ged_view_context_create();
    if (!headless_view || !ged_view_context_host_attach(gedp, headless_view))
	FAIL("headless snapshot test should create an attached view context");
    ged_view_active_ctx_set(gedp, headless_view);
    const char *headless_paths[] = {"pair.c"};
    struct ged_scene_draw_request headless_draw;
    ged_scene_draw_request_init(&headless_draw);
    headless_draw.view = headless_view;
    headless_draw.paths = headless_paths;
    headless_draw.path_count = 1;
    headless_draw.realization.mode = GED_SCENE_REALIZE_EAGER;
    struct ged_scene_result *headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &headless_draw,
	    headless_result), headless_result, 1))
	FAIL("headless draw intent should commit without a renderer");
    struct ged_scene_path_request headless_path;
    ged_scene_path_request_init(&headless_path);
    headless_path.view = headless_view;
    headless_path.path = "box.s";
    headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_visibility_set(gedp,
	    &headless_path, 0, headless_result), headless_result, 0))
	FAIL("exact compact path matching must not fall back to a leaf basename");
    headless_path.path = "pair.c/ball.s";
    headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_opacity_set(gedp, &headless_path,
	    0.25, headless_result), headless_result, 1))
	FAIL("headless compact opacity should be retained before realization");
    headless_path.path = "ball.s";
    headless_path.match = GED_SCENE_PATH_MATCH_OBJECT;
    headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_highlight_set(gedp, &headless_path,
	    1, headless_result), headless_result, 1))
	FAIL("headless object highlight should be retained before realization");
    struct ged_scene_erase_request headless_erase;
    ged_scene_erase_request_init(&headless_erase);
    headless_erase.view = headless_view;
    headless_erase.path = "pair.c/box.s";
    headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp, &headless_erase,
	    headless_result), headless_result, 1))
	FAIL("headless subpath erase should commit a visibility frontier");
    BObolViewController *controller = view.obolViewController();
    if (!controller || render_source_count(controller) != 0)
	FAIL("headless semantic commits must not construct renderer sources");
    if (!ged_view_context_obol_endpoint_set(headless_view,
	    view.displayEndpoint(), 0))
	FAIL("delayed endpoint attachment should accept a semantic snapshot");
    SoBRLDatabaseSource *headless_source =
	source_for_path(controller, "pair.c");
    if (!headless_source || headless_source->getCompactInstanceCount() != 2 ||
	    visible_compact_occurrence_count(headless_source) != 1)
	FAIL("delayed endpoint snapshot should preserve the headless erase frontier");
    BObolCompactOccurrence headless_ball;
    if (!compact_occurrence_for_path(headless_source, "pair.c/ball.s",
	    &headless_ball) || !headless_ball.summary.highlighted ||
	    fabsf(headless_ball.summary.transparency - 0.75f) > 1.0e-6f) {
	fprintf(stderr, "headless ball found=%d highlighted=%d transparency=%g\n",
	    compact_occurrence_for_path(headless_source, "pair.c/ball.s", NULL),
	    headless_ball.summary.highlighted ? 1 : 0,
	    headless_ball.summary.transparency);
	FAIL("delayed endpoint snapshot should preserve compact highlight and opacity");
    }

    /* Addressing one compact occurrence must remain available to editing and
     * picking clients without realizing sibling scene records. */
    struct ged_scene_path_request occurrence_request;
    ged_scene_path_request_init(&occurrence_request);
    occurrence_request.view = headless_view;
    occurrence_request.path = "pair.c/ball.s";
    occurrence_request.match = GED_SCENE_PATH_MATCH_EXACT;
    ged_scene_occurrence_ref occurrence = ged_scene_occurrence_resolve(gedp,
	&occurrence_request);
    struct ged_scene_occurrence_info occurrence_info;
    if (ged_scene_occurrence_ref_is_null(occurrence) ||
	!ged_scene_occurrence_get(gedp, occurrence, &occurrence_info) ||
	!occurrence_info.path ||
	!BU_STR_EQUAL(test_skip_leading_slash(occurrence_info.path),
	    "pair.c/ball.s") || render_source_count(controller) != 1)
	FAIL("exact compact occurrence resolution should retain one source and one semantic path");

    occurrence_request.path = "pair.c";
    occurrence_request.match = GED_SCENE_PATH_MATCH_SUBTREE;
    occurrence = ged_scene_occurrence_resolve(gedp, &occurrence_request);
    if (ged_scene_occurrence_ref_is_null(occurrence) ||
	!ged_scene_occurrence_get(gedp, occurrence, &occurrence_info) ||
	!occurrence_info.path ||
	!BU_STR_EQUAL(test_skip_leading_slash(occurrence_info.path),
	    "pair.c/ball.s"))
	FAIL("compact subtree occurrence resolution should skip hidden children in deterministic order");

    occurrence_request.path = "pair.c/missing.s";
    occurrence_request.match = GED_SCENE_PATH_MATCH_EXACT;
    if (!ged_scene_occurrence_ref_is_null(
	    ged_scene_occurrence_resolve(gedp, &occurrence_request)))
	FAIL("missing compact paths must not manufacture occurrence records");

    ged_scene_path_request_init(&headless_path);
    headless_path.view = headless_view;
    headless_path.path = "pair.c/ball.s";
    headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_visibility_set(gedp,
	    &headless_path, 0, headless_result), headless_result, 1))
	FAIL("compact visibility should update retained scene intent");
    headless_source = source_for_path(controller, "pair.c");
    if (!headless_source || visible_compact_occurrence_count(headless_source))
	FAIL("attached endpoint should apply compact visibility without leaf records");
    if (!compact_occurrence_for_path(headless_source, "pair.c/ball.s",
	    &headless_ball) || !headless_ball.summary.highlighted ||
	    fabsf(headless_ball.summary.transparency - 0.75f) > 1.0e-6f)
	FAIL("compact visibility must preserve prior highlight and opacity state");
    (void)ged_view_context_obol_endpoint_set(headless_view, NULL, 0);
    struct ged_scene_clear_request headless_clear;
    ged_scene_clear_request_init(&headless_clear);
    headless_clear.view = headless_view;
    headless_clear.scope = GED_SCENE_CLEAR_VIEW;
    headless_result = ged_scene_result_create();
    if (ged_scene_clear(gedp, &headless_clear, headless_result) != GED_SCENE_OK)
	FAIL("headless snapshot test should clear its view-scoped draw intent");
    ged_scene_result_destroy(headless_result);
    ged_view_context_free(headless_view);
    controller->clearDatabaseSources();

    struct ged_view_context *ged_view_ctx =
	ged_view_context_from_bv(view.viewContext());
	ged_view_active_ctx_set(gedp, ged_view_ctx);
	(void)ged_view_context_host_attach(gedp, ged_view_ctx);
    controller->clearDatabaseSources();

    struct draw_observer_sync_context obs = {&view, 0, 0};
    ged_scene_observer_token observerToken =
	ged_scene_observer_add(gedp, test_draw_observer, &obs);
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

    /* A nested erase/redraw is a visibility-frontier edit of the retained
     * top-level source.  It deliberately does not rebuild that source, but it
     * must invalidate the completed framebuffer or Qt will replay obsolete
     * pixels forever.  Exercise the actual GED observer bridge used by qged. */
    const char *observer_pair_paths[] = {"pair.c"};
    struct ged_scene_draw_request observer_pair_draw;
    ged_scene_draw_request_init(&observer_pair_draw);
    observer_pair_draw.view = ged_view_ctx;
    observer_pair_draw.paths = observer_pair_paths;
    observer_pair_draw.path_count = 1;
    observer_pair_draw.realization.mode = GED_SCENE_REALIZE_EAGER;
    struct ged_scene_result *observer_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &observer_pair_draw,
	    observer_result), observer_result, 1))
	FAIL("observer retained-frontier setup draw should succeed");
    SoBRLDatabaseSource *pairSource = source_for_path(controller, "pair.c");
    if (!pairSource || visible_compact_occurrence_count(pairSource) != 2)
	FAIL("observer retained-frontier setup should expose both occurrences");

    controller->clearRenderRequest();
    uint64_t frontierSerial = controller->renderRequestSerialGet();
    obs.calls = 0;
    obs.changed = 0;
    struct ged_scene_erase_request observer_nested_erase;
    ged_scene_erase_request_init(&observer_nested_erase);
    observer_nested_erase.view = ged_view_ctx;
    observer_nested_erase.path = "pair.c/box.s";
    observer_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp,
	    &observer_nested_erase, observer_result), observer_result, 1))
	FAIL("observer nested erase should succeed");
    SbBool frontierCapacityRelevant = TRUE;
    if (obs.calls <= 0 || obs.changed <= 0 ||
	visible_compact_occurrence_count(pairSource) != 1 ||
	controller->renderRequestSerialGet() <= frontierSerial ||
	!controller->consumeRenderRequest(NULL,
	    &frontierCapacityRelevant) || frontierCapacityRelevant)
	FAIL("observer nested erase should request a presentation-only frame");

    frontierSerial = controller->renderRequestSerialGet();
    obs.calls = 0;
    obs.changed = 0;
    const char *observer_nested_paths[] = {"pair.c/box.s"};
    struct ged_scene_draw_request observer_nested_draw;
    ged_scene_draw_request_init(&observer_nested_draw);
    observer_nested_draw.view = ged_view_ctx;
    observer_nested_draw.paths = observer_nested_paths;
    observer_nested_draw.path_count = 1;
    observer_nested_draw.realization.mode = GED_SCENE_REALIZE_EAGER;
    observer_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &observer_nested_draw,
	    observer_result), observer_result, 1))
	FAIL("observer nested redraw should succeed");
    frontierCapacityRelevant = TRUE;
    if (obs.calls <= 0 || obs.changed <= 0 ||
	visible_compact_occurrence_count(pairSource) != 2 ||
	controller->renderRequestSerialGet() <= frontierSerial ||
	!controller->consumeRenderRequest(NULL,
	    &frontierCapacityRelevant) || frontierCapacityRelevant)
	FAIL("observer nested redraw should request a presentation-only frame");

    struct ged_scene_erase_request observer_pair_erase;
    ged_scene_erase_request_init(&observer_pair_erase);
    observer_pair_erase.view = ged_view_ctx;
    observer_pair_erase.path = "pair.c";
    observer_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp,
	    &observer_pair_erase, observer_result), observer_result, 1) ||
	render_source_count(controller) != 0)
	FAIL("observer retained-frontier test should cleanly erase its root");
    if (ged_scene_observer_remove(gedp, observerToken) != 1)
	FAIL("GED draw observer should unregister after qtcad Obol sync test");
    controller->clearDatabaseSources();

    const char *box_paths[] = {"box.s"};
    struct ged_scene_draw_request draw_box;
    ged_scene_draw_request_init(&draw_box);
    draw_box.view = ged_view_ctx;
    draw_box.paths = box_paths;
    draw_box.path_count = 1;
    draw_box.realization.mode = GED_SCENE_REALIZE_EAGER;
    draw_box.style.color_override = 1;
    draw_box.style.color[0] = 10;
    draw_box.style.color[1] = 20;
    draw_box.style.color[2] = 30;
    draw_box.style.line_width = 5;
    struct ged_scene_result *scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &draw_box, scene_result),
	    scene_result, 1))
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
    struct ged_scene_occurrence_info box_record;
    if (!test_shape_record_by_path(gedp, "box.s", &box_record))
	FAIL("GED draw should expose neutral source state for box.s");
    if ((source->visible.getValue() ? 1 : 0) != box_record.visible ||
	    (source->highlighted.getValue() ? 1 : 0) != box_record.highlighted ||
	    source->lineWidth.getValue() != box_record.line_width ||
	    fabs(source->transparency.getValue() - (1.0 - box_record.opacity)) >
	    1.0e-6)
	FAIL("Obol draw sync should seed source display state from semantic occurrence state");
    if (!source->materialColorValid.getValue())
	FAIL("Obol draw sync should publish a valid database material color");
    if (source_for_path(controller, "box.s") != source)
	FAIL("Obol draw sync should publish box.s as the authoritative source owner");

    controller->getViewport()->viewAll();
    controller->requestRender("draw-sync-visible");
    QCoreApplication::processEvents();
    QImage visibleImage;
    view.get_viewport_image(visibleImage);
    if (visibleImage.isNull() || nonblack_pixel_count(visibleImage) < 10)
	FAIL("Obol-synced GED draw should be visible through qtcad capture");

    if (ged_event_notify_object_modified(gedp, "box.s", 1, NULL) != GED_EVENT_OK)
	FAIL("GED stale-source event should refresh Obol source state");
    struct ged_scene_occurrence_info stale_record;
    if (!test_shape_record_by_path(gedp, "box.s", &stale_record))
	FAIL("GED stale-source event should retain the semantic occurrence");
    source = source_for_path(controller, "box.s");
    if (!source)
	FAIL("Obol stale-source sync should retain the box source");
    uint32_t expectedSourceRevision = source->sourceRevision.getValue();
    uint32_t expectedInputsRevision = source->inputsRevision.getValue();
    if (!expectedSourceRevision)
	FAIL("Obol stale-source event should advance a nonzero source revision");
    if (source->sourceRevision.getValue() != expectedSourceRevision ||
	    source->inputsRevision.getValue() != expectedInputsRevision ||
	    source->realizedSourceRevision.getValue() != expectedSourceRevision ||
	    source->realizedInputsRevision.getValue() != expectedInputsRevision ||
	    source->stale.getValue() ||
	    source->staleReason.getValue() != SoBRLDatabaseSource::STALE_NONE)
	FAIL("Obol stale-source sync should realize current GED source revisions");
    BObolRealizedShapeSummary shape_summary;
    if (source->getRealizedShapeSummaryCount() <= 0 ||
	    !source->getRealizedShapeSummary(0, shape_summary) ||
	    shape_summary.ownerSourceRevision != expectedSourceRevision ||
	    shape_summary.ownerInputsRevision != expectedInputsRevision ||
	    shape_summary.ownerSourceStale)
	FAIL("Obol realized summary should carry current GED source lineage");

    struct ged_scene_erase_request erase_box;
    ged_scene_erase_request_init(&erase_box);
    erase_box.view = ged_view_ctx;
    erase_box.path = "box.s";
    scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp, &erase_box, scene_result),
	    scene_result, 1))
	FAIL("GED erase should remove an Obol database source");
    if (render_source_count(controller) != 0)
	FAIL("Obol draw sync should remove erased database sources");
    if (scene_display_summary_by_path(controller, "box.s",
	    BObolSceneTreeSummary::NODE_GROUP, NULL))
	FAIL("Obol draw sync should prune empty GED draw groups after erase");

    const char *ball_paths[] = {"ball.s"};
    struct ged_scene_draw_request draw_ball;
    ged_scene_draw_request_init(&draw_ball);
    draw_ball.view = ged_view_ctx;
    draw_ball.paths = ball_paths;
    draw_ball.path_count = 1;
    draw_ball.realization.mode = GED_SCENE_REALIZE_EAGER;
    draw_ball.style.draw_mode = GED_SCENE_DRAW_SHADED;
    scene_result = ged_scene_result_create();
    int drew_ball = scene_result_matches(
	ged_scene_draw(gedp, &draw_ball, scene_result), scene_result, 1);
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
    struct ged_scene_draw_request draw_both;
    ged_scene_draw_request_init(&draw_both);
    draw_both.view = ged_view_ctx;
    draw_both.paths = paths;
    draw_both.path_count = 2;
    draw_both.realization.mode = GED_SCENE_REALIZE_EAGER;
    draw_both.style.draw_mode = GED_SCENE_DRAW_WIRE;
    draw_both.style.mixed_modes = 1;
    scene_result = ged_scene_result_create();
    int drew_both = scene_result_matches(
	ged_scene_draw(gedp, &draw_both, scene_result), scene_result, 1);
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
    struct ged_scene_redraw_request redraw_all;
    ged_scene_redraw_request_init(&redraw_all);
    redraw_all.view = ged_view_ctx;
    scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_redraw(gedp, &redraw_all,
	    scene_result), scene_result, 1))
	FAIL("GED redraw should invalidate the retained Obol render");
    if (render_source_count(controller) != 3 ||
	    source_for_path_mode(controller, "box.s",
		SoBRLDatabaseSource::WIREFRAME) != boxWireBeforeRedraw ||
	    source_for_path_mode(controller, "ball.s",
		SoBRLDatabaseSource::WIREFRAME) != ballWireBeforeRedraw ||
	    source_for_path_mode(controller, "ball.s",
		SoBRLDatabaseSource::SHADED) != ballShadedBeforeRedraw)
	FAIL("Obol redraw should preserve retained database sources");

    struct ged_scene_clear_request clear_all;
    ged_scene_clear_request_init(&clear_all);
    scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_clear(gedp, &clear_all,
	    scene_result), scene_result, 1))
	FAIL("GED clear should clear Obol database sources");
    if (render_source_count(controller) != 0)
	FAIL("Obol draw sync should clear all database sources");

    const char *nested_paths[] = {"pair.c/box.s"};
    struct ged_scene_draw_request draw_nested;
    ged_scene_draw_request_init(&draw_nested);
    draw_nested.view = ged_view_ctx;
    draw_nested.paths = nested_paths;
    draw_nested.path_count = 1;
    draw_nested.realization.mode = GED_SCENE_REALIZE_EAGER;
    scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &draw_nested,
	    scene_result), scene_result, 1))
	FAIL("nested GED draw should sync an Obol database source");
    if (render_source_count(controller) != 1 ||
	    !source_for_path(controller, "pair.c/box.s"))
	FAIL("nested Obol draw sync should retain one full-path database source");
    source = source_for_path(controller, "pair.c/box.s");
    if (!source ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME)
	FAIL("full-path Obol draw sync should publish the nested source owner");

    struct ged_scene_erase_request erase_nested;
    ged_scene_erase_request_init(&erase_nested);
    erase_nested.view = ged_view_ctx;
    erase_nested.path = "pair.c/box.s";
    scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp, &erase_nested,
	    scene_result), scene_result, 1))
	FAIL("nested GED erase should remove an Obol database source");
    if (render_source_count(controller) != 0 ||
	    scene_display_summary_by_path(controller, "pair.c/box.s",
		BObolSceneTreeSummary::NODE_GROUP, NULL) ||
	    scene_display_summary_by_path(controller, "pair.c",
		BObolSceneTreeSummary::NODE_GROUP, NULL))
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
