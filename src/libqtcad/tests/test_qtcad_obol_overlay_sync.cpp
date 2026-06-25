/*      T E S T _ Q T C A D _ O B O L _ O V E R L A Y _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bg/line_layer.h"
#include "brlobol/hud_label_overlay.h"
#include "brlobol/line_layer_overlay.h"
#include "brlobol/view_controller.h"
#include "brlobol/vlist_shape.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/str.h"
#include "ged.h"
#include "QgObolOverlaySyncPrivate.h"
#include "qtcad/QgView.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/annex/HUD/nodes/SoHUDLabel.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

#include <QApplication>
#include <QImage>

#include <stdio.h>
#include <string.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

struct overlay_sync_context {
    QgView *view;
};

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
test_overlay_handler(struct ged *gedp,
	const char *name,
	const struct bg_line_layer_builder *builder,
	void *ctx)
{
    struct overlay_sync_context *syncCtx =
	static_cast<struct overlay_sync_context *>(ctx);
    return syncCtx ? qg_obol_sync_line_layer_overlay(gedp, name, builder,
	    syncCtx->view) : 0;
}

static int
test_hud_label_handler(struct ged *gedp,
	const struct ged_diagnostic_hud_label *label,
	void *ctx)
{
    struct overlay_sync_context *syncCtx =
	static_cast<struct overlay_sync_context *>(ctx);
    return syncCtx ? qg_obol_sync_hud_label_overlay(gedp, label,
	    syncCtx->view) : 0;
}

static int
make_overlay_sync_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;
    int ret = mk_id(wdbp, "qtcad Obol overlay sync test") == 0;
    wdb_close(wdbp);
    return ret;
}

static int
add_line_point(struct bg_line_layer_builder *builder,
	int r,
	int g,
	int b,
	double x,
	double y,
	double z,
	int command)
{
    point_t p;
    VSET(p, x, y, z);
    return bg_line_layer_builder_add(builder, r, g, b, p, command);
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_overlay_sync_tmp.g";
    if (!make_overlay_sync_db(dbpath))
	FAIL("failed to create qtcad Obol overlay-sync test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol overlay-sync test database");

    struct bg_line_layer_builder *builder = bg_line_layer_builder_create();
    if (!builder)
	FAIL("failed to allocate qtcad diagnostic line-layer builder");
    if (!add_line_point(builder, 255, 0, 0, -40.0, 0.0, 0.0, BG_GEOMETRY_LINE_MOVE) ||
	    !add_line_point(builder, 255, 0, 0, 40.0, 0.0, 0.0, BG_GEOMETRY_LINE_DRAW) ||
	    !add_line_point(builder, 0, 255, 255, 0.0, -40.0, 0.0, BG_GEOMETRY_LINE_MOVE) ||
	    !add_line_point(builder, 0, 255, 255, 0.0, 40.0, 0.0, BG_GEOMETRY_LINE_DRAW))
	FAIL("failed to populate qtcad diagnostic line-layer builder");

    QgView view(NULL, QgView_SW);
    view.resize(180, 140);
    qg_legacy_view_ged_active_set(gedp, view.view());

    QgView otherView(NULL, QgView_SW);
    otherView.resize(180, 140);
    if (qg_obol_sync_line_layer_overlay(gedp, "qtcad::wrong-view",
		builder, &otherView) != 0)
	FAIL("qtcad Obol overlay sync should ignore mismatched GED views");
    struct ged_diagnostic_hud_label wrongLabel = GED_DIAGNOSTIC_HUD_LABEL_INIT;
    wrongLabel.label_id = "qtcad::wrong-label";
    wrongLabel.text = "wrong view";
    if (qg_obol_sync_hud_label_overlay(gedp, &wrongLabel, &otherView) != 0)
	FAIL("qtcad Obol HUD label sync should ignore mismatched GED views");

    struct overlay_sync_context ctx;
    ctx.view = &view;
    if (ged_diagnostic_line_layer_handler_available(gedp))
	FAIL("fresh GED should not already have a diagnostic line-layer handler");
    ged_diagnostic_line_layer_handler_set(gedp, &test_overlay_handler, &ctx);
    if (!ged_diagnostic_line_layer_handler_available(gedp))
	FAIL("GED diagnostic line-layer handler should register");
    if (!ged_diagnostic_line_layer_publish(gedp, "qtcad::diagnostic", builder))
	FAIL("GED diagnostic line-layer publish should sync into qtcad Obol");
    if (ged_diagnostic_hud_label_handler_available(gedp))
	FAIL("fresh GED should not already have a diagnostic HUD label handler");
    ged_diagnostic_hud_label_handler_set(gedp, &test_hud_label_handler, &ctx);
    if (!ged_diagnostic_hud_label_handler_available(gedp))
	FAIL("GED diagnostic HUD label handler should register");
    struct ged_diagnostic_hud_label label = GED_DIAGNOSTIC_HUD_LABEL_INIT;
    label.label_id = "qtcad::diagnostic::summary";
    label.text = "2 diagnostic segments";
    label.position[0] = 6.0;
    label.position[1] = 22.0;
    label.color[0] = 255;
    label.color[1] = 220;
    label.color[2] = 0;
    label.font_size = 11.0;
    label.source_id = 42;
    if (!ged_diagnostic_hud_label_publish(gedp, &label))
	FAIL("GED diagnostic HUD label publish should sync into qtcad Obol");

    BRLObolViewController *controller = view.obolViewController();
    if (!controller || !controller->getSceneRoot() ||
	    !controller->getSceneRoot()->isOfType(SoGroup::getClassTypeId()))
	FAIL("qtcad Obol overlay sync should have a scene root group");
    SoGroup *root = static_cast<SoGroup *>(controller->getSceneRoot());
    if (root->getNumChildren() != 2)
	FAIL("qtcad Obol overlay sync should add scene line and label overlays");
    SoNode *node = root->getChild(0);
    if (!node || !node->isOfType(SoBRLLineLayerOverlay::getClassTypeId()))
	FAIL("qtcad Obol overlay sync should add a line-layer overlay node");
    SoBRLLineLayerOverlay *overlay = static_cast<SoBRLLineLayerOverlay *>(node);
    if (!BU_STR_EQUAL(overlay->overlayId.getValue().getString(), "qtcad::diagnostic") ||
	    overlay->getLayerShapeCount() != 2 ||
	    overlay->getPointCount() != 4)
	FAIL("qtcad Obol overlay sync should preserve overlay identity and layers");

    SoBRLVListShape *redLayer = overlay->getLayerShape(0);
    SoBRLVListShape *cyanLayer = overlay->getLayerShape(1);
    if (!redLayer || !cyanLayer ||
	    redLayer->getSegmentCount() != 1 ||
	    cyanLayer->getSegmentCount() != 1 ||
	    !BU_STR_EQUAL(redLayer->sourceType.getValue().getString(), "line-layer") ||
	    !BU_STR_EQUAL(redLayer->sourcePath.getValue().getString(),
		"qtcad::diagnostic/rgb_255_000_000") ||
	    !redLayer->selectable.getValue() ||
	    !cyanLayer->selectable.getValue())
	FAIL("qtcad Obol overlay sync should realize selectable color-layer geometry");
    SoNode *labelNode = root->getChild(1);
    if (!labelNode || !labelNode->isOfType(SoBRLHUDLabelOverlay::getClassTypeId()))
	FAIL("qtcad Obol overlay sync should add a HUD label overlay node");
    SoBRLHUDLabelOverlay *hudOverlay =
	static_cast<SoBRLHUDLabelOverlay *>(labelNode);
    SoHUDLabel *hudLabel = hudOverlay->getHUDLabel();
    if (!hudLabel ||
	    !BU_STR_EQUAL(hudOverlay->labelId.getValue().getString(),
		"qtcad::diagnostic::summary") ||
	    hudOverlay->sourceId.getValue() != 42 ||
	    !BU_STR_EQUAL(hudLabel->string[0].getString(),
		"2 diagnostic segments") ||
	    hudLabel->color.getValue()[0] < 0.99f ||
	    hudLabel->color.getValue()[1] < 0.85f ||
	    hudLabel->color.getValue()[2] > 0.01f ||
	    hudLabel->fontSize.getValue() < 10.9f)
	FAIL("qtcad Obol HUD label sync should preserve label fields");

    controller->getViewport()->viewAll();
    controller->requestRender("overlay-sync-visible");
    QImage visibleImage;
    view.get_viewport_image(visibleImage);
    if (visibleImage.isNull() || lit_pixel_count(visibleImage) < 20)
	FAIL("qtcad Obol overlay sync should be visible through qtcad capture");

    ged_diagnostic_line_layer_handler_set(gedp, NULL, NULL);
    if (ged_diagnostic_line_layer_publish(gedp, "qtcad::diagnostic", builder))
	FAIL("GED diagnostic line-layer publish should stop after handler clear");
    ged_diagnostic_hud_label_handler_set(gedp, NULL, NULL);
    if (ged_diagnostic_hud_label_publish(gedp, &label))
	FAIL("GED diagnostic HUD label publish should stop after handler clear");

    bg_line_layer_builder_free(builder);
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
