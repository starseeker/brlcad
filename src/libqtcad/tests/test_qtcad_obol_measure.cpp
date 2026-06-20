/*          T E S T _ Q T C A D _ O B O L _ M E A S U R E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/line_layer_overlay.h"
#include "brlobol/view_controller.h"
#include "bu/app.h"
#include "bu/color.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/str.h"
#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "qtcad/QgMeasureFilter.h"
#include "qtcad/QgObolDrawSync.h"
#include "qtcad/QgObolMeasure.h"
#include "qtcad/QgView.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

#include <QApplication>
#include <QMouseEvent>

#include <stdio.h>
#include <string>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static int
make_measure_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0;
    wdb_close(wdbp);
    return ret;
}

static int
apply_and_sync(struct ged *gedp,
	QgView *view,
	struct ged_draw_transaction *txn)
{
    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int draw_ret = ged_draw_apply_transaction(gedp, txn, &result);
    int changed = qg_obol_sync_draw_transaction(gedp, txn, &result, view);
    ged_draw_transaction_result_free(&result);

    return draw_ret >= 0 && changed != 0;
}

static QMouseEvent
mouse_event(QEvent::Type type,
	int x,
	int y,
	Qt::MouseButton button,
	Qt::MouseButtons buttons)
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    return QMouseEvent(type, QPoint(x, y), button, buttons, Qt::NoModifier);
#else
    return QMouseEvent(type, QPointF(x, y), QPointF(x, y), button, buttons,
	    Qt::NoModifier);
#endif
}

static SoBRLLineLayerOverlay *
find_measure_overlay(BRLObolViewController *controller,
	const char *overlayId)
{
    if (!controller || !controller->getSceneRoot() ||
	    !controller->getSceneRoot()->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    SoGroup *root = static_cast<SoGroup *>(controller->getSceneRoot());
    for (int i = 0; i < root->getNumChildren(); i++) {
	SoNode *node = root->getChild(i);
	if (!node || !node->isOfType(SoBRLLineLayerOverlay::getClassTypeId()))
	    continue;
	SoBRLLineLayerOverlay *overlay =
	    static_cast<SoBRLLineLayerOverlay *>(node);
	if (BU_STR_EQUAL(overlay->overlayId.getValue().getString(), overlayId))
	    return overlay;
    }
    return NULL;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_measure_tmp.g";
    if (!make_measure_db(dbpath))
	FAIL("failed to create qtcad Obol measure test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol measure test database");

    QgView view(NULL, QgView_SW, NULL);
    view.resize(180, 140);
    gedp->ged_gvp = view.view();

    BRLObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();
    controller->setViewportSize(180, 140);

    struct bsg_appearance_settings shaded_settings = BSG_APPEARANCE_SETTINGS_INIT;
    shaded_settings.draw_mode = BSG_DRAW_MODE_SHADED;
    struct ged_draw_transaction draw_box =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "box.s");
    draw_box.view = view.view();
    draw_box.appearance = &shaded_settings;
    if (!apply_and_sync(gedp, &view, &draw_box))
	FAIL("GED shaded draw should sync box source into Obol");

    controller->getViewport()->viewAll();

    SbVec3f pickedPoint;
    std::string pickedPath;
    if (!qg_obol_measure_pick_point(&view, 90, 70, pickedPoint,
	    &pickedPath))
	FAIL("qtcad Obol measure helper should pick the centered box");
    if (pickedPath != "/box.s")
	FAIL("qtcad Obol measure helper should preserve picked path identity");
    if (view.dmp())
	FAIL("qtcad Obol measure pick should not initialize the legacy display manager");

    SbVec3f guidePoints[2] = {
	pickedPoint,
	SbVec3f(pickedPoint[0] + 0.5f, pickedPoint[1], pickedPoint[2])
    };
    struct bu_color green = BU_COLOR_GREEN;
    if (!qg_obol_measure_update_overlay(&view, "tool:measurement",
	    guidePoints, 2, &green))
	FAIL("qtcad Obol measure helper should publish an Obol overlay");
    if (view.dmp())
	FAIL("qtcad Obol measure overlay should not initialize the legacy display manager");

    SoBRLLineLayerOverlay *overlay =
	find_measure_overlay(controller, "tool:measurement");
    if (!overlay || overlay->getLayerShapeCount() != 1 ||
	    overlay->getPointCount() != 2 ||
	    overlay->selectable.getValue())
	FAIL("qtcad Obol measure overlay should be a non-selectable line layer");

    if (!qg_obol_measure_clear_overlay(&view, "tool:measurement") ||
	    find_measure_overlay(controller, "tool:measurement"))
	FAIL("qtcad Obol measure helper should clear its overlay");
    if (view.dmp())
	FAIL("qtcad Obol measure clear should not initialize the legacy display manager");

    QMeasure3DFilter filter;
    filter.set_view(view.view());
    filter.set_view_widget(&view);
    filter.update_color(&green);

    QMouseEvent press = mouse_event(QEvent::MouseButtonPress, 90, 70,
	    Qt::LeftButton, Qt::LeftButton);
    if (!filter.eventFilter(NULL, &press) || filter.mode != 1)
	FAIL("qtcad 3D measure filter should start with an Obol hit point");
    overlay = find_measure_overlay(controller, "tool:measurement");
    if (!overlay || overlay->getPointCount() != 1)
	FAIL("qtcad 3D measure filter should start with a single-point Obol overlay");
    if (view.dmp())
	FAIL("qtcad 3D measure filter start should not initialize the legacy display manager");

    QMouseEvent move = mouse_event(QEvent::MouseMove, 115, 70,
	    Qt::NoButton, Qt::LeftButton);
    if (!filter.eventFilter(NULL, &move))
	FAIL("qtcad 3D measure filter should update the Obol guide line");
    if (view.dmp())
	FAIL("qtcad 3D measure filter update should not initialize the legacy display manager");

    QMouseEvent release = mouse_event(QEvent::MouseButtonRelease, 115, 70,
	    Qt::LeftButton, Qt::LeftButton);
    if (!filter.eventFilter(NULL, &release) || filter.mode != 2)
	FAIL("qtcad 3D measure filter should finalize the first Obol measure segment");
    if (filter.length1() <= 0.0)
	FAIL("qtcad 3D Obol measure segment should report a positive length");

    overlay = find_measure_overlay(controller, "tool:measurement");
    if (!overlay || overlay->getPointCount() != 2)
	FAIL("qtcad 3D measure filter should publish the finalized segment as an Obol overlay");
    if (view.dmp())
	FAIL("qtcad 3D measure filter finalize should not initialize the legacy display manager");

    QMouseEvent rightPress = mouse_event(QEvent::MouseButtonPress, 90, 70,
	    Qt::RightButton, Qt::RightButton);
    if (!filter.eventFilter(NULL, &rightPress) ||
	    find_measure_overlay(controller, "tool:measurement"))
	FAIL("qtcad 3D measure filter should clear the Obol overlay on right click");
    if (view.dmp())
	FAIL("qtcad 3D measure filter clear should not initialize the legacy display manager");

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
