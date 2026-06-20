/*              T E S T _ Q T C A D _ O B O L _ P I C K . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/view_controller.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "qtcad/QgObolDrawSync.h"
#include "qtcad/QgObolPick.h"
#include "qtcad/QgSelectFilter.h"
#include "qtcad/QgView.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>

#include <QApplication>
#include <QMouseEvent>

#include <stdio.h>
#include <string.h>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static int
make_pick_db(const char *dbpath)
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

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_pick_tmp.g";
    if (!make_pick_db(dbpath))
	FAIL("failed to create qtcad Obol pick test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol pick test database");

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

    std::vector<QgObolPickRecord> picks;
    if (qg_obol_pick_point(&view, 90, 70, 8.0f, false, picks) != 1)
	FAIL("qtcad Obol point pick should find the centered box");
    if (picks[0].path != "/box.s" ||
	    picks[0].primitiveKind != QgObolPickRecord::FACE ||
	    picks[0].primitiveIndex < 0 ||
	    picks[0].sourceName != "box.s" ||
	    picks[0].sourceType != "arb8")
	FAIL("qtcad Obol pick record should preserve BRL-CAD pick detail identity");
    if (view.dmp())
	FAIL("qtcad Obol point pick should not initialize the legacy display manager");

    QgSelectPntFilter filter;
    filter.set_view(view.view());
    filter.set_view_widget(&view);
    QMouseEvent release = mouse_event(QEvent::MouseButtonRelease, 90, 70,
	    Qt::LeftButton, Qt::LeftButton);
    if (!filter.eventFilter(NULL, &release))
	FAIL("qtcad select point filter should accept the release event");
    const std::vector<std::string> &paths = filter.selected_paths();
    if (paths.size() != 1 || paths[0] != "box.s")
	FAIL("qtcad select point filter should expose normalized Obol pick paths");
    if (view.dmp())
	FAIL("qtcad select point filter should not initialize the legacy display manager for Obol picks");

    std::vector<QgObolPickRecord> rectPicks;
    if (qg_obol_pick_rect(&view, 70, 50, 110, 90, 8.0f, false, rectPicks) <= 0)
	FAIL("qtcad Obol rectangle pick should find the centered box");
    bool foundRectPath = false;
    for (const QgObolPickRecord &pick : rectPicks) {
	if (pick.path == "/box.s") {
	    foundRectPath = true;
	    break;
	}
    }
    if (!foundRectPath)
	FAIL("qtcad Obol rectangle pick should preserve BRL-CAD path identity");
    if (view.dmp())
	FAIL("qtcad Obol rectangle pick should not initialize the legacy display manager");

    QgSelectBoxFilter boxFilter;
    boxFilter.set_view(view.view());
    boxFilter.set_view_widget(&view);
    QMouseEvent boxPress = mouse_event(QEvent::MouseButtonPress, 70, 50,
	    Qt::LeftButton, Qt::LeftButton);
    QMouseEvent boxMove = mouse_event(QEvent::MouseMove, 110, 90,
	    Qt::NoButton, Qt::LeftButton);
    QMouseEvent boxRelease = mouse_event(QEvent::MouseButtonRelease, 110, 90,
	    Qt::LeftButton, Qt::LeftButton);
    if (!boxFilter.eventFilter(NULL, &boxPress) ||
	    !boxFilter.eventFilter(NULL, &boxMove) ||
	    !boxFilter.eventFilter(NULL, &boxRelease))
	FAIL("qtcad select box filter should accept the Obol rectangle workflow");
    const std::vector<std::string> &boxPaths = boxFilter.selected_paths();
    if (boxPaths.size() != 1 || boxPaths[0] != "box.s")
	FAIL("qtcad select box filter should expose normalized Obol pick paths");
    if (view.dmp())
	FAIL("qtcad select box filter should not initialize the legacy display manager for Obol picks");

    QgSelectRayFilter rayFilter;
    rayFilter.dbip = gedp->dbip;
    rayFilter.set_view(view.view());
    rayFilter.set_view_widget(&view);
    QMouseEvent rayRelease = mouse_event(QEvent::MouseButtonRelease, 90, 70,
	    Qt::LeftButton, Qt::LeftButton);
    if (!rayFilter.eventFilter(NULL, &rayRelease))
	FAIL("qtcad select ray filter should accept the Obol-backed release event");
    const std::vector<std::string> &rayPaths = rayFilter.selected_paths();
    if (rayPaths.size() != 1 || rayPaths[0] != "box.s")
	FAIL("qtcad select ray filter should expose normalized Obol pick paths");
    if (view.dmp())
	FAIL("qtcad select ray filter should not initialize the legacy display manager for Obol picks");

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
