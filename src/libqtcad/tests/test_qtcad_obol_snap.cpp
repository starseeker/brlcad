/*              T E S T _ Q T C A D _ O B O L _ S N A P . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

extern "C" {
#include "bsg.h"
#include "bsg/view_state.h"
}

#include "brlobol/view_controller.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "qtcad/QgObolDrawSync.h"
#include "qtcad/QgObolSnap.h"
#include "qtcad/QgView.h"
#include "qtcad/QgViewFilter.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>

#include <QApplication>
#include <QMouseEvent>

#include <math.h>
#include <stdio.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

class SnapProbeFilter : public QgViewFilter {
public:
    bool sync(QEvent *event)
    {
	return view_sync(event) != NULL;
    }
};

static int
nearly_equal(float a, float b)
{
    return fabsf(a - b) < 0.001f;
}

static int
near_point(const SbVec3f &p, float x, float y, float z)
{
    return nearly_equal(p[0], x) && nearly_equal(p[1], y) &&
	nearly_equal(p[2], z);
}

static int
make_snap_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0;

    point_t omin = { 9.0,  9.0,  9.0};
    point_t omax = {11.0, 11.0, 11.0};
    ret = ret && mk_rpp(wdbp, "snap_only.s", omin, omax) == 0;

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
left_move_at(int x, int y)
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    return QMouseEvent(QEvent::MouseMove, QPoint(x, y),
	    Qt::NoButton, Qt::NoButton, Qt::NoModifier);
#else
    return QMouseEvent(QEvent::MouseMove, QPointF(x, y),
	    QPointF(x, y), Qt::NoButton, Qt::NoButton, Qt::NoModifier);
#endif
}

static void
set_center_query(struct bsg_view *v, fastf_t x, fastf_t y, fastf_t z)
{
    v->gv_width = 200;
    v->gv_height = 200;
    v->gv_size = 2.0;
    MAT_IDN(v->gv_view2model);
    MAT_DELTAS(v->gv_view2model, x, y, z);
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_snap_tmp.g";
    if (!make_snap_db(dbpath))
	FAIL("failed to create qtcad Obol snap test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol snap test database");

    QgView view(NULL, QgView_SW, NULL);
    view.resize(200, 200);
    gedp->ged_gvp = view.view();

    BRLObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();
    controller->setViewportSize(200, 200);

    struct ged_draw_transaction draw_box =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "box.s");
    draw_box.view = view.view();
    if (!apply_and_sync(gedp, &view, &draw_box))
	FAIL("GED wire draw should sync box source into Obol");

    QgObolSnapRecord record;
    if (!qg_obol_snap_point(&view, SbVec3f(1.02f, 1.02f, 1.02f),
	    0.1f, QgObolSnapRecord::ENDPOINT, record))
	FAIL("qtcad Obol snap helper should find a box endpoint");
    if (record.path != "/box.s" ||
	    record.kind != QgObolSnapRecord::ENDPOINT ||
	    record.primitiveIndex < 0 ||
	    !near_point(record.point, 1.0f, 1.0f, 1.0f))
	FAIL("qtcad Obol snap helper should preserve endpoint identity");
    if (view.dmp())
	FAIL("qtcad Obol snap helper should not initialize the legacy display manager");

    if (controller->replaceDatabaseSource("snap_only.s", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 2) <= 0 ||
	    !controller->realizePending())
	FAIL("test should add an Obol-only snap source");
    if (view.dmp())
	FAIL("qtcad Obol snap setup should not initialize the legacy display manager");

    struct bsg_view *bv = view.view();
    set_center_query(bv, 11.02, 11.02, 11.02);
    bsg_view_set_snap_source_flags(bv, BSG_SNAP_DB);
    bsg_view_set_snap_lines(bv, 1);

    SnapProbeFilter filter;
    filter.set_view(bv);
    filter.set_view_widget(&view);
    QMouseEvent move = left_move_at(100, 100);
    if (!filter.sync(&move))
	FAIL("qtcad view filter should accept a snap probe mouse event");
    if (!nearly_equal((float)bv->gv_point[X], 11.0f) ||
	    !nearly_equal((float)bv->gv_point[Y], 11.0f) ||
	    !nearly_equal((float)bv->gv_point[Z], 11.0f))
	FAIL("qtcad view filter should refine database snapping through Obol");
    if (view.dmp())
	FAIL("qtcad Obol snap refinement should not initialize the legacy display manager");

    set_center_query(bv, 11.02, 11.02, 11.02);
    bsg_view_set_snap_source_flags(bv, BSG_SNAP_VIEW);
    QMouseEvent viewScopedMove = left_move_at(100, 100);
    if (!filter.sync(&viewScopedMove))
	FAIL("qtcad view filter should accept a view-scoped snap probe event");
    if (nearly_equal((float)bv->gv_point[X], 11.0f) &&
	    nearly_equal((float)bv->gv_point[Y], 11.0f) &&
	    nearly_equal((float)bv->gv_point[Z], 11.0f))
	FAIL("qtcad view filter should leave view-scoped snapping on the BSG path");
    if (view.dmp())
	FAIL("qtcad view-scoped snap fallback should not initialize the legacy display manager");

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
