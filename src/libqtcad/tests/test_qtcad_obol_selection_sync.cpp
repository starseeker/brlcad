/*    T E S T _ Q T C A D _ O B O L _ S E L E C T I O N _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/view_controller.h"
#include "brlobol/vlist_shape.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "ged/selection_state.h"
#include "qtcad/QgLegacyViewBsg.h"
#include "qtcad/QgObolDrawSync.h"
#include "qtcad/QgObolSelectionSync.h"
#include "qtcad/QgView.h"
#include "raytrace.h"
#include "wdb.h"

#include <QApplication>

#include <stdio.h>
#include <string.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static int
make_selection_sync_db(const char *dbpath)
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

static SoBRLDatabaseSource *
find_source(BRLObolViewController *controller, const char *path)
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

static int
source_has_selected_shapes(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;

    for (int i = 0; i < source->getRealizedShapeCount(); i++) {
	SoBRLVListShape *shape = source->getRealizedShape(i);
	if (shape && shape->selected.getValue())
	    return 1;
    }
    for (int i = 0; i < source->getRealizedMeshCount(); i++) {
	SoBRLMeshShape *mesh = source->getRealizedMesh(i);
	if (mesh && mesh->selected.getValue())
	    return 1;
    }
    return 0;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_selection_sync_tmp.g";
    if (!make_selection_sync_db(dbpath))
	FAIL("failed to create qtcad Obol selection-sync test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol selection-sync test database");

    QgView view(NULL, QgView_SW, NULL);
    view.resize(180, 140);
    gedp->ged_gvp = qg_legacy_view_to_bsg(view.view());

    BRLObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();

    struct ged_draw_transaction draw_box =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "box.s");
    draw_box.view = qg_legacy_view_to_bsg(view.view());
    if (!apply_and_sync(gedp, &view, &draw_box))
	FAIL("GED wire draw should sync box source into Obol");

    struct bsg_appearance_settings shaded_settings = BSG_APPEARANCE_SETTINGS_INIT;
    shaded_settings.draw_mode = BSG_DRAW_MODE_SHADED;
    struct ged_draw_transaction draw_ball =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "ball.s");
    draw_ball.view = qg_legacy_view_to_bsg(view.view());
    draw_ball.appearance = &shaded_settings;
    if (!apply_and_sync(gedp, &view, &draw_ball))
	FAIL("GED shaded draw should sync ball source into Obol");

    SoBRLDatabaseSource *box = find_source(controller, "box.s");
    SoBRLDatabaseSource *ball = find_source(controller, "ball.s");
    if (!box || box->getRealizedShapeCount() <= 0)
	FAIL("box should have realized Obol wire geometry");
    if (!ball || ball->getRealizedMeshCount() <= 0)
	FAIL("ball should have realized Obol mesh geometry");

    if (qg_obol_sync_selection_state(gedp, &view, nullptr) != 0)
	FAIL("empty GED selection should not change unselected Obol geometry");
    if (source_has_selected_shapes(box) || source_has_selected_shapes(ball))
	FAIL("fresh Obol geometry should start unselected");

    if (!ged_selection_select_path(gedp, nullptr, "box.s", 1))
	FAIL("GED should select box path");
    if (!qg_obol_sync_selection_state(gedp, &view, nullptr))
	FAIL("box selection should update Obol selection fields");
    if (!source_has_selected_shapes(box) || source_has_selected_shapes(ball))
	FAIL("Obol selection should select only box geometry");
    if (qg_obol_sync_selection_state(gedp, &view, nullptr) != 0)
	FAIL("repeating selection sync should be stable");

    if (!ged_selection_select_path(gedp, nullptr, "ball.s", 1))
	FAIL("GED should select ball path");
    if (!qg_obol_sync_selection_state(gedp, &view, nullptr))
	FAIL("ball selection should update Obol mesh selection fields");
    if (!source_has_selected_shapes(box) || !source_has_selected_shapes(ball))
	FAIL("Obol selection should retain box and add ball geometry");

    if (!ged_selection_clear(gedp, nullptr))
	FAIL("GED should clear default selection");
    ged_selection_recompute(gedp, nullptr);
    if (!qg_obol_sync_selection_state(gedp, &view, nullptr))
	FAIL("clearing GED selection should clear Obol selection fields");
    if (source_has_selected_shapes(box) || source_has_selected_shapes(ball))
	FAIL("Obol selection fields should clear with GED selection");

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
