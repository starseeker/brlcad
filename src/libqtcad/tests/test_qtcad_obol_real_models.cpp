/*        T E S T _ Q T C A D _ O B O L _ R E A L _ M O D E L S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/export_action.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/view_controller.h"
#include "brlobol/vlist_shape.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/str.h"
#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "qtcad/QgObolDrawSync.h"
#include "qtcad/QgObolMeasure.h"
#include "qtcad/QgObolPick.h"
#include "qtcad/QgObolSnap.h"
#include "qtcad/QgView.h"

#include <Inventor/SoViewport.h>

#include <QApplication>
#include <QImage>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

struct model_case {
    const char *file;
    const char *root;
    int bsgDrawMode;
    int obolDrawMode;
    int minWireShapes;
    int minWireSegments;
    int minMeshShapes;
    int minMeshTriangles;
    int exerciseInteractions;
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
model_path(const char *modelFile, char *dbpath, size_t dbpathLen)
{
    const char *brlcadRoot = getenv("BRLCAD_ROOT");
    if (brlcadRoot)
	snprintf(dbpath, dbpathLen, "%s/share/db/%s", brlcadRoot, modelFile);
    else
	snprintf(dbpath, dbpathLen, "share/db/%s", modelFile);

    return bu_file_exists(dbpath, NULL);
}

static int
total_segment_count(SoBRLDatabaseSource *source)
{
    int ret = 0;
    if (!source)
	return ret;

    for (int i = 0; i < source->getRealizedShapeCount(); i++) {
	SoBRLVListShape *shape = source->getRealizedShape(i);
	if (shape)
	    ret += shape->getSegmentCount();
    }
    return ret;
}

static int
total_triangle_count(SoBRLDatabaseSource *source)
{
    int ret = 0;
    if (!source)
	return ret;

    for (int i = 0; i < source->getRealizedMeshCount(); i++) {
	SoBRLMeshShape *mesh = source->getRealizedMesh(i);
	if (mesh)
	    ret += mesh->getTriangleCount();
    }
    return ret;
}

static float
distance_between(const SbVec3f &a, const SbVec3f &b)
{
    return (a - b).length();
}

static int
longest_export_line(const SoBRLExportAction &exportAction,
	SoBRLExportAction::LineRecord &line)
{
    float longest = -FLT_MAX;
    int found = 0;

    for (int i = 0; i < exportAction.getLineCount(); i++) {
	const SoBRLExportAction::LineRecord &candidate = exportAction.getLine(i);
	float length = distance_between(candidate.a, candidate.b);
	if (length > longest) {
	    longest = length;
	    line = candidate;
	    found = 1;
	}
    }

    return found && longest > 0.0f;
}

static int
exercise_real_model_interactions(const struct model_case &testCase,
	QgView &view,
	BRLObolViewController *controller)
{
    std::vector<QgObolPickRecord> picks;
    int pickCount = qg_obol_pick_point(&view,
	    view.width() / 2, view.height() / 2,
	    80.0f, true, picks);
    if (pickCount <= 0 || picks.empty() ||
	    picks[0].path.empty() ||
	    picks[0].sourceName.empty() ||
	    picks[0].primitiveIndex < 0) {
	fprintf(stderr, "%s:%s qtcad Obol pick workflow did not return BRL-CAD identity\n",
		testCase.file, testCase.root);
	return 0;
    }

    SoBRLExportAction exportAction;
    exportAction.apply(controller->getViewport()->getRoot());
    SoBRLExportAction::LineRecord line;
    if (!longest_export_line(exportAction, line)) {
	fprintf(stderr, "%s:%s did not export a measurable Obol line\n",
		testCase.file, testCase.root);
	return 0;
    }

    SbVec3f midpoint = (line.a + line.b) * 0.5f;
    QgObolSnapRecord snap;
    if (!qg_obol_snap_point(&view, midpoint, 0.5f,
	    QgObolSnapRecord::MIDPOINT, snap) ||
	    snap.kind != QgObolSnapRecord::MIDPOINT ||
	    snap.path.empty() ||
	    snap.primitiveIndex < 0 ||
	    snap.distance > 0.001f) {
	fprintf(stderr, "%s:%s qtcad Obol snap workflow did not find a real-model midpoint\n",
		testCase.file, testCase.root);
	return 0;
    }

    SbVec3f measurePoints[2] = {line.a, line.b};
    if (!qg_obol_measure_update_overlay(&view, "real-model::measurement",
	    measurePoints, 2, NULL)) {
	fprintf(stderr, "%s:%s qtcad Obol measure workflow did not publish an overlay\n",
		testCase.file, testCase.root);
	return 0;
    }

    SoBRLExportAction measuredExport;
    measuredExport.apply(controller->getViewport()->getRoot());
    if (measuredExport.getLineCount() <= exportAction.getLineCount()) {
	fprintf(stderr, "%s:%s qtcad Obol measure overlay was not visible to Obol export\n",
		testCase.file, testCase.root);
	return 0;
    }

    if (!qg_obol_measure_clear_overlay(&view, "real-model::measurement")) {
	fprintf(stderr, "%s:%s qtcad Obol measure workflow did not clear its overlay\n",
		testCase.file, testCase.root);
	return 0;
    }

    SoBRLExportAction clearedExport;
    clearedExport.apply(controller->getViewport()->getRoot());
    if (clearedExport.getLineCount() != exportAction.getLineCount()) {
	fprintf(stderr, "%s:%s qtcad Obol measure overlay remained after clear\n",
		testCase.file, testCase.root);
	return 0;
    }

    return 1;
}

static int
sync_draw_case(const struct model_case &testCase)
{
    char dbpath[MAXPATHLEN] = {0};
    if (!model_path(testCase.file, dbpath, sizeof(dbpath))) {
	fprintf(stderr, "missing qtcad Obol workflow model: %s\n", dbpath);
	return 0;
    }

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp) {
	fprintf(stderr, "failed to open qtcad Obol workflow model: %s\n", dbpath);
	return 0;
    }

    QgView view(NULL, QgView_SW, NULL);
    view.resize(220, 170);
    gedp->ged_gvp = view.view();

    BRLObolViewController *controller = view.obolViewController();
    if (!controller) {
	ged_close(gedp);
	return 0;
    }
    controller->clearDatabaseSources();

    struct bsg_appearance_settings settings = BSG_APPEARANCE_SETTINGS_INIT;
    settings.draw_mode = testCase.bsgDrawMode;

    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, testCase.root);
    txn.view = view.view();
    txn.appearance = &settings;

    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int drawRet = ged_draw_apply_transaction(gedp, &txn, &result);
    int changed = qg_obol_sync_draw_transaction(gedp, &txn, &result, &view);
    if (drawRet < 0) {
	fprintf(stderr, "%s:%s GED draw failed: %s\n", testCase.file,
		testCase.root, bu_vls_cstr(&result.errors));
	ged_draw_transaction_result_free(&result);
	ged_close(gedp);
	return 0;
    }
    ged_draw_transaction_result_free(&result);

    if (!changed || controller->getDatabaseSourceCount() != 1) {
	fprintf(stderr, "%s:%s did not create one Obol database source\n",
		testCase.file, testCase.root);
	ged_close(gedp);
	return 0;
    }

    SoBRLDatabaseSource *source = controller->getDatabaseSource(0);
    if (!source ||
	    !BU_STR_EQUAL(source->path.getValue().getString(), testCase.root) ||
	    source->drawMode.getValue() != testCase.obolDrawMode ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED) {
	fprintf(stderr, "%s:%s produced an invalid Obol database source: %s\n",
		testCase.file, testCase.root,
		source ? source->realizationDiagnostic.getValue().getString() : "none");
	ged_close(gedp);
	return 0;
    }

    if (testCase.obolDrawMode == SoBRLDatabaseSource::SHADED) {
	int meshCount = source->getRealizedMeshCount();
	int triangleCount = total_triangle_count(source);
	if (meshCount < testCase.minMeshShapes ||
		triangleCount < testCase.minMeshTriangles) {
	    fprintf(stderr, "%s:%s shaded Obol geometry too small: meshes=%d triangles=%d\n",
		    testCase.file, testCase.root, meshCount, triangleCount);
	    ged_close(gedp);
	    return 0;
	}
    } else {
	int shapeCount = source->getRealizedShapeCount();
	int segmentCount = total_segment_count(source);
	if (shapeCount < testCase.minWireShapes ||
		segmentCount < testCase.minWireSegments) {
	    fprintf(stderr, "%s:%s wire Obol geometry too small: shapes=%d segments=%d\n",
		    testCase.file, testCase.root, shapeCount, segmentCount);
	    ged_close(gedp);
	    return 0;
	}
    }

    controller->getViewport()->viewAll();
    controller->requestRender("real-model-visible");
    QImage visibleImage;
    view.get_viewport_image(visibleImage);
    int litPixels = lit_pixel_count(visibleImage);
    if (visibleImage.isNull() || litPixels < 20) {
	fprintf(stderr, "%s:%s qtcad Obol capture too dark: width=%d height=%d lit=%d\n",
		testCase.file, testCase.root,
		visibleImage.width(), visibleImage.height(), litPixels);
	ged_close(gedp);
	return 0;
    }
    if (controller->isRenderRequested()) {
	fprintf(stderr, "%s:%s qtcad Obol real-model capture did not consume the render request\n",
		testCase.file, testCase.root);
	ged_close(gedp);
	return 0;
    }
    if (view.dmp()) {
	fprintf(stderr, "%s:%s qtcad Obol real-model capture initialized the legacy display manager\n",
		testCase.file, testCase.root);
	ged_close(gedp);
	return 0;
    }

    if (testCase.exerciseInteractions &&
	    !exercise_real_model_interactions(testCase, view, controller)) {
	ged_close(gedp);
	return 0;
    }

    ged_close(gedp);
    return 1;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const struct model_case cases[] = {
	{"pinewood.g", "pinewood", BSG_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 2, 41, 0, 0, 0},
	{"pinewood.g", "pinewood", BSG_DRAW_MODE_SHADED,
	    SoBRLDatabaseSource::SHADED, 0, 0, 21, 501, 0},
	{"havoc.g", "havoc", BSG_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 10, 100, 0, 0, 0},
	{"m35.g", "all.g", BSG_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 100, 1000, 0, 0, 1}
    };

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
	if (!sync_draw_case(cases[i]))
	    FAIL("qtcad Obol real-model draw workflow should pass");
    }

    if (BU_STR_EQUAL(getenv("BRLOBOL_QTCAD_GENERIC_TWIN"), "1")) {
	const struct model_case genericTwinCase = {
	    "faa/Generic_Twin.g", "all", BSG_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 100, 1000, 0, 0, 0
	};
	if (!sync_draw_case(genericTwinCase))
	    FAIL("qtcad Obol Generic_Twin maturity workflow should pass");
    }

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
