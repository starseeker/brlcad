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
#include "bu/time.h"
#include "ged.h"
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgObolDrawSync.h"
#include "qtcad/QgObolMeasure.h"
#include "qtcad/QgObolPick.h"
#include "qtcad/QgObolSnap.h"
#include "qtcad/QgView.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoGroup.h>

#include <QApplication>
#include <QImage>

#include <float.h>
#include <math.h>
#include <stdint.h>
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
    const char *name;
    const char *file;
    const char *root;
    int legacyDrawMode;
    int obolDrawMode;
    int minWireShapes;
    int minWireSegments;
    int minMeshShapes;
    int minMeshTriangles;
    int exerciseInteractions;
    int pickAllStress;
};

struct geometry_counts {
    int shapeCount;
    int segmentCount;
    int meshCount;
    int triangleCount;
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
should_run_case(int argc, char **argv, const char *name)
{
    if (argc <= 1)
	return 1;

    for (int i = 1; i < argc; i++) {
	if (BU_STR_EQUAL(argv[i], "all") || BU_STR_EQUAL(argv[i], name))
	    return 1;
    }

    return 0;
}

static int
case_was_requested(int argc, char **argv, const char *name)
{
    for (int i = 1; i < argc; i++) {
	if (BU_STR_EQUAL(argv[i], name))
	    return 1;
    }

    return 0;
}

static int
timing_enabled()
{
    return BU_STR_EQUAL(getenv("BRLOBOL_QTCAD_REAL_MODEL_TIMING"), "1");
}

static void
print_timing(const struct model_case &testCase, const char *phase, int64_t start)
{
    if (!timing_enabled())
	return;
    double elapsed = (double)(bu_gettime() - start) / 1000000.0;
    fprintf(stderr, "TIMING %s %s %.3f sec\n", testCase.name, phase, elapsed);
}

static void
accumulate_geometry_counts(SoNode *node, struct geometry_counts &counts)
{
    if (!node)
	return;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	counts.shapeCount++;
	counts.segmentCount += shape->getSegmentCount();
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *mesh = static_cast<SoBRLMeshShape *>(node);
	counts.meshCount++;
	counts.triangleCount += mesh->getTriangleCount();
	return;
    }

    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    accumulate_geometry_counts(group->getChild(i), counts);
    }
}

static struct geometry_counts
realized_geometry_counts(SoBRLDatabaseSource *source)
{
    struct geometry_counts counts = {0, 0, 0, 0};
    if (!source)
	return counts;

    for (int i = 0; i < source->getNumChildren(); i++)
	accumulate_geometry_counts(source->getChild(i), counts);
    return counts;
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
    int64_t phaseStart = bu_gettime();
    std::vector<QgObolPickRecord> picks;
    bool pickAll = testCase.pickAllStress ? true : false;
    int pickCount = qg_obol_pick_point(&view,
	    view.width() / 2, view.height() / 2,
	    80.0f, pickAll, picks);
    print_timing(testCase, pickAll ? "pick-all-stress" : "pick", phaseStart);
    if (pickCount <= 0 || picks.empty() ||
	    picks[0].path.empty() ||
	    picks[0].sourceName.empty() ||
	    picks[0].primitiveIndex < 0) {
	fprintf(stderr, "%s:%s qtcad Obol pick workflow did not return BRL-CAD identity\n",
		testCase.file, testCase.root);
	return 0;
    }
    if (pickAll) {
	if (pickCount < 2) {
	    fprintf(stderr, "%s:%s qtcad Obol pick-all stress did not return multiple hits\n",
		    testCase.file, testCase.root);
	    return 0;
	}
	for (size_t i = 1; i < picks.size(); i++) {
	    if (picks[i].distance < picks[i - 1].distance) {
		fprintf(stderr, "%s:%s qtcad Obol pick-all stress returned unordered hits\n",
			testCase.file, testCase.root);
		return 0;
	    }
	}
    }

    phaseStart = bu_gettime();
    SoBRLExportAction exportAction;
    exportAction.apply(controller->getViewport()->getRoot());
    SoBRLExportAction::LineRecord line;
    print_timing(testCase, "export-longest-line", phaseStart);
    if (!longest_export_line(exportAction, line)) {
	fprintf(stderr, "%s:%s did not export a measurable Obol line\n",
		testCase.file, testCase.root);
	return 0;
    }

    phaseStart = bu_gettime();
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
    print_timing(testCase, "snap", phaseStart);

    phaseStart = bu_gettime();
    SbVec3f measurePoints[2] = {line.a, line.b};
    if (!qg_obol_measure_update_overlay(&view, "real-model::measurement",
	    measurePoints, 2, NULL)) {
	fprintf(stderr, "%s:%s qtcad Obol measure workflow did not publish an overlay\n",
		testCase.file, testCase.root);
	return 0;
    }

    SoBRLExportAction measuredExport;
    measuredExport.apply(controller->getViewport()->getRoot());
    print_timing(testCase, "measure-update-export", phaseStart);
    if (measuredExport.getLineCount() <= exportAction.getLineCount()) {
	fprintf(stderr, "%s:%s qtcad Obol measure overlay was not visible to Obol export\n",
		testCase.file, testCase.root);
	return 0;
    }

    phaseStart = bu_gettime();
    if (!qg_obol_measure_clear_overlay(&view, "real-model::measurement")) {
	fprintf(stderr, "%s:%s qtcad Obol measure workflow did not clear its overlay\n",
		testCase.file, testCase.root);
	return 0;
    }

    SoBRLExportAction clearedExport;
    clearedExport.apply(controller->getViewport()->getRoot());
    print_timing(testCase, "measure-clear-export", phaseStart);
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
    int64_t totalStart = bu_gettime();
    char dbpath[MAXPATHLEN] = {0};
    if (!model_path(testCase.file, dbpath, sizeof(dbpath))) {
	fprintf(stderr, "missing qtcad Obol workflow model: %s\n", dbpath);
	return 0;
    }

    struct ged *gedp = ged_open("db", dbpath, 1);
    print_timing(testCase, "ged-open", totalStart);
    if (!gedp) {
	fprintf(stderr, "failed to open qtcad Obol workflow model: %s\n", dbpath);
	return 0;
    }

    QgView view(NULL, QgView_SW);
    view.resize(220, 170);
    qg_legacy_view_ged_active_set(gedp, view.view());

    BRLObolViewController *controller = view.obolViewController();
    if (!controller) {
	ged_close(gedp);
	return 0;
    }
    controller->clearDatabaseSources();

    qg_legacy_view_draw_appearance *appearance =
	qg_legacy_view_draw_appearance_create(testCase.legacyDrawMode);

    qg_legacy_view_draw_transaction txn;
    qg_legacy_view_draw_transaction_init(&txn,
	    QG_LEGACY_VIEW_DRAW_TXN_DRAW, testCase.root);
    qg_legacy_view_draw_transaction_view_set(&txn, view.view());
    qg_legacy_view_draw_transaction_appearance_set(&txn, appearance);

    qg_legacy_view_draw_transaction_result result;
    qg_legacy_view_draw_result_init(&result);
    int64_t phaseStart = bu_gettime();
    int drawRet = qg_legacy_view_draw_transaction_apply(gedp, &txn, &result);
    print_timing(testCase, "ged-draw-transaction", phaseStart);
    phaseStart = bu_gettime();
    int changed = qg_obol_sync_draw_transaction(gedp, &txn, &result, &view);
    print_timing(testCase, "obol-sync-transaction", phaseStart);
    qg_legacy_view_draw_appearance_destroy(appearance);
    if (drawRet < 0) {
	const char *drawErrors = qg_legacy_view_draw_result_errors(&result);
	fprintf(stderr, "%s:%s GED draw failed: %s\n", testCase.file,
		testCase.root, drawErrors ? drawErrors : "");
	qg_legacy_view_draw_result_free(&result);
	ged_close(gedp);
	return 0;
    }
    qg_legacy_view_draw_result_free(&result);

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

    phaseStart = bu_gettime();
    struct geometry_counts counts = realized_geometry_counts(source);
    if (testCase.obolDrawMode == SoBRLDatabaseSource::SHADED) {
	if (counts.meshCount < testCase.minMeshShapes ||
		counts.triangleCount < testCase.minMeshTriangles) {
	    fprintf(stderr, "%s:%s shaded Obol geometry too small: meshes=%d triangles=%d\n",
		    testCase.file, testCase.root, counts.meshCount,
		    counts.triangleCount);
	    ged_close(gedp);
	    return 0;
	}
    } else {
	if (counts.shapeCount < testCase.minWireShapes ||
		counts.segmentCount < testCase.minWireSegments) {
	    fprintf(stderr, "%s:%s wire Obol geometry too small: shapes=%d segments=%d\n",
		    testCase.file, testCase.root, counts.shapeCount,
		    counts.segmentCount);
	    ged_close(gedp);
	    return 0;
	}
    }
    print_timing(testCase, "geometry-count-check", phaseStart);

    controller->getViewport()->viewAll();
    controller->requestRender("real-model-visible");
    phaseStart = bu_gettime();
    QImage visibleImage;
    view.get_viewport_image(visibleImage);
    int litPixels = lit_pixel_count(visibleImage);
    print_timing(testCase, "render-readback", phaseStart);
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
    if (view.legacyBackendInitialized()) {
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

    print_timing(testCase, "total", totalStart);
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
	{"pinewood_wire", "pinewood.g", "pinewood", QG_LEGACY_VIEW_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 2, 41, 0, 0, 0, 0},
	{"pinewood_shaded", "pinewood.g", "pinewood", QG_LEGACY_VIEW_DRAW_MODE_SHADED,
	    SoBRLDatabaseSource::SHADED, 0, 0, 21, 501, 0, 0},
	{"havoc_wire", "havoc.g", "havoc", QG_LEGACY_VIEW_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 10, 100, 0, 0, 0, 0},
	{"m35_wire_interactions", "m35.g", "all.g", QG_LEGACY_VIEW_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 100, 1000, 0, 0, 1, 0},
	{"m35_wire_pick_all_stress", "m35.g", "all.g", QG_LEGACY_VIEW_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 100, 1000, 0, 0, 1, 1}
    };

    int ran = 0;
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
	if (!should_run_case(argc, argv, cases[i].name))
	    continue;
	ran = 1;
	if (!sync_draw_case(cases[i]))
	    FAIL("qtcad Obol real-model draw workflow should pass");
    }

    if (BU_STR_EQUAL(getenv("BRLOBOL_QTCAD_GENERIC_TWIN"), "1")) {
	const struct model_case genericTwinCase = {
	    "generic_twin_wire", "faa/Generic_Twin.g", "all", QG_LEGACY_VIEW_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 100, 1000, 0, 0, 0, 0
	};
	if (should_run_case(argc, argv, genericTwinCase.name)) {
	    ran = 1;
	    if (!sync_draw_case(genericTwinCase))
		FAIL("qtcad Obol Generic_Twin maturity workflow should pass");
	}
    } else if (case_was_requested(argc, argv, "generic_twin_wire")) {
	    FAIL("qtcad Obol Generic_Twin maturity workflow should pass");
    }

    if (!ran)
	FAIL("unknown qtcad Obol real-model case requested");

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
