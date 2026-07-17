/*        T E S T _ Q T C A D _ O B O L _ R E A L _ M O D E L S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BExportAction.h"
#include "BObol/BMeshShape.h"
#include "BObol/BViewController.h"
#include "BObol/BVListShape.h"
#include "bu/color.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/str.h"
#include "bu/time.h"
#include "ged.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "icv.h"
#include "QgObolDrawSyncPrivate.h"
#include "qtcad/QgObolMeasure.h"
#include "qtcad/QgObolPick.h"
#include "qtcad/QgObolSnap.h"
#include "qtcad/QgView.h"
#include "rt/db_attr.h"
#include "rt/db_io.h"
#include "rt/mater.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoGroup.h>

#include <QApplication>
#include <QCoreApplication>
#include <QImage>

#include <float.h>
#include <fstream>
#include <math.h>
#include <stdint.h>
#include <set>
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
    int gedDrawMode;
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
image_byte_diff(const QImage &a, const QImage &b)
{
    if (a.size() != b.size())
	return -1;
    QImage ar = a.convertToFormat(QImage::Format_RGBA8888);
    QImage br = b.convertToFormat(QImage::Format_RGBA8888);
    int different = 0;
    for (int y = 0; y < ar.height(); y++) {
	const unsigned char *ap = ar.constScanLine(y);
	const unsigned char *bp = br.constScanLine(y);
	for (int x = 0; x < ar.width() * 4; x++)
	    different += ap[x] != bp[x];
    }
    return different;
}

static fastf_t
image_ssim(const QImage &a, const QImage &b)
{
    if (a.size() != b.size())
	return -1.0;
    QImage ar = a.convertToFormat(QImage::Format_RGBA8888);
    QImage br = b.convertToFormat(QImage::Format_RGBA8888);
    icv_image_t *ai = icv_create((size_t)ar.width(),
	(size_t)ar.height(), ICV_COLOR_SPACE_RGB);
    icv_image_t *bi = icv_create((size_t)br.width(),
	(size_t)br.height(), ICV_COLOR_SPACE_RGB);
    if (!ai || !bi) {
	if (ai)
	    icv_destroy(ai);
	if (bi)
	    icv_destroy(bi);
	return -1.0;
    }
    for (int y = 0; y < ar.height(); y++) {
	const unsigned char *ap = ar.constScanLine(y);
	const unsigned char *bp = br.constScanLine(y);
	for (int x = 0; x < ar.width(); x++) {
	    size_t offset = ((size_t)y * (size_t)ar.width() +
		(size_t)x) * 3;
	    for (int c = 0; c < 3; c++) {
		ai->data[offset + (size_t)c] = ap[x * 4 + c] / 255.0;
		bi->data[offset + (size_t)c] = bp[x * 4 + c] / 255.0;
	    }
	}
    }
    fastf_t score = icv_adiff(ai, bi, ICV_DIFF_SSIM);
    icv_destroy(ai);
    icv_destroy(bi);
    return score;
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
    return BU_STR_EQUAL(getenv("BOBOL_QTCAD_REAL_MODEL_TIMING"), "1");
}

static int
system_gl_enabled()
{
    return BU_STR_EQUAL(getenv("BOBOL_QTCAD_REAL_MODEL_GL"), "1");
}

static BObolViewController::SoftwareWireMode
software_wire_mode()
{
    const char *value = getenv("BOBOL_QTCAD_SOFTWARE_WIRE");
    if (BU_STR_EQUAL(value, "fast"))
	return BObolViewController::SOFTWARE_WIRE_FAST;
    if (BU_STR_EQUAL(value, "quality"))
	return BObolViewController::SOFTWARE_WIRE_QUALITY;
    return BObolViewController::SOFTWARE_WIRE_AUTO;
}

static const char *
software_wire_mode_name(BObolViewController::SoftwareWireMode mode)
{
    switch (mode) {
	case BObolViewController::SOFTWARE_WIRE_FAST:
	    return "fast";
	case BObolViewController::SOFTWARE_WIRE_QUALITY:
	    return "quality";
	default:
	    return "auto";
    }
}

static void
print_timing(const struct model_case &testCase, const char *phase, int64_t start)
{
    if (!timing_enabled())
	return;
    double elapsed = (double)(bu_gettime() - start) / 1000000.0;
    fprintf(stderr, "TIMING %s %s %.3f sec\n", testCase.name, phase, elapsed);
}

static int
copy_file(const char *src_path, const char *dst_path)
{
    if (!src_path || !dst_path)
	return 0;

    std::ifstream src(src_path, std::ios::binary);
    std::ofstream dst(dst_path, std::ios::binary | std::ios::trunc);
    if (!src || !dst)
	return 0;

    dst << src.rdbuf();
    return !src.bad() && dst.good();
}

static int
path_has_component_suffix(const char *path, const char *suffix)
{
    if (!path || !suffix || !suffix[0])
	return 0;

    auto split_components = [](const char *str) {
	std::vector<std::string> components;
	std::string cur;
	while (str && *str) {
	    if (*str == '/') {
		if (!cur.empty()) {
		    components.push_back(cur);
		    cur.clear();
		}
	    } else {
		cur.push_back(*str);
	    }
	    str++;
	}
	if (!cur.empty())
	    components.push_back(cur);
	return components;
    };

    auto strip_instance = [](const std::string &component) {
	size_t atPos = component.rfind('@');
	if (atPos == std::string::npos || atPos == 0 ||
		atPos + 1 >= component.size())
	    return component;
	for (size_t i = atPos + 1; i < component.size(); i++) {
	    if (component[i] < '0' || component[i] > '9')
		return component;
	}
	return component.substr(0, atPos);
    };

    std::vector<std::string> pathComponents = split_components(path);
    std::vector<std::string> suffixComponents = split_components(suffix);
    if (suffixComponents.empty() ||
	    pathComponents.size() < suffixComponents.size())
	return 0;

    size_t offset = pathComponents.size() - suffixComponents.size();
    for (size_t i = 0; i < suffixComponents.size(); i++) {
	if (strip_instance(pathComponents[offset + i]) !=
		strip_instance(suffixComponents[i]))
	    return 0;
    }
    return 1;
}

static int
summary_from_compact_instance(SoBRLDatabaseSource *source,
	const BObolCompactInstanceSummary &instance,
	BObolDatabaseSourceSummary &summary)
{
    if (!source || !instance.valid || !source->getSummary(summary) ||
	!summary.valid)
	return 0;

    summary.path = instance.path;
    summary.materialColorValid = instance.materialColorValid;
    summary.materialColor = instance.materialColor;
    summary.databaseMetadataValid = TRUE;
    summary.databaseRegionId = instance.regionId;
    summary.databaseAirCode = instance.airCode;
    summary.databaseMaterialId = instance.materialId;
    summary.databaseLos = instance.los;
    summary.databaseMaterialColorValid = instance.materialColorValid;
    summary.databaseMaterialColor = instance.materialColor;
    return 1;
}

static int
source_material_matches_rgb(const BObolDatabaseSourceSummary &summary,
	const unsigned char rgb[3])
{
    if (!summary.valid || !summary.materialColorValid)
	return 0;

    return fabsf(summary.materialColor[0] - ((float)rgb[0] / 255.0f)) < 1.0e-5f &&
	fabsf(summary.materialColor[1] - ((float)rgb[1] / 255.0f)) < 1.0e-5f &&
	fabsf(summary.materialColor[2] - ((float)rgb[2] / 255.0f)) < 1.0e-5f;
}

static int
source_material_matches_db_color(struct ged *gedp,
	const BObolDatabaseSourceSummary &summary,
	unsigned char expectedRgb[3])
{
    if (expectedRgb) {
	expectedRgb[0] = 0;
	expectedRgb[1] = 0;
	expectedRgb[2] = 0;
    }
    if (!gedp || !gedp->dbip || !summary.valid)
	return 0;

    const char *path = summary.path.getString();
    while (path && *path == '/')
	path++;
    if (!path || !path[0])
	return 0;

    struct db_full_path fp;
    db_full_path_init(&fp);
    if (db_string_to_path(&fp, gedp->dbip, path) != 0) {
	db_free_full_path(&fp);
	return 0;
    }

    struct bu_color color = BU_COLOR_INIT_ZERO;
    unsigned char rgb[3] = {0, 0, 0};
    db_full_path_color(&color, &fp, gedp->dbip);
    int ret = bu_color_to_rgb_chars(&color, rgb);
    db_free_full_path(&fp);
    if (!ret)
	return 0;
    if (expectedRgb) {
	expectedRgb[0] = rgb[0];
	expectedRgb[1] = rgb[1];
	expectedRgb[2] = rgb[2];
    }
    return source_material_matches_rgb(summary, rgb);
}

static int
all_source_materials_match_db_colors(struct ged *gedp,
	SoBRLSceneController *controller)
{
    if (!gedp || !gedp->dbip || !controller)
	return 0;

    const int sourceCount = controller->getDatabaseSourceCount();
    int materialCount = 0;
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	BObolDatabaseSourceSummary summary;
	if (!source || !source->getSummary(summary) || !summary.valid)
	    return 0;
	if (source->hasCompactInstanceIndex()) {
	    const int instanceCount = source->getCompactInstanceCount();
	    for (int j = 0; j < instanceCount; j++) {
		BObolCompactInstanceHandle handle;
		BObolCompactInstanceSummary instance;
		if (!source->getCompactInstanceHandle(j, handle) ||
		    !source->getCompactInstanceSummary(handle, instance) ||
		    !summary_from_compact_instance(source, instance, summary))
		    return 0;
		SbColor expected;
		if (!bobol_database_source_path_material_color(gedp->dbip,
			summary.path.getString(), expected) ||
		    !summary.materialColorValid ||
		    fabsf(summary.materialColor[0] - expected[0]) >= 1.0e-5f ||
		    fabsf(summary.materialColor[1] - expected[1]) >= 1.0e-5f ||
		    fabsf(summary.materialColor[2] - expected[2]) >= 1.0e-5f) {
		    fprintf(stderr,
			"material sweep/reference mismatch: source=%d path=%s "
			"valid=%d actual=(%.9g %.9g %.9g) expected=(%.9g %.9g %.9g)\n",
			i, summary.path.getString(),
			(int)summary.materialColorValid,
			summary.materialColor[0], summary.materialColor[1],
			summary.materialColor[2], expected[0], expected[1], expected[2]);
		    return 0;
		}
		materialCount++;
	    }
	    continue;
	}

	SbColor expected;
	if (!bobol_database_source_path_material_color(gedp->dbip,
		summary.path.getString(), expected) ||
	    !summary.materialColorValid ||
	    fabsf(summary.materialColor[0] - expected[0]) >= 1.0e-5f ||
	    fabsf(summary.materialColor[1] - expected[1]) >= 1.0e-5f ||
	    fabsf(summary.materialColor[2] - expected[2]) >= 1.0e-5f) {
	    fprintf(stderr,
		"material sweep/reference mismatch: source=%d path=%s "
		"valid=%d actual=(%.9g %.9g %.9g) expected=(%.9g %.9g %.9g)\n",
		i, summary.path.getString(),
		(int)summary.materialColorValid,
		summary.materialColor[0], summary.materialColor[1],
		summary.materialColor[2], expected[0], expected[1], expected[2]);
	    return 0;
	}
	materialCount++;
    }
    return materialCount > 0;
}

static SoBRLDatabaseSource *
find_source_by_path_suffix(SoBRLSceneController *controller,
	const char *suffix,
	BObolDatabaseSourceSummary &summary)
{
    summary = BObolDatabaseSourceSummary();
    if (!controller || !suffix)
	return NULL;

    const int sourceCount = controller->getDatabaseSourceCount();
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (!source || !source->getSummary(summary) || !summary.valid)
	    continue;
	if (source->hasCompactInstanceIndex()) {
	    const int instanceCount = source->getCompactInstanceCount();
	    for (int j = 0; j < instanceCount; j++) {
		BObolCompactInstanceHandle handle;
		BObolCompactInstanceSummary instance;
		if (!source->getCompactInstanceHandle(j, handle) ||
		    !source->getCompactInstanceSummary(handle, instance) ||
		    !summary_from_compact_instance(source, instance, summary))
		    continue;
		if (path_has_component_suffix(summary.path.getString(), suffix))
		    return source;
	    }
	    continue;
	}
	if (path_has_component_suffix(summary.path.getString(), suffix))
	    return source;
    }

    summary = BObolDatabaseSourceSummary();
    return NULL;
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

static struct geometry_counts
realized_geometry_counts(BObolViewController *controller,
			 int expectedDrawMode,
			 int *realizedSources,
			 int *modeMismatches)
{
    struct geometry_counts counts = {0, 0, 0, 0};
    if (realizedSources)
	*realizedSources = 0;
    if (modeMismatches)
	*modeMismatches = 0;
    if (!controller)
	return counts;

    const int sourceCount = controller->getDatabaseSourceCount();
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (!source)
	    continue;
	if (source->drawMode.getValue() != expectedDrawMode) {
	    if (modeMismatches)
		(*modeMismatches)++;
	    continue;
	}
	if (source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED)
	    continue;
	if (realizedSources)
	    (*realizedSources)++;
	struct geometry_counts sourceCounts = realized_geometry_counts(source);
	counts.shapeCount += sourceCounts.shapeCount;
	counts.segmentCount += sourceCounts.segmentCount;
	counts.meshCount += sourceCounts.meshCount;
	counts.triangleCount += sourceCounts.triangleCount;
    }

    /* Retained batches deliberately decouple source records from emitted
     * shape nodes.  Validate their public export records when geometry is no
     * longer parented directly below each database source. */
    if (!counts.segmentCount && !counts.triangleCount &&
	controller->getViewport() && controller->getViewport()->getRoot()) {
	SoBRLExportAction export_action;
	export_action.apply(controller->getViewport()->getRoot());
	std::set<std::string> line_paths;
	std::set<std::string> triangle_paths;
	for (int i = 0; i < export_action.getLineCount(); i++) {
	    const SoBRLExportAction::LineRecord &record = export_action.getLine(i);
	    line_paths.insert(record.path.getString());
	}
	for (int i = 0; i < export_action.getTriangleCount(); i++) {
	    const SoBRLExportAction::TriangleRecord &record =
		export_action.getTriangle(i);
	    triangle_paths.insert(record.path.getString());
	}
	counts.shapeCount = static_cast<int>(line_paths.size());
	counts.segmentCount = export_action.getLineCount();
	counts.meshCount = static_cast<int>(triangle_paths.size());
	counts.triangleCount = export_action.getTriangleCount();
    }

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
	BObolViewController *controller)
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

    QgView view(NULL, system_gl_enabled() ? QgViewType::GL : QgViewType::SW);
    view.resize(220, 170);
    if (system_gl_enabled()) {
	view.show();
	QCoreApplication::processEvents();
    }
	ged_view_active_ctx_set(gedp, view.viewContext());
	(void)ged_view_context_host_attach(gedp, view.viewContext());
    if (!ged_view_context_display_endpoint_set(
	    view.viewContext(), view.displayEndpoint(), 0)) {
	ged_close(gedp);
	return 0;
    }

    BObolViewController *controller = view.obolViewController();
    if (!controller) {
	ged_close(gedp);
	return 0;
    }
    const BObolViewController::SoftwareWireMode wireMode =
	software_wire_mode();
    controller->setSoftwareWireMode(wireMode);
    if (timing_enabled()) {
	fprintf(stderr, "CONFIG %s renderer=%s software_wire=%s size=220x170\n",
	    testCase.name, system_gl_enabled() ? "system-gl" : "osmesa",
	    software_wire_mode_name(wireMode));
    }
    controller->setViewportSize(220, 170);
    controller->clearDatabaseSources();
    controller->requestRender("real-model-empty-baseline");
    QCoreApplication::processEvents();
    QImage emptyImage;
    view.get_viewport_image(emptyImage);

    struct ged_draw_appearance_settings appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    appearance.draw_mode = testCase.gedDrawMode;

    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, testCase.root);
    txn.view = view.viewContext();
    txn.appearance = &appearance;

    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int64_t phaseStart = bu_gettime();
    int drawRet = ged_draw_apply_transaction(gedp, &txn, &result);
    print_timing(testCase, "ged-draw-transaction", phaseStart);
    phaseStart = bu_gettime();
    int changed = qg_obol_sync_ged_draw_transaction(gedp, &txn, &result, &view);
    print_timing(testCase, "obol-sync-transaction", phaseStart);
    if (drawRet < 0) {
	const char *drawErrors = bu_vls_cstr(&result.errors);
	fprintf(stderr, "%s:%s GED draw failed: %s\n", testCase.file,
		testCase.root, drawErrors ? drawErrors : "");
	ged_draw_transaction_result_free(&result);
	ged_close(gedp);
	return 0;
    }
    ged_draw_transaction_result_free(&result);

    /* QgView endpoints use the same deferred publication boundary as the
     * interactive canvas.  Drain it explicitly before inspecting geometry. */
    (void)controller->realizePending();
    BObolViewController *geometryController =
	ged_draw_obol_controller(gedp);
    if (!geometryController)
	geometryController = controller;
    (void)geometryController->realizePending();
    const int sourceCount = geometryController->getDatabaseSourceCount();
    if (!changed || sourceCount <= 0) {
	fprintf(stderr, "%s:%s did not create Obol database sources\n",
		testCase.file, testCase.root);
	ged_close(gedp);
	return 0;
    }

    int realizedSources = 0;
    int modeMismatches = 0;
    phaseStart = bu_gettime();
    struct geometry_counts counts = realized_geometry_counts(geometryController,
	testCase.obolDrawMode, &realizedSources, &modeMismatches);
    if (realizedSources <= 0 || modeMismatches > 0) {
	fprintf(stderr,
		"%s:%s produced invalid Obol database sources: sources=%d realized=%d mode_mismatches=%d expected_mode=%d\n",
		testCase.file, testCase.root, sourceCount, realizedSources,
		modeMismatches, testCase.obolDrawMode);
	ged_close(gedp);
	return 0;
    }

    if (testCase.obolDrawMode == SoBRLDatabaseSource::SHADED) {
	if (counts.meshCount < testCase.minMeshShapes ||
		counts.triangleCount < testCase.minMeshTriangles) {
	    fprintf(stderr,
		    "%s:%s shaded Obol geometry too small: sources=%d realized=%d meshes=%d triangles=%d\n",
		    testCase.file, testCase.root, sourceCount, realizedSources,
		    counts.meshCount, counts.triangleCount);
	    ged_close(gedp);
	    return 0;
	}
    } else {
	if (counts.shapeCount < testCase.minWireShapes ||
		counts.segmentCount < testCase.minWireSegments) {
	    fprintf(stderr,
		    "%s:%s wire Obol geometry too small: sources=%d realized=%d shapes=%d segments=%d\n",
		    testCase.file, testCase.root, sourceCount, realizedSources,
		    counts.shapeCount, counts.segmentCount);
	    ged_close(gedp);
	    return 0;
	}
    }
    print_timing(testCase, "geometry-count-check", phaseStart);

    controller->getViewport()->viewAll();
    controller->requestRender("real-model-visible");
    QCoreApplication::processEvents();
    phaseStart = bu_gettime();
    QImage visibleImage;
    view.get_viewport_image(visibleImage);
    int litPixels = lit_pixel_count(visibleImage);
    QImage comparisonEmpty = emptyImage.size() == visibleImage.size() ?
	emptyImage : emptyImage.scaled(visibleImage.size(), Qt::IgnoreAspectRatio,
	    Qt::SmoothTransformation);
    int visibleByteDiff = image_byte_diff(comparisonEmpty, visibleImage);
    fastf_t visibleSsim = image_ssim(comparisonEmpty, visibleImage);
    print_timing(testCase, "render-readback", phaseStart);
    if (visibleImage.isNull() || litPixels < 20 || visibleByteDiff < 100 ||
	visibleSsim >= 0.9999) {
	fprintf(stderr, "%s:%s qtcad Obol capture did not show geometry: width=%d height=%d lit=%d diff=%d ssim=%.9g\n",
		testCase.file, testCase.root,
		visibleImage.width(), visibleImage.height(), litPixels,
		visibleByteDiff, visibleSsim);
	ged_close(gedp);
	return 0;
    }
    if (controller->isRenderRequested()) {
	fprintf(stderr, "%s:%s qtcad Obol real-model capture did not consume the render request\n",
		testCase.file, testCase.root);
	ged_close(gedp);
	return 0;
    }

    struct ged_draw_transaction redrawTxn =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    redrawTxn.view = view.viewContext();
    struct ged_draw_transaction_result redrawResult;
    ged_draw_transaction_result_init(&redrawResult);
    phaseStart = bu_gettime();
    int redrawRet = ged_draw_apply_transaction(gedp, &redrawTxn,
	&redrawResult);
    (void)qg_obol_sync_ged_draw_transaction(gedp, &redrawTxn,
	&redrawResult, &view);
    int64_t redrawUs = bu_gettime() - phaseStart;
    print_timing(testCase, "redraw-transaction", phaseStart);
    ged_draw_transaction_result_free(&redrawResult);
    QCoreApplication::processEvents();
    QImage redrawnImage;
    view.get_viewport_image(redrawnImage);

    struct ged_draw_transaction_result secondRedrawResult;
    ged_draw_transaction_result_init(&secondRedrawResult);
    int64_t secondRedrawStart = bu_gettime();
    int secondRedrawRet = ged_draw_apply_transaction(gedp, &redrawTxn,
	&secondRedrawResult);
    (void)qg_obol_sync_ged_draw_transaction(gedp, &redrawTxn,
	&secondRedrawResult, &view);
    int64_t secondRedrawUs = bu_gettime() - secondRedrawStart;
    ged_draw_transaction_result_free(&secondRedrawResult);
    QCoreApplication::processEvents();
    QImage secondRedrawnImage;
    view.get_viewport_image(secondRedrawnImage);

    fastf_t redrawSsim = image_ssim(redrawnImage, secondRedrawnImage);
    int redrawByteDiff = image_byte_diff(redrawnImage, secondRedrawnImage);
    int redrawLitPixels = lit_pixel_count(secondRedrawnImage);
    int maxRasterDiff = redrawnImage.width() * redrawnImage.height() * 4 / 100;
    struct geometry_counts redrawCounts = realized_geometry_counts(geometryController,
	testCase.obolDrawMode, NULL, NULL);
    if (redrawRet < 0 || secondRedrawRet < 0 ||
	redrawUs > 10000000 || secondRedrawUs > 10000000 ||
	geometryController->getDatabaseSourceCount() != sourceCount ||
	redrawByteDiff < 0 || redrawByteDiff > maxRasterDiff ||
	redrawCounts.shapeCount != counts.shapeCount ||
	redrawCounts.segmentCount != counts.segmentCount ||
	redrawCounts.meshCount != counts.meshCount ||
	redrawCounts.triangleCount != counts.triangleCount) {
	fprintf(stderr,
		"%s:%s retained Obol redraw failed: ret=%d/%d elapsed=%.3f/%.3f sec sources=%d/%d image_diff=%d/%d ssim=%.9g lit=%d/%d geometry=%d,%d,%d,%d/%d,%d,%d,%d\n",
		testCase.file, testCase.root, redrawRet, secondRedrawRet,
		(double)redrawUs / 1000000.0,
		(double)secondRedrawUs / 1000000.0,
	geometryController->getDatabaseSourceCount(), sourceCount,
		redrawByteDiff, maxRasterDiff, redrawSsim,
		redrawLitPixels, lit_pixel_count(redrawnImage),
		redrawCounts.shapeCount, redrawCounts.segmentCount,
		redrawCounts.meshCount, redrawCounts.triangleCount,
		counts.shapeCount, counts.segmentCount,
		counts.meshCount, counts.triangleCount);
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

static int
sync_material_refresh_to_view(struct ged *gedp, QgView &view)
{
    /* Direct librt mutations bypass GED's material-change event.  Advance the
     * synchronization stamp explicitly before requesting a recolor sweep. */
    ged_draw_bump_material_revision(gedp);
    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_REFRESH_MATERIAL_COLORS, NULL);
    txn.view = view.viewContext();

    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    result.status = 1;
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    int changed = scene ? ged_draw_obol_scene_sync_transaction(gedp, &txn,
	    &result, scene) : 0;
    if (changed)
	view.need_update(QG_VIEW_REFRESH);
    ged_draw_transaction_result_free(&result);
    return changed;
}

static int
replace_global_color_table(struct ged *gedp, const char *table)
{
    if (!gedp || !gedp->dbip || !table)
	return 0;

    if (db5_update_attribute(DB5_GLOBAL_OBJECT_NAME, "regionid_colortable",
	    table, gedp->dbip) != 0)
	return 0;

    db_mater_free(gedp->dbip);
    std::vector<char> tableCopy(strlen(table) + 1);
    memcpy(tableCopy.data(), table, strlen(table) + 1);
    db5_import_color_table(gedp->dbip, tableCopy.data());
    return 1;
}

static int
exercise_m35_color_table_mutation(void)
{
    int64_t totalStart = bu_gettime();
    int64_t phaseStart = totalStart;
    auto report_phase = [&phaseStart](const char *name) {
	int64_t now = bu_gettime();
	fprintf(stderr, "m35_color_table_mutation:%s %.6f sec\n", name,
		(double)(now - phaseStart) / 1000000.0);
	phaseStart = now;
    };
    char src_db[MAXPATHLEN] = {0};
    if (!model_path("m35.g", src_db, sizeof(src_db))) {
	fprintf(stderr, "missing qtcad Obol m35 color-table source model: %s\n",
		src_db);
	return 0;
    }

    char tmp_db[MAXPATHLEN] = {0};
    FILE *tmp_fp = bu_temp_file(tmp_db, sizeof(tmp_db));
    if (!tmp_fp) {
	fprintf(stderr, "failed to allocate qtcad Obol m35 color-table temp database\n");
	return 0;
    }
    fclose(tmp_fp);

    if (!copy_file(src_db, tmp_db)) {
	fprintf(stderr, "failed to copy m35 color-table source database %s to %s\n",
		src_db, tmp_db);
	bu_file_delete(tmp_db);
	return 0;
    }

    struct ged *gedp = ged_open("db", tmp_db, 1);
    if (!gedp) {
	fprintf(stderr, "failed to open qtcad Obol m35 color-table copy: %s\n",
		tmp_db);
	bu_file_delete(tmp_db);
	return 0;
    }

    QgView view(NULL, QgViewType::SW);
    view.resize(220, 170);
	ged_view_active_ctx_set(gedp, view.viewContext());
	(void)ged_view_context_host_attach(gedp, view.viewContext());
    if (!ged_view_context_display_endpoint_set(
	    view.viewContext(), view.displayEndpoint(), 0)) {
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }

    BObolViewController *controller = view.obolViewController();
    if (!controller) {
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }
    controller->clearDatabaseSources();
    report_phase("setup");

    struct ged_draw_appearance_settings appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    appearance.draw_mode = GED_DRAW_MODE_WIRE;

    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "all.g");
    txn.view = view.viewContext();
    txn.appearance = &appearance;

    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int drawRet = ged_draw_apply_transaction(gedp, &txn, &result);
    report_phase("draw-transaction");
    int changed = qg_obol_sync_ged_draw_transaction(gedp, &txn, &result,
	    &view);
    if (drawRet < 0 || !changed) {
	fprintf(stderr, "m35 color-table draw/sync failed: draw=%d changed=%d errors=%s\n",
		drawRet, changed, bu_vls_cstr(&result.errors));
	ged_draw_transaction_result_free(&result);
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }
    ged_draw_transaction_result_free(&result);
    report_phase("initial-draw");

    SoBRLSceneController *source_controller =
	ged_draw_obol_scene_controller(gedp);

    BObolDatabaseSourceSummary canary_summary;
    SoBRLDatabaseSource *canary = find_source_by_path_suffix(source_controller,
	    "r850/s850", canary_summary);
    unsigned char expectedRgb[3] = {0, 0, 0};
    if (!canary || !source_material_matches_db_color(gedp, canary_summary,
	    expectedRgb)) {
	fprintf(stderr,
		"m35 color-table canary initial color mismatch: found=%d path=%s valid=%d expected=(%u %u %u) color=(%.9g %.9g %.9g)\n",
		canary ? 1 : 0,
		canary ? canary_summary.path.getString() : "",
		canary ? (int)canary_summary.materialColorValid : 0,
		expectedRgb[0], expectedRgb[1], expectedRgb[2],
		canary ? canary_summary.materialColor[0] : 0.0f,
		canary ? canary_summary.materialColor[1] : 0.0f,
		canary ? canary_summary.materialColor[2] : 0.0f);
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }

    /* Aggregate metadata describes the draw root, not every compact leaf.
     * A root fallback must not erase the individual full-path colors. */
    (void)canary->setDatabaseMetadataState(TRUE, 0, 0, 0, 0, TRUE,
	SbColor(1.0f, 1.0f, 1.0f), SbString("aggregate-test"));
    if (!all_source_materials_match_db_colors(gedp, source_controller)) {
	fprintf(stderr,
		"m35 aggregate metadata overrode compact occurrence colors\n");
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }

    const char *global_color_av[6] = {
	"color", "0", "15000", "20", "30", "40"
    };
    if (ged_exec_color(gedp, 6, global_color_av) != BRLCAD_OK) {
	fprintf(stderr, "m35 color-table global color mutation did not sync\n");
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }
    report_phase("global-color-refresh");

    if (!all_source_materials_match_db_colors(gedp, source_controller)) {
	fprintf(stderr, "m35 color-table cached sweep disagrees with full-path colors\n");
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }
    report_phase("global-color-reference-validation");

    canary = find_source_by_path_suffix(source_controller, "r850/s850",
	    canary_summary);
    if (!canary || !source_material_matches_db_color(gedp, canary_summary,
	    expectedRgb)) {
	fprintf(stderr,
		"m35 color-table canary global update mismatch: found=%d path=%s valid=%d expected=(%u %u %u) color=(%.9g %.9g %.9g)\n",
		canary ? 1 : 0,
		canary ? canary_summary.path.getString() : "",
		canary ? (int)canary_summary.materialColorValid : 0,
		expectedRgb[0], expectedRgb[1], expectedRgb[2],
		canary ? canary_summary.materialColor[0] : 0.0f,
		canary ? canary_summary.materialColor[1] : 0.0f,
		canary ? canary_summary.materialColor[2] : 0.0f);
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }

    const char *new_id_color_av[6] = {
	"color", "16001", "16001", "80", "90", "100"
    };
    if (ged_exec_color(gedp, 6, new_id_color_av) != BRLCAD_OK) {
	fprintf(stderr, "m35 color-table new region-id color mutation did not sync\n");
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }
    report_phase("new-id-color-refresh");

    const char *item_av[3] = {"item", "r850", "16001"};
    if (ged_exec_item(gedp, 3, item_av) != BRLCAD_OK) {
	fprintf(stderr, "m35 color-table region-id mutation did not sync\n");
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }
    report_phase("item-refresh");

    canary = find_source_by_path_suffix(source_controller, "r850/s850",
	    canary_summary);
    if (!canary || canary_summary.databaseRegionId != 16001 ||
	!source_material_matches_db_color(gedp, canary_summary, expectedRgb)) {
	fprintf(stderr,
		"m35 color-table canary region-id update mismatch: found=%d path=%s region=%d valid=%d expected=(%u %u %u) color=(%.9g %.9g %.9g)\n",
		canary ? 1 : 0,
		canary ? canary_summary.path.getString() : "",
		canary ? canary_summary.databaseRegionId : -1,
		canary ? (int)canary_summary.materialColorValid : 0,
		expectedRgb[0], expectedRgb[1], expectedRgb[2],
		canary ? canary_summary.materialColor[0] : 0.0f,
		canary ? canary_summary.materialColor[1] : 0.0f,
		canary ? canary_summary.materialColor[2] : 0.0f);
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }

    const char *direct_table =
	"{0 15000 60 70 80} {16001 16001 120 130 140} {17002 17002 200 210 220} ";
    if (!replace_global_color_table(gedp, direct_table) ||
	    !sync_material_refresh_to_view(gedp, view)) {
	fprintf(stderr, "m35 color-table direct _GLOBAL mutation did not sync\n");
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }
    report_phase("direct-global-refresh");

    canary = find_source_by_path_suffix(source_controller, "r850/s850",
	    canary_summary);
    if (!canary || !source_material_matches_db_color(gedp, canary_summary,
	    expectedRgb)) {
	fprintf(stderr,
		"m35 color-table canary direct _GLOBAL update mismatch: found=%d path=%s valid=%d expected=(%u %u %u) color=(%.9g %.9g %.9g)\n",
		canary ? 1 : 0,
		canary ? canary_summary.path.getString() : "",
		canary ? (int)canary_summary.materialColorValid : 0,
		expectedRgb[0], expectedRgb[1], expectedRgb[2],
		canary ? canary_summary.materialColor[0] : 0.0f,
		canary ? canary_summary.materialColor[1] : 0.0f,
		canary ? canary_summary.materialColor[2] : 0.0f);
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }

    if (db5_update_attribute("r850", "region_id", "17002", gedp->dbip) != 0 ||
	    !sync_material_refresh_to_view(gedp, view)) {
	fprintf(stderr, "m35 color-table direct region-id attribute mutation did not sync\n");
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }
    report_phase("direct-region-refresh");

    canary = find_source_by_path_suffix(source_controller, "r850/s850",
	    canary_summary);
    if (!canary || canary_summary.databaseRegionId != 17002 ||
	!source_material_matches_db_color(gedp, canary_summary, expectedRgb)) {
	fprintf(stderr,
		"m35 color-table canary direct region-id attr update mismatch: found=%d path=%s region=%d valid=%d expected=(%u %u %u) color=(%.9g %.9g %.9g)\n",
		canary ? 1 : 0,
		canary ? canary_summary.path.getString() : "",
		canary ? canary_summary.databaseRegionId : -1,
		canary ? (int)canary_summary.materialColorValid : 0,
		expectedRgb[0], expectedRgb[1], expectedRgb[2],
		canary ? canary_summary.materialColor[0] : 0.0f,
		canary ? canary_summary.materialColor[1] : 0.0f,
		canary ? canary_summary.materialColor[2] : 0.0f);
	ged_close(gedp);
	bu_file_delete(tmp_db);
	return 0;
    }

    print_timing({"m35_color_table_mutation", "m35.g", "all.g",
	    GED_DRAW_MODE_WIRE, SoBRLDatabaseSource::WIREFRAME,
	    0, 0, 0, 0, 0, 0}, "total", totalStart);

    ged_close(gedp);
    bu_file_delete(tmp_db);
    return 1;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const struct model_case cases[] = {
	{"pinewood_wire", "pinewood.g", "pinewood", GED_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 2, 41, 0, 0, 0, 0},
	{"pinewood_shaded", "pinewood.g", "pinewood", GED_DRAW_MODE_SHADED,
	    SoBRLDatabaseSource::SHADED, 0, 0, 21, 501, 0, 0},
	{"havoc_wire", "havoc.g", "havoc", GED_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 10, 100, 0, 0, 0, 0},
	{"m35_wire_interactions", "m35.g", "all.g", GED_DRAW_MODE_WIRE,
	    SoBRLDatabaseSource::WIREFRAME, 100, 1000, 0, 0, 1, 0},
	{"m35_wire_pick_all_stress", "m35.g", "all.g", GED_DRAW_MODE_WIRE,
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

    if (should_run_case(argc, argv, "m35_color_table_mutation")) {
	ran = 1;
	if (!exercise_m35_color_table_mutation())
	    FAIL("qtcad Obol m35 color-table mutation workflow should pass");
    }

    if (BU_STR_EQUAL(getenv("BOBOL_QTCAD_GENERIC_TWIN"), "1")) {
	const struct model_case genericTwinCase = {
	    "generic_twin_wire", "faa/Generic_Twin.g", "all", GED_DRAW_MODE_WIRE,
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
