/*                    A N N O T A T E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file annotate.cpp
 *
 * Command-level annotation creation and in-scene drawing regression tests.
 */

#include "common.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#include <bu.h>
#include <BObol/BDatabaseSource.h>
#include <BObol/BDisplayEndpoint.h>
#include <BObol/BExportAction.h>
#include <BObol/BSceneController.h>
#include <BObol/BViewController.h>
#include <BObol/BViewQuery.h>
#include <Inventor/nodes/SoCamera.h>
#include <Obol/cad/CadGeometry.h>
#include <ged/display_obol_private.h>
#include <ged/view.h>
#include <ged/scene.h>
#include <ged.h>
#include <rt/geom.h>
#include <rt/primitives/annot.h>
#include <wdb.h>

#include "view_test_util.h"
#include "../../ged_private.h"

#define ADIFF_THRESHOLD 0.99

static constexpr int ANNOTATE_IMAGE_SIZE = 1024;
static constexpr const char *ANNOTATE_TEXT_HEIGHT = "30";
static constexpr const char *ANNOTATE_DIMENSION_OFFSET = "65";
static constexpr const char *DIRECT_TEXT_HEIGHT = "60";
static constexpr const char *DIRECT_DIMENSION_OFFSET = "120";
static constexpr const char *PRIMITIVE_TEXT_HEIGHT = "40";
static constexpr fastf_t PRIMITIVE_DIMENSION_OFFSET_SCALE = 1.5;
static constexpr fastf_t PRIMITIVE_ANGULAR_OFFSET_SCALE = 3.0;
static constexpr unsigned int ANNOTATE_DRAIN_ATTEMPTS = 2000;
static constexpr unsigned int ANNOTATE_DRAIN_SLEEP_MILLISECONDS = 1;
static constexpr fastf_t OUTPUT_FILL_HALF_SIZE_SCALE = 0.15;
static constexpr const char *OUTPUT_PS = "annotate-output.ps";
static constexpr const char *OUTPUT_PLOT = "annotate-output.plot";
static constexpr const char *OUTPUT_PNG = "annotate-output.png";
static constexpr const char *OUTPUT_PNG_PLAIN = "annotate-output-plain.png";
static constexpr const char *OUTPUT_PNG_WITHOUT_FILL =
    "annotate-output-without-fill.png";

extern "C" void dm_refresh(struct ged *);
extern "C" int img_cmp(int, struct ged *, const char *, bool, bool, int, fastf_t,
	const char *, const char *);
extern "C" int unpack_apng(const char *, const char *, const char *, const char *);


static void
verify_annotation_coloring(struct ged *gedp)
{
    const char *name = "color-test-annotation";
    const char *group = "color-test-group";
    const char *create_argv[] = {
	"annotate", "text", "--no-draw", "--at", "0 0 0", name, "test", NULL
    };
    if (ged_exec_annotate(gedp, 7, create_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to create color test annotation: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));

    SbColor actual;
    const SbColor white(1.0f, 1.0f, 1.0f);
    if (!bobol_database_source_path_material_color(gedp->dbip, name, actual) ||
	!actual.equals(white, SMALL_FASTF))
	bu_exit(EXIT_FAILURE, "Annotation inherited the region color table\n");

    const char *group_argv[] = {"g", group, name, NULL};
    const char *color_argv[] = {"attr", "set", group, "color", "200/50/25", NULL};
    if (ged_exec(gedp, 3, group_argv) != BRLCAD_OK ||
	ged_exec(gedp, 5, color_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to create inherited annotation color fixture\n");
    const SbColor inherited(200.0f / 255.0f, 50.0f / 255.0f, 25.0f / 255.0f);
    if (!bobol_database_source_path_material_color(gedp->dbip,
	"color-test-group/color-test-annotation", actual) ||
	!actual.equals(inherited, SMALL_FASTF))
	bu_exit(EXIT_FAILURE, "Annotation lost its inherited color\n");

    struct wmember members;
    BU_LIST_INIT(&members.l);
    if (!mk_addmember("component", &members.l, NULL, WMOP_UNION) ||
	mk_comb(wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT), "color-test-region",
	    &members.l, 1, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0) != 0)
	bu_exit(EXIT_FAILURE, "Unable to create color-table test region\n");

    const struct mater *material = db_mater_head(gedp->dbip);
    while (material != MATER_NULL && (material->mt_low > 0 || material->mt_high < 0))
	material = material->mt_forw;
    if (material == MATER_NULL)
	bu_exit(EXIT_FAILURE, "Missing region color-table fixture\n");
    const SbColor region_color(material->mt_r / 255.0f,
	material->mt_g / 255.0f, material->mt_b / 255.0f);
    if (!bobol_database_source_path_material_color(gedp->dbip, "color-test-region", actual) ||
	!actual.equals(region_color, SMALL_FASTF))
	bu_exit(EXIT_FAILURE, "Region color-table behavior changed\n");

    const char *kill_argv[] = {"kill", group, name, "color-test-region", NULL};
    if (ged_exec(gedp, 4, kill_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to remove annotation color fixture\n");
}


static void
verify_screen_annotation(BObolViewController *controller)
{
    bool found_screen_plane = false;
    const SbVec2s size = controller->getViewportRegion().getViewportSizePixels();
    const SbViewVolume camera = controller->getCamera()->getViewVolume(
	static_cast<float>(size[0]) / size[1]);
    for (SoBRLDatabaseSource *source : controller->getRenderDatabaseSources()) {
	if (!source || !BU_STR_EQUAL(source->path.getValue().getString(),
		"component-screen-note"))
	    continue;
	BObolCompactOccurrence occurrence;
	if (!source->getCompactOccurrence(0, occurrence) || !occurrence.geometry ||
	    !occurrence.geometry->displayPlane)
	    bu_exit(EXIT_FAILURE, "Screen annotation lost its display-plane coordinates\n");
	const auto &plane = *occurrence.geometry->displayPlane;
	SbMatrix placement = occurrence.geometryTransform;
	placement.multRight(occurrence.localTransform);
	SbVec3f anchor;
	placement.multVecMatrix(plane.anchor, anchor);
	camera.projectToScreen(anchor, anchor);
	const SbVec3f center = Obol::cadPartGeometryBounds(*occurrence.geometry).getCenter();
	SoBRLExportAction export_action;
	export_action.applyViewport(*controller->getViewport());
	bool exported = false;
	bool exported_authored_width = false;
	for (int line_index = 0; line_index < export_action.getLineCount(); ++line_index) {
	    const auto &line = export_action.getLine(line_index);
	    if (!BU_STR_EQUAL(line.sourceName.getString(), "component-screen-note"))
		continue;
	    /* The screen leader is authored with --line-width 2 below. */
	    if (EQUAL(line.lineWidth, 2.0f))
		exported_authored_width = true;
	    const auto &wire = *occurrence.geometry->wire;
	    if (line.primitiveIndex < 0 ||
		static_cast<size_t>(line.primitiveIndex) >= wire.segmentPoints.size() / 2)
		bu_exit(EXIT_FAILURE, "Exported annotation lost its segment identity\n");
	    const SbVec3f &offset = wire.segmentPoints[2 * line.primitiveIndex];
	    SbVec3f projected;
	    camera.projectToScreen(line.a, projected);
	    const float expected_x = anchor[0] * size[0] + offset[0] * RT_ANNOT_DISPLAY_PIXELS_PER_MM;
	    const float expected_y = anchor[1] * size[1] + offset[1] * RT_ANNOT_DISPLAY_PIXELS_PER_MM;
	    const float pixel_tolerance = 0.05f;
	    if (fabs(projected[0] * size[0] - expected_x) > pixel_tolerance ||
		fabs(projected[1] * size[1] - expected_y) > pixel_tolerance)
		bu_exit(EXIT_FAILURE, "Exported screen annotation missed its pixel offset\n");
	    exported = true;
	}
	if (!exported)
	    bu_exit(EXIT_FAILURE, "Viewport export omitted the screen annotation\n");
	if (!exported_authored_width)
	    bu_exit(EXIT_FAILURE, "Viewport export lost the screen annotation's authored width\n");
	/* Independent pixel arithmetic: the rectangle must follow the label's
	 * screen offsets, rather than selecting only its model-space anchor. */
	const float x = 2.0f * (anchor[0] + center[0] *
	    RT_ANNOT_DISPLAY_PIXELS_PER_MM / size[0]) - 1.0f;
	const float y = 2.0f * (anchor[1] + center[1] *
	    RT_ANNOT_DISPLAY_PIXELS_PER_MM / size[1]) - 1.0f;
	const float half_width = 2.0f / size[0];
	const float half_height = 2.0f / size[1];
	std::vector<BObolViewPickRecord> records;
	if (source->queryCompactRectangle(SbMatrix::identity(), camera.getMatrix(),
		size, x - half_width, y - half_height, x + half_width,
		y + half_height, records) != 1)
	    bu_exit(EXIT_FAILURE, "Screen annotation rectangle selection missed its display bounds: "
		"visible=%d selectable=%d anchor=(%g,%g,%g) center=(%g,%g,%g) rectangle=(%g,%g) viewport=%d,%d\n",
		source->visible.getValue(), occurrence.summary.selectable,
		anchor[0], anchor[1], anchor[2], center[0], center[1], center[2], x, y, size[0], size[1]);
	/* Reusing the action without a viewport must restore viewless placement. */
	export_action.apply(source);
	SbVec3f model_point = plane.anchor + occurrence.geometry->wire->segmentPoints.front();
	placement.multVecMatrix(model_point, model_point);
	const float model_tolerance = 0.001f;
	if (export_action.getLineCount() == 0 ||
	    !export_action.getLine(0).a.equals(model_point, model_tolerance))
	    bu_exit(EXIT_FAILURE, "Viewless annotation export retained a previous camera\n");
	found_screen_plane = true;
    }
    if (!found_screen_plane)
	bu_exit(EXIT_FAILURE, "Missing retained screen annotation\n");
}


static BObolViewController *
settled_drawing_controller(struct ged *gedp)
{
    struct ged_view_context *view = ged_view_active_ctx(gedp);
    if (!draw_test_obol_progressive_drain(gedp, view,
	ANNOTATE_DRAIN_ATTEMPTS, ANNOTATE_DRAIN_SLEEP_MILLISECONDS))
	bu_exit(EXIT_FAILURE, "Annotation drawing did not settle\n");
    bobol_display_endpoint_t *endpoint = ged_view_context_obol_endpoint_get(view);
    BObolViewController *controller = endpoint ?
	static_cast<BObolViewController *>(bobol_display_endpoint_controller(endpoint)) : NULL;
    if (!controller || !controller->getRenderSceneRoot())
	bu_exit(EXIT_FAILURE, "Missing annotation scene\n");
    return controller;
}


static void
verify_component_materials(struct ged *gedp)
{
    BObolSceneController scene(settled_drawing_controller(gedp)->getRenderSceneRoot());
    /* The old image controls used the region table's brown for every leaf.
     * Explicit gray, white and black materials must survive alongside that
     * table fallback. */
    struct {
	const char *path;
	SbColor color;
	bool found = false;
    } samples[] = {
	{"/component/bed/r850/s850", SbColor(210.0f / 255.0f, 146.0f / 255.0f, 1.0f / 255.0f)},
	{"/component/bed/r851/s851", SbColor(180.0f / 255.0f, 180.0f / 255.0f, 180.0f / 255.0f)},
	{"/component/bed/r872/s872", SbColor(1.0f, 1.0f, 1.0f)},
	{"/component/cab/r828/s828", SbColor(10.0f / 255.0f, 10.0f / 255.0f, 10.0f / 255.0f)},
	{"/component/cab/r691/s691", SbColor(0.0f, 0.0f, 0.0f)}
    };
    for (int i = 0; i < scene.getRealizedShapeSummaryCount(); ++i) {
	BObolRealizedShapeSummary shape;
	if (!scene.getRealizedShapeSummary(i, shape) || !shape.valid || !shape.visible)
	    continue;
	for (auto &sample : samples) {
	    if (!BU_STR_EQUAL(shape.path.getString(), sample.path))
		continue;
	    const SbColor color = shape.colorOverride ? shape.color : shape.materialColor;
	    if (!shape.materialColorValid ||
		!color.equals(sample.color, SMALL_FASTF) || shape.segmentCount <= 0)
		bu_exit(EXIT_FAILURE, "Annotated model lost material or geometry for %s\n", sample.path);
	    sample.found = true;
	}
    }
    for (const auto &sample : samples)
	if (!sample.found)
	    bu_exit(EXIT_FAILURE, "Annotated model omitted %s\n", sample.path);
}


static void
verify_radial_dimensions(struct ged *gedp, const point_t center, fastf_t radius)
{
    BObolViewController *controller = settled_drawing_controller(gedp);
    SoBRLExportAction export_action;
    export_action.applyViewport(*controller->getViewport());
    const SbVec3f origin(center[X], center[Y], center[Z]);
    const SbVec3f radial(radius, 0.0f, 0.0f);
    const SbVec3f diameter(0.0f, radius, 0.0f);
    struct {
	const char *source;
	SbVec3f a;
	SbVec3f b;
	bool found = false;
    } spans[] = {
	{"sphere-radius", origin, origin + radial},
	{"sphere-diameter", origin - diameter, origin + diameter},
	{"sphere-angle", origin, origin + radial * PRIMITIVE_ANGULAR_OFFSET_SCALE},
	{"sphere-angle", origin, origin + diameter * PRIMITIVE_ANGULAR_OFFSET_SCALE}
    };
    /* Millimeter tolerance accommodates float storage at the fixture's model
     * coordinates, while detecting a second anchor translation. */
    const float model_tolerance = 0.01f;
    const float pixel_tolerance = 0.05f;
    const SbVec2s size = controller->getViewportRegion().getViewportSizePixels();
    const SbViewVolume camera = controller->getCamera()->getViewVolume(
	static_cast<float>(size[0]) / size[1]);
    SbBox3f bounds;
    for (int i = 0; i < export_action.getLineCount(); ++i) {
	const auto &line = export_action.getLine(i);
	for (auto &span : spans) {
	    if (BU_STR_EQUAL(line.sourceName.getString(), span.source) &&
		((line.a.equals(span.a, model_tolerance) && line.b.equals(span.b, model_tolerance)) ||
		 (line.a.equals(span.b, model_tolerance) && line.b.equals(span.a, model_tolerance))))
		span.found = true;
	}
	for (const SbVec3f &point : {line.a, line.b}) {
	    bounds.extendBy(point);
	    SbVec3f projected;
	    camera.projectToScreen(point, projected);
	    if (!std::isfinite(projected[0]) || !std::isfinite(projected[1]) ||
		projected[0] * size[0] < -pixel_tolerance ||
		projected[0] * size[0] > size[0] + pixel_tolerance ||
		projected[1] * size[1] < -pixel_tolerance ||
		projected[1] * size[1] > size[1] + pixel_tolerance)
		bu_exit(EXIT_FAILURE, "Autoview clipped radial dimension source %s\n",
		    line.sourceName.getString());
	}
    }
    for (const auto &span : spans)
	if (!span.found)
	    bu_exit(EXIT_FAILURE, "Radial dimension %s missed its authored span\n", span.source);
    /* Fit the delivered geometry's center. Per-object cube padding in the
     * old display list displaced this center even when every stroke fit. */
    SbVec3f projected_center;
    camera.projectToScreen(bounds.getCenter(), projected_center);
    if (fabs((projected_center[0] - 0.5f) * size[0]) > pixel_tolerance ||
	fabs((projected_center[1] - 0.5f) * size[1]) > pixel_tolerance)
	bu_exit(EXIT_FAILURE, "Autoview displaced the radial scene's geometry center\n");
}


static void
verify_visible_annotations(struct ged *gedp)
{
    BObolViewController *controller = settled_drawing_controller(gedp);
    BObolSceneController scene(controller->getRenderSceneRoot());

    /* Inspect delivered geometry: accepted draw intents alone did not catch
     * hiding one annotation cancelling other annotations' pending draws. */
    const struct {
	const char *name;
	int count;
	SbColor color;
    } expected[] = {
	{"component-dim", 0, SbColor(1.0f, 220.0f / 255.0f, 0.0f)},
	{"component-obb-dim", 3, SbColor(1.0f, 80.0f / 255.0f, 1.0f)}, /* Three axes. */
	{"component-note", 1, SbColor(0.0f, 1.0f, 1.0f)},
	{"component-screen-note", 1, SbColor(80.0f / 255.0f, 1.0f, 80.0f / 255.0f)}
    };
    const int shape_count = scene.getRealizedShapeSummaryCount();
    for (const auto &annotation : expected) {
	int count = 0;
	for (int i = 0; i < shape_count; ++i) {
	    BObolRealizedShapeSummary shape;
	    if (!scene.getRealizedShapeSummary(i, shape) || !shape.valid ||
		!shape.visible || !BU_STR_EQUAL(shape.ownerSourcePath.getString(), annotation.name))
		continue;
	    if (BU_STR_EQUAL(shape.geometryKind.getString(), "overview-aabb") ||
		BU_STR_EQUAL(shape.recordRole.getString(), "lod-overview") ||
		!shape.selectable || shape.segmentCount <= 0)
		bu_exit(EXIT_FAILURE, "Annotation %s retained preview state: kind=%s role=%s selectable=%d segments=%d\n",
		    annotation.name, shape.geometryKind.getString(), shape.recordRole.getString(), shape.selectable, shape.segmentCount);
	    const SbColor color = shape.colorOverride ? shape.color : shape.materialColor;
	    if ((!shape.colorOverride && !shape.materialColorValid) ||
		!color.equals(annotation.color, SMALL_FASTF))
		bu_exit(EXIT_FAILURE, "Annotation %s lost its color after hide/show or update: "
		    "got %.0f/%.0f/%.0f, expected %.0f/%.0f/%.0f\n", annotation.name,
		    color[0] * 255.0f, color[1] * 255.0f, color[2] * 255.0f,
		    annotation.color[0] * 255.0f, annotation.color[1] * 255.0f,
		    annotation.color[2] * 255.0f);
	    ++count;
	}
	if (count != annotation.count)
	    bu_exit(EXIT_FAILURE, "Annotation %s delivered %d shapes, expected %d\n",
		annotation.name, count, annotation.count);
    }

    verify_screen_annotation(controller);

}


static void
capture_image(struct ged *gedp, int frame)
{
    struct bu_vls name = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&name, "annotate%03d.png", frame);
    dm_refresh(gedp);
    const char *argv[] = {"screengrab", bu_vls_cstr(&name), NULL};
    if (ged_exec_screengrab(gedp, 2, argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to capture annotation frame %d\n", frame);
    bu_vls_free(&name);
}


static std::string
point_arg(const point_t point, fastf_t base2local)
{
    struct bu_vls value = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&value, "%.12g %.12g %.12g", point[X] * base2local,
	point[Y] * base2local, point[Z] * base2local);
    std::string result(bu_vls_cstr(&value));
    bu_vls_free(&value);
    return result;
}


static std::string
annotation_text(struct ged *gedp, const char *name)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (dp == RT_DIR_NULL ||
	rt_db_get_internal(&intern, dp, gedp->dbip, NULL) != ID_ANNOT)
	bu_exit(EXIT_FAILURE, "Unable to read annotation text from %s\n", name);
    struct rt_annot_internal *annotation =
	static_cast<struct rt_annot_internal *>(intern.idb_ptr);
    std::string label;
    for (size_t i = 0; i < annotation->ant.count; ++i) {
	uint32_t magic = *static_cast<uint32_t *>(annotation->ant.segments[i]);
	if (magic == ANN_TSEG_MAGIC) {
	    struct txt_seg *text = static_cast<struct txt_seg *>(annotation->ant.segments[i]);
	    label = bu_vls_cstr(&text->label);
	    break;
	}
    }
    rt_db_free_internal(&intern);
    return label;
}


static double
annotation_value(struct ged *gedp, const char *name)
{
    const std::string label = annotation_text(gedp, name);
    const size_t separator = label.find(':');
    const char *number = label.c_str() +
	(separator == std::string::npos ? 0 : separator + 1);
    char *end = NULL;
    const double value = strtod(number, &end);
    if (end == number)
	bu_exit(EXIT_FAILURE, "Unable to parse annotation value '%s'\n", label.c_str());
    return value;
}


static void
annotation_anchor(point_t anchor, struct ged *gedp, const char *name)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (dp == RT_DIR_NULL ||
	rt_db_get_internal(&intern, dp, gedp->dbip, NULL) != ID_ANNOT)
	bu_exit(EXIT_FAILURE, "Unable to read annotation anchor from %s\n", name);
    struct rt_annot_internal *annotation =
	static_cast<struct rt_annot_internal *>(intern.idb_ptr);
    VMOVE(anchor, annotation->V);
    rt_db_free_internal(&intern);
}


static bool
annotation_is_screen_space(struct ged *gedp, const char *name)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (dp == RT_DIR_NULL ||
	rt_db_get_internal(&intern, dp, gedp->dbip, NULL) != ID_ANNOT)
	bu_exit(EXIT_FAILURE, "Unable to read annotation placement from %s\n", name);
    struct rt_annot_internal *annotation =
	static_cast<struct rt_annot_internal *>(intern.idb_ptr);
    const bool screen_space = !(annotation->flags & RT_ANNOT_MODEL_SPACE);
    rt_db_free_internal(&intern);
    return screen_space;
}


static void
verify_help(struct ged *gedp, int argc, const char **argv, const char *expected)
{
    const int ret = ged_exec_annotate(gedp, argc, argv);
    if (ret != GED_HELP || !strstr(bu_vls_cstr(gedp->ged_result_str), expected))
	bu_exit(EXIT_FAILURE, "annotate help failed for '%s': %s\n", argv[1],
	    bu_vls_cstr(gedp->ged_result_str));
}


static void
verify_help_output(struct ged *gedp)
{
    const char *root_argv[] = {"annotate", "--help", NULL};
    verify_help(gedp, 2, root_argv, "Available subcommands:");
    const char *targeted_leader_argv[] = {"annotate", "help", "leader", NULL};
    verify_help(gedp, 3, targeted_leader_argv,
	"Create a text callout with a leader");
    const char *leader_argv[] = {"annotate", "leader", "--help", NULL};
    verify_help(gedp, 3, leader_argv, "--dpi");
    const char *dimension_argv[] = {
	"annotate", "dimension", "linear", "--help", NULL
    };
    verify_help(gedp, 4, dimension_argv,
	"Measures the distance between --from and --to.");
    const char *dimension_root_argv[] = {"annotate", "dimension", "--help", NULL};
    verify_help(gedp, 3, dimension_root_argv, "Available subcommands:");
    const char *update_argv[] = {"annotate", "update", "--help", NULL};
    verify_help(gedp, 3, update_argv, "--view-only");
}


static std::string
drawing_intents(struct ged *gedp)
{
    struct bu_vls paths = BU_VLS_INIT_ZERO;
    (void)ged_scene_paths_append(gedp, ged_view_active_ctx(gedp),
	GED_SCENE_DRAW_DEFAULT, GED_SCENE_PATHS_DRAW_INTENTS, &paths);
    const std::string result(bu_vls_cstr(&paths));
    bu_vls_free(&paths);
    return result;
}


static void
verify_geometry_update(struct ged *gedp)
{
    const char *zap_argv[] = {"zap", NULL};
    if (ged_exec_zap(gedp, 1, zap_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to isolate annotation update scene\n");
    const std::string hidden_intents = drawing_intents(gedp);
    const char *source_name = "annotate-update-source";
    const char *dimension_name = "annotate-update-dim";
    const char *leader_name = "annotate-update-leader";
    point_t center = VINIT_ZERO;
    struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
    if (mk_sph(wdbp, source_name, center, 10.0))
	bu_exit(EXIT_FAILURE, "Unable to create autodim update source\n");

    /* Input discovery must work from draw intent without a render pass. */
    const char *draw_source_argv[] = {"draw", source_name, NULL};
    const char *inferred_argv[] = {"annotate", "autodim", "--no-draw",
	"--axes", "x", "--precision", "1", "annotate-inferred-dim", NULL};
    if (ged_exec_draw(gedp, 2, draw_source_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 8, inferred_argv) != BRLCAD_OK ||
	!NEAR_EQUAL(annotation_value(gedp, "annotate-inferred-dim.x"), 20.0 * gedp->dbip->dbi_base2local, 0.11))
	bu_exit(EXIT_FAILURE, "Autodim could not discover its input drawing intent: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    const char *erase_source_argv[] = {"erase", source_name, NULL};
    const char *kill_inferred_argv[] = {"kill", "annotate-inferred-dim",
	"annotate-inferred-dim.x", NULL};
    if (ged_exec_erase(gedp, 2, erase_source_argv) != BRLCAD_OK ||
	ged_exec(gedp, 3, kill_inferred_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to clear inferred annotation fixture\n");

    const char *create_argv[] = {
	"annotate", "autodim", "--no-draw", "--axes", "x,y", "--precision", "1",
	dimension_name, source_name, NULL
    };
    if (ged_exec_annotate(gedp, 9, create_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Initial autodim update value is incorrect\n");
    const double initial_value = annotation_value(gedp, "annotate-update-dim.x");

    const char *leader_argv[] = {
	"annotate", "leader", "--no-draw", "--for", source_name, "--at",
	"60 0 0", leader_name, "UPDATED LEADER", NULL
    };
    if (ged_exec_annotate(gedp, 10, leader_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to create leader update source: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    point_t initial_target;
    annotation_anchor(initial_target, gedp, leader_name);

    const char *colored_annotations[] = {dimension_name, leader_name};
    /* Earlier annotation writers stored the standard color under its rgb
     * alias. Exercise those records as well as normal attr-command edits. */
    for (const char *name : colored_annotations)
	if (db5_update_attribute(name, "rgb", "80/120/160", gedp->dbip) < 0)
	    bu_exit(EXIT_FAILURE, "Unable to store legacy color for %s\n", name);
    const auto verify_colors = [&](const SbColor &expected) {
	for (const char *name : colored_annotations) {
	    SbColor actual;
	    if (!bobol_database_source_path_material_color(gedp->dbip, name, actual) ||
		!actual.equals(expected, SMALL_FASTF))
		bu_exit(EXIT_FAILURE, "Annotation update lost stored color for %s\n", name);
	}
    };

    struct directory *source_dp = db_lookup(gedp->dbip, source_name, LOOKUP_QUIET);
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (source_dp == RT_DIR_NULL ||
	rt_db_get_internal(&intern, source_dp, gedp->dbip, NULL) != ID_ELL)
	bu_exit(EXIT_FAILURE, "Unable to read autodim update source\n");
    struct rt_ell_internal *sphere = static_cast<struct rt_ell_internal *>(intern.idb_ptr);
    VSET(sphere->a, 20.0, 0.0, 0.0);
    VSET(sphere->b, 0.0, 20.0, 0.0);
    VSET(sphere->c, 0.0, 0.0, 20.0);
    if (rt_db_put_internal(source_dp, gedp->dbip, &intern) < 0)
	bu_exit(EXIT_FAILURE, "Unable to resize autodim update source\n");

    const char *update_leader_argv[] = {"annotate", "update", leader_name, NULL};
    if (ged_exec_annotate(gedp, 3, update_leader_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Leader did not update after a geometry change: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    point_t updated_target;
    annotation_anchor(updated_target, gedp, leader_name);
    if (DIST_PNT_PNT(updated_target, initial_target) <= SMALL_FASTF)
	bu_exit(EXIT_FAILURE, "Leader target did not track resized geometry\n");

    const char *view_update_argv[] = {
	"annotate", "update", "--view-only", dimension_name, NULL
    };
    if (ged_exec_annotate(gedp, 4, view_update_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Autodim view-only update failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    const double cached_value = annotation_value(gedp, "annotate-update-dim.x");
    if (!NEAR_EQUAL(cached_value, initial_value, 0.01))
	bu_exit(EXIT_FAILURE, "Autodim view-only update changed a cached measurement\n");

    const char *update_argv[] = {"annotate", "update", dimension_name, NULL};
    if (ged_exec_annotate(gedp, 3, update_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Autodim did not update after a geometry change: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    const double updated_value = annotation_value(gedp, "annotate-update-dim.x");
    if (!NEAR_EQUAL(updated_value, initial_value * 2.0, 0.11))
	bu_exit(EXIT_FAILURE, "Autodim value did not track resized geometry\n");

    if (drawing_intents(gedp) != hidden_intents)
	bu_exit(EXIT_FAILURE, "Updating hidden annotations changed the drawing intents\n");
    verify_colors(SbColor(80.0f / 255.0f, 120.0f / 255.0f, 160.0f / 255.0f));
    for (const char *name : colored_annotations) {
	const char *color_argv[] = {"attr", "set", name,
	    db5_standard_attribute(ATTR_COLOR), "120/30/90", NULL};
	if (ged_exec(gedp, 5, color_argv) != BRLCAD_OK)
	    bu_exit(EXIT_FAILURE, "Unable to edit annotation color for %s\n", name);
    }

    const char *parent_argv[] = {"g", "annotate-update-parent", dimension_name,
	leader_name, NULL};
    const char *draw_parent_argv[] = {"draw", "-C", "9/80/150",
	"annotate-update-parent", NULL};
    if (ged_exec(gedp, 4, parent_argv) != BRLCAD_OK ||
	ged_exec_draw(gedp, 4, draw_parent_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to draw nested annotation update fixture\n");
    const std::string nested_intents = drawing_intents(gedp);
    if (ged_exec_annotate(gedp, 4, view_update_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 3, update_leader_argv) != BRLCAD_OK ||
	drawing_intents(gedp) != nested_intents)
	bu_exit(EXIT_FAILURE, "Updating nested annotations changed their drawing roots\n");
    verify_colors(SbColor(120.0f / 255.0f, 30.0f / 255.0f, 90.0f / 255.0f));
    dm_refresh(gedp);
    const char *erase_parent_argv[] = {"erase", "annotate-update-parent", NULL};
    if (ged_exec_erase(gedp, 2, erase_parent_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 3, update_argv) != BRLCAD_OK ||
	drawing_intents(gedp) != hidden_intents)
	bu_exit(EXIT_FAILURE, "Updating erased annotations made them visible\n");

    struct directory *missing_member = db_lookup(gedp->dbip, "annotate-update-dim.y",
	LOOKUP_QUIET);
    if (missing_member == RT_DIR_NULL || db_delete(gedp->dbip, missing_member) ||
	db_dirdelete(gedp->dbip, missing_member))
	bu_exit(EXIT_FAILURE, "Unable to prepare missing-autodim-member update test\n");
    if (ged_exec_annotate(gedp, 3, update_argv) == BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Autodim update unexpectedly accepted a missing member\n");
    const double member_failure_value = annotation_value(gedp, "annotate-update-dim.x");
    if (!NEAR_EQUAL(member_failure_value, updated_value, 0.01))
	bu_exit(EXIT_FAILURE,
	    "Failed autodim member update replaced an existing annotation\n");

    struct directory *dimension_dp = db_lookup(gedp->dbip, dimension_name, LOOKUP_QUIET);
    struct bu_attribute_value_set avs = BU_AVS_INIT_ZERO;
    if (dimension_dp == RT_DIR_NULL ||
	db5_get_attributes(gedp->dbip, &avs, dimension_dp) ||
	bu_avs_add(&avs, "annotate:sources", "annotate-missing-source") < 0 ||
	db5_update_attributes(dimension_dp, &avs, gedp->dbip)) {
	bu_avs_free(&avs);
	bu_exit(EXIT_FAILURE, "Unable to prepare failed-autodim update test\n");
    }
    bu_avs_free(&avs);
    if (ged_exec_annotate(gedp, 3, update_argv) == BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Autodim update unexpectedly accepted a missing source\n");
    const double retained_value = annotation_value(gedp, "annotate-update-dim.x");
    if (!NEAR_EQUAL(retained_value, updated_value, 0.01))
	bu_exit(EXIT_FAILURE,
	    "Failed autodim update replaced the existing annotation\n");

    struct directory *leader_dp = db_lookup(gedp->dbip, "annotate-update-leader",
	LOOKUP_QUIET);
    struct bu_attribute_value_set leader_avs = BU_AVS_INIT_ZERO;
    if (leader_dp == RT_DIR_NULL ||
	db5_get_attributes(gedp->dbip, &leader_avs, leader_dp) ||
	bu_avs_add(&leader_avs, "annotate:sources", "annotate-missing-source") < 0 ||
	db5_update_attributes(leader_dp, &leader_avs, gedp->dbip)) {
	bu_avs_free(&leader_avs);
	bu_exit(EXIT_FAILURE, "Unable to prepare failed-leader update test\n");
    }
    bu_avs_free(&leader_avs);
    if (ged_exec_annotate(gedp, 3, update_leader_argv) == BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Leader update unexpectedly accepted a missing source\n");
    point_t retained_target;
    annotation_anchor(retained_target, gedp, "annotate-update-leader");
    if (DIST_PNT_PNT(retained_target, updated_target) > SMALL_FASTF)
	bu_exit(EXIT_FAILURE,
	    "Failed leader update replaced the existing annotation\n");
}


static bool
make_output_fill_annotation(struct ged *gedp, const point_t anchor,
	fastf_t halfSize, const char *name, bool styledLine)
{
    if (!gedp || !gedp->dbip || !name)
	return false;
    struct rt_annot_internal annot = {};
    annot.magic = RT_ANNOT_INTERNAL_MAGIC;
    annot.flags = RT_ANNOT_MODEL_SPACE;
    VMOVE(annot.V, anchor);
    VSET(annot.u_vec, 1.0, 0.0, 0.0);
    VSET(annot.v_vec, 0.0, 0.0, 1.0);
    point2d_t vertices[] = {
	{-halfSize, -halfSize}, {halfSize, -halfSize},
	{halfSize, halfSize}, {-halfSize, halfSize}
    };
    int indices[] = {0, 1, 2, 3};
    int ends[] = {4};
    struct fill_seg fill = {};
    fill.magic = ANN_FSEG_MAGIC;
    fill.loop_count = 1;
    fill.point_count = 4;
    fill.loop_ends = ends;
    fill.points = indices;
    fill.legacy_start = 0;
    fill.legacy_count = 4;
    struct line_seg outlines[4] = {};
    struct line_seg line = {};
    line.magic = CURVE_LSEG_MAGIC;
    line.start = 0;
    line.end = 2;
    void *segments[6] = {};
    int reverse[6] = {};
    for (int i = 0; i < 4; ++i) {
	outlines[i].magic = CURVE_LSEG_MAGIC;
	outlines[i].start = i;
	outlines[i].end = (i + 1) % 4;
	segments[i] = &outlines[i];
    }
    segments[4] = &line;
    segments[5] = &fill;
    struct rt_annot_seg_style styles[6] = {};
    if (styledLine) {
	styles[4].flags = RT_ANNOT_STYLE_WIDTH | RT_ANNOT_STYLE_COLOR;
	styles[4].line_pattern = RT_ANNOT_LINE_DASHED;
	styles[4].line_width = 3.0;
	styles[4].color[1] = 255;
	styles[4].color[2] = 255;
	styles[4].color[3] = 255;
    }
    styles[5].flags = RT_ANNOT_STYLE_COLOR;
    styles[5].color[0] = 255;
    styles[5].color[3] = 255;
    annot.vert_count = 4;
    annot.verts = vertices;
    annot.ant.count = 6;
    annot.ant.segments = segments;
    annot.ant.reverse = reverse;
    annot.styles = styles;
    return mk_annot(wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT),
	name, &annot) == 0;
}


static std::string
read_output_file(const char *path)
{
    std::ifstream stream(path, std::ios::binary);
    return stream ? std::string(std::istreambuf_iterator<char>(stream),
	std::istreambuf_iterator<char>()) : std::string();
}


static void
verify_output_consumers(struct ged *gedp, const point_t anchor,
	fastf_t fillHalfSize)
{
    const char *styledName = "annotate-output-styled";
    const char *plainName = "annotate-output-plain";
    const char *zap_argv[] = {"zap", NULL};
    if (ged_exec_zap(gedp, 1, zap_argv) != BRLCAD_OK ||
	!make_output_fill_annotation(gedp, anchor, fillHalfSize, styledName, true) ||
	!make_output_fill_annotation(gedp, anchor, fillHalfSize, plainName, false))
	bu_exit(EXIT_FAILURE, "Unable to create annotation output fixture\n");
    const char *draw_argv[] = {"draw", styledName, NULL};
    const char *autoview_argv[] = {"autoview", NULL};
    const char *ae_argv[] = {"ae", "45", "35", NULL};
    if (ged_exec_draw(gedp, 2, draw_argv) != BRLCAD_OK ||
	ged_exec_autoview(gedp, 1, autoview_argv) != BRLCAD_OK ||
	ged_exec_ae(gedp, 3, ae_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to draw annotation output fixture\n");
    (void)settled_drawing_controller(gedp);

    const char *ps_argv[] = {"ps", "-l", "3", OUTPUT_PS, NULL};
    const char *plot_argv[] = {"plot", OUTPUT_PLOT, NULL};
    const char *png_argv[] = {"png", "-c", "12/34/56", "-s", "256",
	OUTPUT_PNG, NULL};
    if (ged_exec(gedp, 4, ps_argv) != BRLCAD_OK ||
	ged_exec(gedp, 2, plot_argv) != BRLCAD_OK ||
	ged_exec(gedp, 6, png_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Annotation output command failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));

    const std::string postscript = read_output_file(OUTPUT_PS);
    const std::string plot = read_output_file(OUTPUT_PLOT);
    const std::string png = read_output_file(OUTPUT_PNG);
    if (postscript.find("9 setlinewidth") == std::string::npos ||
	postscript.find("[24 24] 0 setdash") == std::string::npos ||
	postscript.find("closepath fill") == std::string::npos ||
	postscript.find("0.000000 1.000000 1.000000 setrgbcolor") ==
	    std::string::npos)
	bu_exit(EXIT_FAILURE,
	    "PostScript output lost annotation width, pattern, color, or fill\n");
    if (plot.find("shortdashed\n") == std::string::npos)
	bu_exit(EXIT_FAILURE, "Plot output lost the declared dashed mapping\n");
    if (png.size() < 8 || png.compare(1, 3, "PNG") != 0)
	bu_exit(EXIT_FAILURE, "PNG output fixture is invalid\n");

    const char *erase_argv[] = {"erase", styledName, NULL};
    const char *draw_plain_argv[] = {"draw", plainName, NULL};
    const char *png_plain_argv[] = {
	"png", "-c", "12/34/56", "-s", "256", OUTPUT_PNG_PLAIN, NULL
    };
    const char *png_without_fill_argv[] = {
	"png", "-c", "12/34/56", "-s", "256", OUTPUT_PNG_WITHOUT_FILL, NULL
    };
    if (ged_exec_erase(gedp, 2, erase_argv) != BRLCAD_OK ||
	ged_exec_draw(gedp, 2, draw_plain_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to replace the styled output fixture\n");
    (void)settled_drawing_controller(gedp);
    if (ged_exec(gedp, 6, png_plain_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to write the plain annotation PNG\n");
    const std::string plainPng = read_output_file(OUTPUT_PNG_PLAIN);
    if (plainPng == png)
	bu_exit(EXIT_FAILURE, "PNG output did not consume annotation line style\n");

    const char *erase_plain_argv[] = {"erase", plainName, NULL};
    if (ged_exec_erase(gedp, 2, erase_plain_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to erase the plain output fixture\n");
    (void)settled_drawing_controller(gedp);
    if (ged_exec(gedp, 6, png_without_fill_argv) != BRLCAD_OK ||
	read_output_file(OUTPUT_PNG_WITHOUT_FILL) == plainPng)
	bu_exit(EXIT_FAILURE, "PNG output did not consume annotation fill geometry\n");

    bu_file_delete(OUTPUT_PS);
    bu_file_delete(OUTPUT_PLOT);
    bu_file_delete(OUTPUT_PNG);
    bu_file_delete(OUTPUT_PNG_PLAIN);
    bu_file_delete(OUTPUT_PNG_WITHOUT_FILL);
}


int
main(int argc, const char **argv)
{
    int generate = 0;
    int keep_images = 0;
    int continue_on_failure = 0;
    int commands_only = 0;
    int ret = BRLCAD_OK;
    struct bu_opt_desc options[5];
    BU_OPT(options[0], "G", "generate", "", NULL, &generate,
	"Generate PNG frames without comparing controls");
    BU_OPT(options[1], "k", "keep", "", NULL, &keep_images,
	"Keep generated PNG frames");
    BU_OPT(options[2], "c", "continue", "", NULL, &continue_on_failure,
	"Continue after an image mismatch");
    BU_OPT(options[3], "n", "commands-only", "", NULL, &commands_only,
	"Run command assertions without image comparison");
    BU_OPT_NULL(options[4]);

    bu_setprogname(argv[0]);
    argc--; argv++;
    int remaining = bu_opt_parse(NULL, argc, argv, options);
    if (remaining != 2)
	bu_exit(EXIT_FAILURE,
	    "Usage: ged_test_annotate [-G] [-k] [-c] [-n] control-directory m35.g\n");
    const char *control_dir = argv[0];
    const char *m35_path = argv[1];
    if (!bu_file_directory(control_dir) || !bu_file_exists(m35_path, NULL))
	bu_exit(EXIT_FAILURE, "Annotation drawing test inputs are unavailable\n");

    const char *working_db = "m35_annotate_tmp.g";
    std::ifstream original(m35_path, std::ios::binary);
    std::ofstream temporary(working_db, std::ios::binary);
    if (!original || !temporary)
	bu_exit(EXIT_FAILURE, "Unable to prepare the annotation test database\n");
    temporary << original.rdbuf();
    if (!temporary)
	bu_exit(EXIT_FAILURE, "Unable to copy the annotation test database\n");
    original.close();
    temporary.close();

    char cache_dir[MAXPATHLEN] = {0};
    char runtime_cache[MAXPATHLEN] = {0};
    bu_dir(cache_dir, MAXPATHLEN, BU_DIR_CURR, "ged_annotate_test_cache", NULL);
    bu_mkdir(cache_dir);
    bu_dir(runtime_cache, MAXPATHLEN, BU_DIR_CURR, "ged_annotate_test_cache",
	"cache", NULL);
    bu_mkdir(runtime_cache);
    bu_setenv("BU_DIR_CACHE", runtime_cache, 1);
    if (!generate && !commands_only && unpack_apng(control_dir, "annotate.apng", cache_dir,
	"annotate"))
	bu_exit(EXIT_FAILURE, "Unable to unpack annotation controls\n");

    struct ged *gedp = ged_open("db", working_db, 1);
    if (!gedp)
	bu_exit(EXIT_FAILURE, "Unable to open annotation test database\n");
    verify_help_output(gedp);
    if (draw_test_obol_view_init(gedp, ged_view_active_ctx(gedp),
	ANNOTATE_IMAGE_SIZE, ANNOTATE_IMAGE_SIZE) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to initialize annotation drawing endpoint\n");
    verify_annotation_coloring(gedp);

    const char *draw_argv[] = {"draw", "component", NULL};
    const char *autoview_argv[] = {"autoview", NULL};
    const char *ae_argv[] = {"ae", "45", "35", NULL};
    if (ged_exec_draw(gedp, 2, draw_argv) != BRLCAD_OK ||
	ged_exec_autoview(gedp, 1, autoview_argv) != BRLCAD_OK ||
	ged_exec_ae(gedp, 3, ae_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to prepare annotation test scene\n");

    const char *autodim_argv[] = {
	"annotate", "autodim", "--color", "255/220/0", "--precision", "1",
	"--text-height", ANNOTATE_TEXT_HEIGHT, "--line-width", "2", "--offset",
	ANNOTATE_DIMENSION_OFFSET,
	"component-dim", "component", NULL
    };
    if (ged_exec_annotate(gedp, 14, autodim_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "autodim failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    (void)ged_exec_autoview(gedp, 1, autoview_argv);
    (void)ged_exec_ae(gedp, 3, ae_argv);
    if (generate && !commands_only)
	capture_image(gedp, 1);
    else if (!commands_only)
	ret += img_cmp(1, gedp, cache_dir, false, !keep_images, continue_on_failure,
	    ADIFF_THRESHOLD, "annotate_clear", "annotate");
    verify_component_materials(gedp);

    const char *hide_aabb_argv[] = {"annotate", "hide", "component-dim", NULL};
    const char *obb_argv[] = {
	"annotate", "autodim", "--bounds", "obb", "--color", "255/80/255",
	"--precision", "1", "--text-height", ANNOTATE_TEXT_HEIGHT,
	"--line-width", "2", "--offset", ANNOTATE_DIMENSION_OFFSET,
	"component-obb-dim", "component", NULL
    };
    if (ged_exec_annotate(gedp, 3, hide_aabb_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 16, obb_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "oriented autodim failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    (void)ged_exec_autoview(gedp, 1, autoview_argv);
    (void)ged_exec_ae(gedp, 3, ae_argv);
    if (generate && !commands_only)
	capture_image(gedp, 2);
    else if (!commands_only)
	ret += img_cmp(2, gedp, cache_dir, false, !keep_images, continue_on_failure,
	    ADIFF_THRESHOLD, "annotate_clear", "annotate");

    const char *updated_ae_argv[] = {"ae", "125", "25", NULL};
    const char *update_obb_argv[] = {
	"annotate", "update", "--view-only", "component-obb-dim", NULL
    };
    if (ged_exec_ae(gedp, 3, updated_ae_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 4, update_obb_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to update autodim for the current view: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    (void)ged_exec_autoview(gedp, 1, autoview_argv);
    if (generate && !commands_only)
	capture_image(gedp, 3);
    else if (!commands_only)
	ret += img_cmp(3, gedp, cache_dir, false, !keep_images, continue_on_failure,
	    ADIFF_THRESHOLD, "annotate_clear", "annotate");

    const char *hide_argv[] = {"annotate", "hide", "component", NULL};
    const char *show_argv[] = {"annotate", "show", "component", NULL};
    if (ged_exec_annotate(gedp, 3, hide_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 3, show_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to toggle annotations for component: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));

    point_t bmin, bmax, label;
    if (rt_obj_bounds(gedp->ged_result_str, gedp->dbip, 1, &draw_argv[1], 1,
	bmin, bmax) & BRLCAD_ERROR)
	bu_exit(EXIT_FAILURE, "Unable to find component bounds\n");
    VMOVE(label, bmax);
    label[X] += (bmax[X] - bmin[X]) * 0.15;
    label[Z] += (bmax[Z] - bmin[Z]) * 0.1;
    std::string label_arg = point_arg(label, gedp->dbip->dbi_base2local);
    const char *leader_argv[] = {
	"annotate", "leader", "--for", "component", "--at",
	label_arg.c_str(), "--color", "0/255/255", "--text-height",
	ANNOTATE_TEXT_HEIGHT, "--line-width", "2", "--bold", "--italic", "component-note",
	"M35 component", NULL
    };
    if (ged_exec_annotate(gedp, 16, leader_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "leader annotation failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    point_t leader_target, bounds_center;
    annotation_anchor(leader_target, gedp, "component-note");
    VADD2SCALE(bounds_center, bmin, bmax, 0.5);
    if (DIST_PNT_PNT(leader_target, bounds_center) <= SMALL_FASTF)
	bu_exit(EXIT_FAILURE, "default leader target remained at the bounds center\n");
    const char *update_leader_argv[] = {
	"annotate", "update", "--view-only", "component-note", NULL
    };
    if (ged_exec_annotate(gedp, 4, update_leader_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "leader view-only update failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    if (ged_exec_annotate(gedp, 3, hide_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 3, show_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 3, hide_aabb_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to toggle associated leader: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    (void)ged_exec_autoview(gedp, 1, autoview_argv);
    point_t target_view, screen_label_view, screen_label;
    MAT4X3PNT(target_view, DRAW_TEST_BV_CONST(ged_view_active_ctx(gedp))->model2view, leader_target);
    VSET(screen_label_view, target_view[X] + 0.45, target_view[Y] + 0.35,
	target_view[Z]);
    MAT4X3PNT(screen_label, DRAW_TEST_BV_CONST(ged_view_active_ctx(gedp))->view2model, screen_label_view);
    const std::string screen_label_arg = point_arg(screen_label,
	gedp->dbip->dbi_base2local);
    const char *screen_leader_argv[] = {
	"annotate", "leader", "--screen-space", "--for", "component", "--at",
	screen_label_arg.c_str(), "--color", "80/255/80", "--dpi", "120",
	"--text-height", "0.2",
	"--line-width", "2",
	"component-screen-note", "VIEW FACING", NULL
    };
    if (ged_exec_annotate(gedp, 17, screen_leader_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "screen-space leader creation failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    if (!annotation_is_screen_space(gedp, "component-screen-note"))
	bu_exit(EXIT_FAILURE, "screen-space leader was stored in model space\n");
    verify_visible_annotations(gedp);
    if (generate && !commands_only)
	capture_image(gedp, 4);
    else if (!commands_only)
	ret += img_cmp(4, gedp, cache_dir, false, !keep_images, continue_on_failure,
	    ADIFF_THRESHOLD, "annotate_clear", "annotate");

    const char *screen_ae_argv[] = {"ae", "-70", "40", NULL};
    const char *screen_update_argv[] = {
	"annotate", "update", "--view-only", "component-obb-dim",
	"component-screen-note", NULL
    };
    if (ged_exec_ae(gedp, 3, screen_ae_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 5, screen_update_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to update annotations for a second view: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    if (!annotation_is_screen_space(gedp, "component-screen-note"))
	bu_exit(EXIT_FAILURE, "screen-space update changed annotation coordinates\n");
    verify_visible_annotations(gedp);
    if (generate && !commands_only)
	capture_image(gedp, 5);
    else if (!commands_only)
	ret += img_cmp(5, gedp, cache_dir, false, !keep_images, continue_on_failure,
	    ADIFF_THRESHOLD, "annotate_clear", "annotate");

    const char *hide_existing_argv[] = {
	"annotate", "hide", "component-obb-dim", "component-note",
	"component-screen-note", NULL
    };
    point_t linear_from, linear_to, ordinate_to, text_at;
    VSET(linear_from, bmin[X], bmax[Y], bmax[Z]);
    VSET(linear_to, bmax[X], bmax[Y], bmax[Z]);
    VSET(ordinate_to, bmin[X], bmin[Y], bmax[Z]);
    VSET(text_at, bmin[X], bmax[Y], bmax[Z] + (bmax[Z] - bmin[Z]) * 0.15);
    const std::string linear_from_arg = point_arg(linear_from, gedp->dbip->dbi_base2local);
    const std::string linear_to_arg = point_arg(linear_to, gedp->dbip->dbi_base2local);
    const std::string ordinate_to_arg = point_arg(ordinate_to, gedp->dbip->dbi_base2local);
    const std::string text_at_arg = point_arg(text_at, gedp->dbip->dbi_base2local);
    const char *direct_text_argv[] = {
	"annotate", "text", "--at", text_at_arg.c_str(), "--plane", "xy", "--frame",
	"--bold", "--text-height", DIRECT_TEXT_HEIGHT, "--color", "255/255/255",
	"component-title", "M35 OVERALL", NULL
    };
    const char *linear_argv[] = {
	"annotate", "dimension", "linear", "--from", linear_from_arg.c_str(), "--to",
	linear_to_arg.c_str(), "--offset", DIRECT_DIMENSION_OFFSET, "--text-height",
	DIRECT_TEXT_HEIGHT, "--color", "0/255/255",
	"--line-style", "dashed", "component-linear", NULL
    };
    const char *ordinate_argv[] = {
	"annotate", "dimension", "ordinate", "--origin", linear_from_arg.c_str(), "--to",
	ordinate_to_arg.c_str(), "--axis", "y", "--offset", DIRECT_DIMENSION_OFFSET,
	"--text-height", DIRECT_TEXT_HEIGHT, "--color", "255/80/255", "component-ordinate", NULL
    };
    const char *direct_ae_argv[] = {"ae", "45", "35", NULL};
    if (ged_exec_annotate(gedp, 5, hide_existing_argv) != BRLCAD_OK ||
	ged_exec_ae(gedp, 3, direct_ae_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 14, direct_text_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 17, linear_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 16, ordinate_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to create direct annotation examples: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    (void)ged_exec_autoview(gedp, 1, autoview_argv);
    (void)ged_exec_ae(gedp, 3, direct_ae_argv);
    if (generate && !commands_only)
	capture_image(gedp, 6);
    else if (!commands_only)
	ret += img_cmp(6, gedp, cache_dir, false, !keep_images, continue_on_failure,
	    ADIFF_THRESHOLD, "annotate_clear", "annotate");
    const fastf_t sphere_radius = std::min(bmax[X] - bmin[X], bmax[Z] - bmin[Z]) * 0.12;
    point_t sphere_center, sphere_rim, sphere_top, angular_from, angular_to;
    VSET(sphere_center, bmax[X] + sphere_radius * 2.5, bmax[Y], bmax[Z]);
    VSET(sphere_rim, sphere_center[X] + sphere_radius, sphere_center[Y], sphere_center[Z]);
    VSET(sphere_top, sphere_center[X], sphere_center[Y] + sphere_radius, sphere_center[Z]);
    VSET(angular_from, sphere_center[X] + sphere_radius, sphere_center[Y], sphere_center[Z]);
    VSET(angular_to, sphere_center[X], sphere_center[Y] + sphere_radius, sphere_center[Z]);
    struct rt_wdb *dimension_wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
    if (mk_sph(dimension_wdbp, "annotation-demo-sphere", sphere_center, sphere_radius))
	bu_exit(EXIT_FAILURE, "Unable to create direct-dimension geometry\n");
    const char *sphere_draw_argv[] = {"draw", "annotation-demo-sphere", NULL};
    const std::string center_arg = point_arg(sphere_center, gedp->dbip->dbi_base2local);
    const std::string rim_arg = point_arg(sphere_rim, gedp->dbip->dbi_base2local);
    const std::string top_arg = point_arg(sphere_top, gedp->dbip->dbi_base2local);
    const std::string angular_from_arg = point_arg(angular_from, gedp->dbip->dbi_base2local);
    const std::string angular_to_arg = point_arg(angular_to, gedp->dbip->dbi_base2local);
    const std::string sphere_offset_arg = std::to_string(sphere_radius *
	PRIMITIVE_DIMENSION_OFFSET_SCALE * gedp->dbip->dbi_base2local);
    const std::string angular_offset_arg = std::to_string(sphere_radius *
	PRIMITIVE_ANGULAR_OFFSET_SCALE * gedp->dbip->dbi_base2local);
    const char *radius_argv[] = {
	"annotate", "dimension", "radius", "--center", center_arg.c_str(), "--to",
	rim_arg.c_str(), "--offset", sphere_offset_arg.c_str(), "--text-height", PRIMITIVE_TEXT_HEIGHT,
	"--color", "255/80/255",
	"sphere-radius", NULL
    };
    const char *diameter_argv[] = {
	"annotate", "dimension", "diameter", "--center", center_arg.c_str(), "--to",
	top_arg.c_str(), "--offset", sphere_offset_arg.c_str(), "--text-height", PRIMITIVE_TEXT_HEIGHT,
	"--color", "255/160/0",
	"sphere-diameter", NULL
    };
    const char *angular_argv[] = {
	"annotate", "dimension", "angular", "--vertex", center_arg.c_str(), "--from",
	angular_from_arg.c_str(), "--to", angular_to_arg.c_str(), "--offset", angular_offset_arg.c_str(),
	"--text-height", PRIMITIVE_TEXT_HEIGHT, "--color", "80/255/80", "sphere-angle", NULL
    };
    if (ged_exec_draw(gedp, 2, sphere_draw_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 15, radius_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 15, diameter_argv) != BRLCAD_OK ||
	ged_exec_annotate(gedp, 17, angular_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to create radial annotation examples: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    const char *erase_component_argv[] = {
	"erase", "component", "component-title", "component-linear", "component-ordinate", NULL
    };
    if (ged_exec_erase(gedp, 5, erase_component_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to isolate direct-dimension geometry\n");
    if (ged_exec_draw(gedp, 2, sphere_draw_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to redraw direct-dimension geometry\n");
    (void)ged_exec_autoview(gedp, 1, autoview_argv);
    const char *primitive_ae_argv[] = {"ae", "0", "90", NULL};
    (void)ged_exec_ae(gedp, 3, primitive_ae_argv);
    if (generate && !commands_only)
	capture_image(gedp, 7);
    else if (!commands_only)
	ret += img_cmp(7, gedp, cache_dir, false, !keep_images, continue_on_failure,
	    ADIFF_THRESHOLD, "annotate_clear", "annotate");
    verify_radial_dimensions(gedp, sphere_center, sphere_radius);

    const fastf_t output_fill_half_size =
	std::min(bmax[X] - bmin[X], bmax[Z] - bmin[Z]) *
	OUTPUT_FILL_HALF_SIZE_SCALE;
    verify_output_consumers(gedp, bounds_center, output_fill_half_size);

    verify_geometry_update(gedp);
    bu_log("PASS annotate command, color, and update assertions\n");

    ged_close(gedp);
    bu_file_delete(working_db);
    return ret ? EXIT_FAILURE : EXIT_SUCCESS;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
