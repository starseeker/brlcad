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

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#include <bu.h>
#include <BObol/BDatabaseSource.h>
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
	"60 0 0", "annotate-update-leader", "UPDATED LEADER", NULL
    };
    if (ged_exec_annotate(gedp, 10, leader_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Unable to create leader update source: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    point_t initial_target;
    annotation_anchor(initial_target, gedp, "annotate-update-leader");

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

    const char *update_leader_argv[] = {
	"annotate", "update", "annotate-update-leader", NULL
    };
    if (ged_exec_annotate(gedp, 3, update_leader_argv) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "Leader did not update after a geometry change: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    point_t updated_target;
    annotation_anchor(updated_target, gedp, "annotate-update-leader");
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

    const char *parent_argv[] = {"g", "annotate-update-parent", dimension_name,
	"annotate-update-leader", NULL};
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
