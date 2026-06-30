/*                O B O L _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file obol_sync.cpp
 *
 * Tests libged's direct GED draw transaction to Obol controller bridge.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/init.h"
#include "brlobol/lod_realization.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/scene_controller.h"
#include "brlobol/scene_group.h"
#include "brlobol/vlist_shape.h"
#include "brlobol/view_controller.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/str.h"
#include "ged.h"
#include "ged/commands.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "rt/db_internal.h"
#include "rt/view.h"
#include "wdb.h"

#include "../../bsg_ged_draw_private.h"

#include <Inventor/SbMatrix.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/nodes/SoSeparator.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static int
make_obol_sync_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    point_t center = {5.0, 0.0, 0.0};
    point_t group_center = {-5.0, 0.0, 0.0};
    point_t draft_center = {0.0, 5.0, 0.0};
    point_t nested_leaf_center = {0.0, -5.0, 0.0};
    point_t nested_sibling_center = {3.0, -5.0, 0.0};
    point_t rename_center = {8.0, 3.0, 0.0};

    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0 &&
	mk_sph(wdbp, "ball.s", center, 1.0) == 0 &&
	mk_sph(wdbp, "group_only.s", group_center, 1.0) == 0 &&
	mk_sph(wdbp, "draft_move.s", draft_center, 1.0) == 0 &&
	mk_sph(wdbp, "nested_leaf.s", nested_leaf_center, 1.0) == 0 &&
	mk_sph(wdbp, "nested_sibling.s", nested_sibling_center, 1.0) == 0 &&
	mk_sph(wdbp, "rename_source.s", rename_center, 1.0) == 0;
    if (ret) {
	struct wmember child_wm;
	struct wmember parent_wm;
	BU_LIST_INIT(&child_wm.l);
	BU_LIST_INIT(&parent_wm.l);
	ret = mk_addmember("nested_leaf.s", &child_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_comb(wdbp, "nested_child.c", &child_wm.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) == 0 &&
	    mk_addmember("nested_child.c", &parent_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_addmember("nested_sibling.s", &parent_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_comb(wdbp, "nested_parent.c", &parent_wm.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) == 0;
    }
    if (ret) {
	struct rt_annot_internal ann;
	memset(&ann, 0, sizeof(ann));
	ann.magic = RT_ANNOT_INTERNAL_MAGIC;
	VSET(ann.V, 50.0, 0.0, 0.0);
	ann.vert_count = 2;
	ann.verts = (point2d_t *)bu_calloc(2, sizeof(point2d_t),
		"obol sync annot verts");
	V2SET(ann.verts[0], 0.0, 0.0);
	V2SET(ann.verts[1], 0.25, 0.5);
	ann.ant.count = 1;
	ann.ant.reverse = (int *)bu_calloc(1, sizeof(int),
		"obol sync annot reverse");
	ann.ant.segments = (void **)bu_calloc(1, sizeof(void *),
		"obol sync annot segments");
	struct line_seg *lsg;
	BU_ALLOC(lsg, struct line_seg);
	lsg->magic = CURVE_LSEG_MAGIC;
	lsg->start = 0;
	lsg->end = 1;
	ann.ant.segments[0] = (void *)lsg;
	ret = mk_annot(wdbp, "annot_line.s", &ann) == 0;
	BU_PUT(lsg, struct line_seg);
	bu_free(ann.ant.segments, "obol sync annot segments");
	bu_free(ann.ant.reverse, "obol sync annot reverse");
	bu_free(ann.verts, "obol sync annot verts");
    }
    wdb_close(wdbp);
    return ret;
}

static const char *
skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static int
path_equal(const char *a, const char *b)
{
    if (!a || !b)
	return 0;
    if (BU_STR_EQUAL(a, b))
	return 1;
    return BU_STR_EQUAL(skip_leading_slash(a), skip_leading_slash(b));
}

static SoBRLDatabaseSource *
source_for_path(SoBRLSceneController *controller, const char *path)
{
    if (!controller || !path)
	return NULL;
    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (source && path_equal(source->path.getValue().getString(), path))
	    return source;
    }
    return NULL;
}

struct record_source_state {
    int found;
    ged_draw_shape_ref ref;
    ged_draw_group_ref group;
    uint64_t sourceRevision;
    uint64_t inputsRevision;
    const char *matchPath;
    int visible;
    int highlighted;
    int drawMode;
    int lineWidth;
    fastf_t transparency;
};

static int
record_source_state_cb(const struct ged_draw_shape_record *record,
	void *userdata)
{
    record_source_state *state =
	static_cast<record_source_state *>(userdata);
    if (!state || !record)
	return 1;
    const char *target = state->matchPath ? state->matchPath : "box.s";
    if (!path_equal(record->display_name, target) &&
	    !path_equal(record->leaf_name, target))
	return 1;

    state->found = 1;
    state->ref = record->ref;
    state->group = record->group;
    state->sourceRevision = record->source_revision;
    state->inputsRevision = record->inputs_revision;
    state->visible = record->visible;
    state->highlighted = record->highlighted;
    state->drawMode = record->draw_mode;
    state->lineWidth = record->line_width;
    state->transparency = record->transparency;
    return 0;
}

struct group_source_state {
    int found;
    ged_draw_group_ref ref;
    const char *path;
    const char *matchPath;
};

static int
group_source_state_cb(const struct ged_draw_group_record *record,
	void *userdata)
{
    group_source_state *state =
	static_cast<group_source_state *>(userdata);
    if (!state || !record || !record->path)
	return 1;
    if (state->matchPath && !BU_STR_EQUAL(state->matchPath, record->path))
	return 1;

    state->found = 1;
    state->ref = record->ref;
    state->path = record->path;
    return 0;
}

int
main(int argc, char **argv)
{
    (void)argc;
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);
    brlobol_init(NULL);

    const char *dbpath = "ged_obol_draw_sync_tmp.g";
    bu_file_delete(dbpath);
    if (!make_obol_sync_db(dbpath))
	FAIL("failed to create GED Obol draw-sync test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open GED Obol draw-sync test database");

    SoSeparator *root = new SoSeparator;
    root->ref();
    SoBRLSceneController scene(root);
    if (!ged_draw_obol_scene_controller_attach(gedp, &scene, 0))
	FAIL("GED Obol scene-controller attachment should succeed");
    if (ged_draw_obol_scene_controller(gedp) != &scene)
	FAIL("GED should return the attached Obol scene controller");
    if (ged_draw_obol_controller(gedp))
	FAIL("direct scene-controller attachment should not report a view controller");

    const char *draw_box[2] = {"draw", "box.s"};
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK)
	FAIL("real GED draw command should succeed");
    if (scene.getDatabaseSourceCount() != 1 ||
	    !source_for_path(&scene, "box.s"))
	FAIL("attached Obol scene controller should mirror GED draw command");
    SoBRLDatabaseSource *box_source = source_for_path(&scene, "box.s");
    if (!box_source ||
	    box_source->realizationStatus.getValue() !=
	    SoBRLDatabaseSource::REALIZED ||
	    box_source->getRealizedShapeCount() <= 0)
	FAIL("mirrored GED draw should realize Obol wire geometry");

    ged_draw_obol_scene_controller_detach(gedp);
    if (ged_draw_obol_scene_controller(gedp) || ged_draw_obol_controller(gedp))
	FAIL("GED Obol scene-controller detach should clear borrowed pointers");
    const char *draw_ball[2] = {"draw", "ball.s"};
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("second real GED draw command should succeed");
    if (scene.getDatabaseSourceCount() != 1 ||
	    !source_for_path(&scene, "box.s") ||
	    source_for_path(&scene, "ball.s"))
	FAIL("detached scene controller should not mirror later GED draw commands");

    SoBRLSceneController *owned_scene =
	ged_draw_obol_scene_controller_ensure(gedp, 1);
    if (!owned_scene || ged_draw_obol_scene_controller(gedp) != owned_scene ||
	    !ged_draw_obol_scene_controller_owned(gedp) ||
	    ged_draw_obol_controller(gedp))
	FAIL("GED should create and report an owned Obol scene controller");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should full-sync current GED draw state");
    const char *draw_annot_line[2] = {"draw", "annot_line.s"};
    if (ged_exec_draw(gedp, 2, draw_annot_line) != BRLCAD_OK)
	FAIL("GED annotation draw should succeed for owned Obol publication");
    record_source_state annot_record = {0, GED_DRAW_SHAPE_REF_NULL_INIT,
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "annot_line.s", 0, 0, 0, 0, 0.0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &annot_record);
    if (!annot_record.found)
	FAIL("GED annotation draw should create a shape record");
    struct ged_draw_shape_geometry_summary annot_geometry;
    memset(&annot_geometry, 0, sizeof(annot_geometry));
    if (!ged_draw_shape_ref_geometry_summary(gedp, annot_record.ref,
	    &annot_geometry) ||
	    !annot_geometry.valid ||
	    !annot_geometry.geometry_name ||
	    !BU_STR_EQUAL(annot_geometry.geometry_name, "annotation") ||
	    annot_geometry.point_count != 2 ||
	    annot_geometry.index_count != 0)
	FAIL("GED annotation geometry summary should read owned Obol annotation VLIST");
    SoBRLDatabaseSource *annot_source =
	source_for_path(owned_scene, "annot_line.s");
    SoBRLVListShape *annot_shape = annot_source ?
	annot_source->getRealizedShape() : NULL;
    if (owned_scene->getDatabaseSourceCount() != 3 ||
	    !annot_shape ||
	    annot_shape->point.getNum() != 2 ||
	    annot_shape->command.getNum() != 2 ||
	    annot_shape->command[0] != SoBRLVListShape::MOVE ||
	    annot_shape->command[1] != SoBRLVListShape::DRAW ||
	    strcmp(annot_shape->sourceType.getValue().getString(),
		"annotation") != 0 ||
	    strcmp(annot_shape->geometryKind.getValue().getString(),
		"annotation") != 0 ||
	    fabs(annot_shape->point[1][0] - 50.25f) > 0.001f ||
	    fabs(annot_shape->point[1][1] - 0.5f) > 0.001f ||
	    fabs(annot_shape->point[1][2]) > 0.001f)
	FAIL("GED annotation draw should publish line segments into the owned Obol source");
    const char *erase_annot_line[2] = {"erase", "annot_line.s"};
    if (ged_exec_erase(gedp, 2, erase_annot_line) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "annot_line.s"))
	FAIL("GED annotation erase should restore the owned Obol source baseline");
    const char *draw_rename_source[2] = {"draw", "rename_source.s"};
    if (ged_exec_draw(gedp, 2, draw_rename_source) != BRLCAD_OK)
	FAIL("GED rename-source draw should succeed");
    if (!owned_scene->findDatabaseSource("rename_source.s") ||
	    owned_scene->getDatabaseSourceCount() != 3)
	FAIL("GED rename-source draw should create an owned Obol source");
    if (owned_scene->setDatabaseSourceState("rename_source.s",
	    TRUE,
	    9191,
	    7,
	    TRUE,
	    FALSE,
	    6,
	    11,
	    0.25f,
	    FALSE,
	    SbColor(1.0f, 1.0f, 1.0f),
	    FALSE,
	    SbColor(1.0f, 1.0f, 1.0f),
	    0) <= 0)
	FAIL("owned Obol rename-source state sentinel update should succeed");
    const char *move_rename_source[3] = {
	"move", "rename_source.s", "renamed_source.s"
    };
    if (ged_exec(gedp, 3, move_rename_source) != BRLCAD_OK)
	FAIL("GED move command should rename the drawn source");
    SoBRLDatabaseSource *renamed_source =
	owned_scene->findDatabaseSource("renamed_source.s");
    BRLObolDatabaseSourceSummary renamed_summary;
    if (source_for_path(owned_scene, "rename_source.s") ||
	    !renamed_source ||
	    owned_scene->getDatabaseSourceCount() != 3 ||
	    !renamed_source->getSummary(renamed_summary) ||
	    strcmp(renamed_summary.path.getString(), "renamed_source.s") != 0 ||
	    renamed_summary.lineWidth != 11 ||
	    renamed_summary.inputsRevision != 7 ||
	    renamed_summary.sourceRevision == 9191)
	FAIL("GED rename transaction should rename the owned Obol source in place");
    SoBRLVListShape *renamed_shape = renamed_source->getRealizedShape();
    if (!renamed_shape ||
	    !path_equal(renamed_shape->sourcePath.getValue().getString(),
		"renamed_source.s") ||
	    !path_equal(renamed_shape->ownerSourcePath.getValue().getString(),
		"renamed_source.s") ||
	    renamed_shape->lineWidth.getValue() != 11)
	FAIL("GED rename transaction should retarget owned Obol realized shape metadata");
    const char *erase_renamed_source[2] = {"erase", "renamed_source.s"};
    if (ged_exec_erase(gedp, 2, erase_renamed_source) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "renamed_source.s"))
	FAIL("GED renamed source erase should restore the owned Obol source baseline");
    if (owned_scene->findGroup("group_only.s"))
	FAIL("owned Obol group creation sentinel should not exist before lookup");
    struct db_full_path group_only_path;
    db_full_path_init(&group_only_path);
    if (db_string_to_path(&group_only_path, gedp->dbip,
	    "group_only.s") != 0)
	FAIL("owned Obol group creation sentinel path should resolve");
    ged_draw_group_ref group_only_ref =
	ged_draw_group_ref_lookup_or_create(gedp, &group_only_path);
    db_free_full_path(&group_only_path);
    if (ged_draw_group_ref_is_null(group_only_ref))
	FAIL("GED group lookup/create should return the sentinel group");
    if (!owned_scene->findGroup("group_only.s") ||
	    owned_scene->getGroupDatabaseSourceCount("group_only.s") != 0)
	FAIL("GED group lookup/create should ensure an empty owned Obol group");
    if (!path_equal(ged_draw_registry_group_ref_semantic_path(gedp,
	    group_only_ref), "group_only.s"))
	FAIL("GED group refs should retain a semantic registry path for Obol lookup");
    if (!owned_scene->ensureGroup("group_only.s/obol_child.s"))
	FAIL("owned Obol group child-count sentinel should be created");
    struct ged_draw_scene_tree_summary group_only_tree;
    memset(&group_only_tree, 0, sizeof(group_only_tree));
    if (!ged_draw_group_ref_tree_summary(gedp, group_only_ref,
	    &group_only_tree) ||
	    group_only_tree.child_count != 1)
	FAIL("GED group tree summaries should prefer owned Obol child counts");
    void *group_only_ctx = ged_draw_group_ref_context(gedp, group_only_ref);
    struct ged_draw_scene_tree_summary group_only_context_tree;
    memset(&group_only_context_tree, 0, sizeof(group_only_context_tree));
    if (!group_only_ctx ||
	    !ged_draw_scene_context_tree_summary(group_only_ctx,
		&group_only_context_tree) ||
	    group_only_context_tree.child_count != 1)
	FAIL("GED scene-context tree summaries should prefer owned Obol child counts");
    void *group_only_child_ctx =
	ged_draw_scene_context_child_at(group_only_ctx, 0);
    struct ged_draw_scene_tree_summary group_only_child_tree;
    memset(&group_only_child_tree, 0, sizeof(group_only_child_tree));
    if (!group_only_child_ctx ||
	    !ged_draw_scene_context_tree_summary(group_only_child_ctx,
		&group_only_child_tree) ||
	    !group_only_child_tree.valid ||
	    !group_only_child_tree.is_group ||
	    !path_equal(ged_draw_scene_context_name(group_only_child_ctx),
		"obol_child.s") ||
	    ged_draw_scene_context_parent(group_only_child_ctx) !=
		group_only_ctx)
	FAIL("GED scene-context child traversal should return owned Obol child contexts");
    if (owned_scene->setGroupDrawIntent("group_only.s",
	    "ged-draw-group:group_only.s", BRLOBOL_LOD_DRAW_WIRE,
	    BRLOBOL_LOD_DRAW_WIRE, TRUE, 0) < 0 ||
	    !ged_draw_group_context_is_overlay(group_only_ctx))
	FAIL("GED group contexts should read owned Obol overlay state");
    if (owned_scene->setGroupDrawIntent("group_only.s",
	    "ged-draw-group:group_only.s", BRLOBOL_LOD_DRAW_WIRE,
	    BRLOBOL_LOD_DRAW_WIRE, FALSE, 0) < 0 ||
	    ged_draw_group_context_is_overlay(group_only_ctx))
	FAIL("GED group contexts should clear owned Obol overlay state");
    if (owned_scene->removeGroup("group_only.s/obol_child.s") <= 0)
	FAIL("owned Obol group child-count sentinel should be removable");
    if (owned_scene->setGroupDrawIntent("group_only.s",
	    "ged-draw-group:group_only.s", BRLOBOL_LOD_DRAW_WIRE,
	    BRLOBOL_LOD_DRAW_WIRE, TRUE, 0) < 0)
	FAIL("owned Obol group overlay erase sentinel update should succeed");
    if (ged_draw_source_erase_groups_by_name_at_root(gedp, "group_only.s"))
	FAIL("GED group-name erase should treat owned Obol overlay state as authoritative");
    memset(&group_only_tree, 0, sizeof(group_only_tree));
    if (!ged_draw_group_ref_tree_summary(gedp, group_only_ref,
	    &group_only_tree) ||
	    !group_only_tree.valid ||
	    !group_only_tree.is_group)
	FAIL("GED source-adapter group summary should preserve overlay groups from owned Obol state");
    if (owned_scene->setGroupDrawIntent("group_only.s",
	    "ged-draw-group:group_only.s", BRLOBOL_LOD_DRAW_WIRE,
	    BRLOBOL_LOD_DRAW_WIRE, FALSE, 0) < 0)
	FAIL("owned Obol group overlay erase sentinel restore should succeed");
    const size_t original_root_child_count =
	ged_draw_source_root_child_count(gedp);
    if (!ged_draw_obol_group_ensure_for_path(gedp,
	    "__obol_root_count_only.s", "__obol_root_count_only.s",
	    GED_DRAW_MODE_WIRE, 0))
	FAIL("GED Obol root child-count sentinel group should be ensured");
    if (owned_scene->getGroupChildCount("/") !=
	    (int)(original_root_child_count + 1) ||
	    ged_draw_source_root_child_count(gedp) !=
	    original_root_child_count + 1)
	FAIL("GED source-root child count should prefer owned Obol root children");
    if (owned_scene->removeGroup("__obol_root_count_only.s") <= 0)
	FAIL("GED Obol root child-count sentinel group should be removable");
    if (owned_scene->setDatabaseSourceState("box.s",
	    TRUE,
	    77,
	    88,
	    FALSE,
	    TRUE,
	    3,
	    7,
	    0.375f,
	    FALSE,
	    SbColor(1.0f, 1.0f, 1.0f),
	    TRUE,
	    SbColor(0.2f, 0.4f, 0.6f),
	    1234) <= 0)
	FAIL("owned Obol source state sentinel update should succeed");
    record_source_state box_record = {0, GED_DRAW_SHAPE_REF_NULL_INIT,
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "box.s", 0, 0, 0, 0, 0.0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &box_record);
    if (!box_record.found ||
	    box_record.sourceRevision != 77 ||
	    box_record.inputsRevision != 88)
	FAIL("GED shape records should read source summary state from the owned Obol controller");
    if (!path_equal(ged_draw_registry_shape_ref_semantic_path(gedp,
	    box_record.ref), "box.s"))
	FAIL("GED shape refs should retain a semantic registry path for Obol lookup");
    void *box_ctx = ged_draw_shape_ref_context(gedp, box_record.ref);
    struct ged_draw_scene_tree_summary box_context_tree;
    memset(&box_context_tree, 0, sizeof(box_context_tree));
    if (!box_ctx ||
	    !ged_draw_shape_context_has_state(box_ctx) ||
	    !ged_draw_scene_context_tree_summary(box_ctx,
		&box_context_tree) ||
	    !box_context_tree.valid ||
	    !box_context_tree.is_group ||
	    box_context_tree.child_count <= 0)
	FAIL("GED shape-ref context should resolve to an owned Obol database-source context");
    if (ged_draw_shape_context_source(box_ctx) != box_ctx)
	FAIL("GED Obol shape-ref context source should be the owned Obol database-source context");
    struct ged_draw_database_source_summary box_context_source;
    memset(&box_context_source, 0, sizeof(box_context_source));
    if (!ged_draw_scene_context_source_summary(box_ctx,
		&box_context_source) ||
	    !box_context_source.valid ||
	    box_context_source.source_revision != 77 ||
	    box_context_source.inputs_revision != 88)
	FAIL("GED scene-context source summaries should read owned Obol source state");
    struct ged_draw_scene_display_summary box_context_display;
    memset(&box_context_display, 0, sizeof(box_context_display));
    if (!ged_draw_scene_context_display_summary(box_ctx,
	    &box_context_display) ||
	    !box_context_display.valid ||
	    box_context_display.visible ||
	    !box_context_display.highlighted ||
	    box_context_display.line_width != 7 ||
	    fabs(box_context_display.transparency - 0.375) > 0.001)
	FAIL("GED scene-context display summaries should read owned Obol display state");
    if (box_record.visible ||
	    !box_record.highlighted ||
	    box_record.drawMode != GED_DRAW_MODE_WIRE ||
	    box_record.lineWidth != 7 ||
	    fabs(box_record.transparency - 0.375) > 0.001)
	FAIL("GED shape records should read display state from the owned Obol controller");
    struct ged_draw_scene_display_summary box_display;
    memset(&box_display, 0, sizeof(box_display));
    if (!ged_draw_shape_ref_display_summary(gedp, box_record.ref,
	    &box_display) ||
	    !box_display.valid ||
	    box_display.visible ||
	    !box_display.highlighted ||
	    box_display.line_style != 3 ||
	    box_display.line_width != 7 ||
	    fabs(box_display.transparency - 0.375) > 0.001 ||
	    !box_display.material_valid ||
	    box_display.material_color[0] != 51 ||
	    box_display.material_color[1] != 102 ||
	    box_display.material_color[2] != 153)
	FAIL("GED shape display summary should read state from the owned Obol controller");
    struct ged_draw_shape_material_summary box_material;
    memset(&box_material, 0, sizeof(box_material));
    if (!ged_draw_shape_ref_material_summary(gedp, box_record.ref,
	    &box_material) ||
	    !box_material.valid ||
	    box_material.material_revision != 1234 ||
	    box_material.material_color[0] != 51 ||
	    box_material.material_color[1] != 102 ||
	    box_material.material_color[2] != 153)
	FAIL("GED shape material summary should read state from the owned Obol controller");
    box_source = owned_scene->findDatabaseSource("/box.s");
    if (!box_source)
	box_source = source_for_path(owned_scene, "box.s");
    if (!box_source)
	FAIL("owned Obol source should be available for geometry summary");
    SoBRLVListShape *box_shape = box_source->getRealizedShape();
    if (!box_shape)
	FAIL("owned Obol source should expose realized VLIST geometry");
    if (!ged_draw_shape_ref_set_evaluated_region(gedp, box_record.ref, 1) ||
	    box_shape->regionId.getValue() != 1)
	FAIL("GED evaluated-region setter should mutate owned Obol shape metadata");
    struct ged_draw_shape_record eval_record;
    memset(&eval_record, 0, sizeof(eval_record));
    if (!ged_draw_shape_record_get(gedp, box_record.ref, &eval_record) ||
	    eval_record.evaluated_region != 1)
	FAIL("GED shape records should read evaluated-region metadata from owned Obol shapes");
    if (!ged_draw_shape_ref_set_evaluated_region(gedp, box_record.ref, 0) ||
	    box_shape->regionId.getValue() != 0)
	FAIL("GED evaluated-region setter should clear owned Obol shape metadata");
    SbVec3f sentinel_points[2] = {
	SbVec3f(11.0f, 0.0f, 0.0f),
	SbVec3f(12.0f, 0.0f, 0.0f)
    };
    int32_t sentinel_commands[2] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW
    };
    box_shape->setLineSet(sentinel_points, sentinel_commands, 2);
    struct ged_draw_view_line_summary box_line;
    memset(&box_line, 0, sizeof(box_line));
    if (!ged_draw_shape_ref_line_summary(gedp, box_record.ref, &box_line) ||
	    !box_line.valid ||
	    box_line.point_count != 2)
	FAIL("GED shape line summary should read realized VLIST state from the owned Obol controller");
    struct ged_draw_view_line_summary box_context_line;
    memset(&box_context_line, 0, sizeof(box_context_line));
    if (!ged_draw_shape_context_line_summary(box_ctx, &box_context_line) ||
	    !box_context_line.valid ||
	    box_context_line.point_count != 2)
	FAIL("GED shape-context line summary should read realized VLIST state from the owned Obol controller");
    void *box_source_ctx = ged_draw_shape_context_source(box_ctx);
    if (!box_source_ctx)
	FAIL("GED shape contexts should expose a source context for Obol traversal");
    void *box_child_ctx = ged_draw_scene_context_child_at(box_source_ctx, 0);
    struct ged_draw_scene_tree_summary box_child_tree;
    memset(&box_child_tree, 0, sizeof(box_child_tree));
    if (!box_child_ctx)
	FAIL("GED source scene-context traversal should create owned Obol realized child contexts");
    if (!ged_draw_scene_context_tree_summary(box_child_ctx,
	    &box_child_tree) ||
	    !box_child_tree.valid)
	FAIL("GED source scene-context traversal should summarize owned Obol realized children");
    if (!box_child_tree.is_shape)
	FAIL("GED source scene-context traversal should classify owned Obol realized children as shapes");
    if (ged_draw_scene_context_parent(box_child_ctx) != box_source_ctx)
	FAIL("GED source scene-context traversal should return owned Obol realized children");
    point_t box_line_point;
    if (!ged_draw_shape_ref_line_point_at(gedp, box_record.ref, 1,
	    box_line_point) ||
	    fabs(box_line_point[0] - 12.0) > 0.001 ||
	    fabs(box_line_point[1]) > 0.001 ||
	    fabs(box_line_point[2]) > 0.001)
	FAIL("GED shape line point readback should read realized VLIST points from the owned Obol controller");
    point_t box_context_line_point;
    if (!ged_draw_shape_context_line_point_at(box_ctx, 1,
	    box_context_line_point) ||
	    fabs(box_context_line_point[0] - 12.0) > 0.001 ||
	    fabs(box_context_line_point[1]) > 0.001 ||
	    fabs(box_context_line_point[2]) > 0.001)
	FAIL("GED shape-context line point readback should read realized VLIST points from the owned Obol controller");
    int box_line_command = -1;
    if (!ged_draw_shape_ref_line_command_at(gedp, box_record.ref, 1,
	    &box_line_command) ||
	    box_line_command != GED_DRAW_VIEW_LINE_DRAW)
	FAIL("GED shape line command readback should read realized VLIST commands from the owned Obol controller");
    int box_context_line_command = -1;
    if (!ged_draw_shape_context_line_command_at(box_ctx, 1,
	    &box_context_line_command) ||
	    box_context_line_command != GED_DRAW_VIEW_LINE_DRAW)
	FAIL("GED shape-context line command readback should read realized VLIST commands from the owned Obol controller");
    point_t box_last_point;
    if (!ged_draw_shape_ref_last_point(gedp, box_record.ref,
	    box_last_point) ||
	    fabs(box_last_point[0] - 12.0) > 0.001 ||
	    fabs(box_last_point[1]) > 0.001 ||
	    fabs(box_last_point[2]) > 0.001)
	FAIL("GED shape last-point readback should read realized VLIST points from the owned Obol controller");
    struct ged_draw_shape_geometry_summary box_geometry;
    memset(&box_geometry, 0, sizeof(box_geometry));
    if (!ged_draw_shape_ref_geometry_summary(gedp, box_record.ref,
	    &box_geometry) ||
	    !box_geometry.valid ||
	    !box_geometry.geometry_name ||
	    !BU_STR_EQUAL(box_geometry.geometry_name, "line-set") ||
	    box_geometry.point_count != 2 ||
	    box_geometry.index_count != 0)
	FAIL("GED shape geometry summary should read realized geometry from the owned Obol controller");
    struct ged_draw_shape_geometry_summary box_context_geometry;
    memset(&box_context_geometry, 0, sizeof(box_context_geometry));
    if (!ged_draw_shape_context_geometry_summary(box_ctx,
	    &box_context_geometry) ||
	    !box_context_geometry.valid ||
	    !box_context_geometry.geometry_name ||
	    !BU_STR_EQUAL(box_context_geometry.geometry_name, "line-set") ||
	    box_context_geometry.point_count != 2 ||
	    box_context_geometry.index_count != 0)
	FAIL("GED shape-context geometry summary should read realized geometry from the owned Obol controller");
	    vect_t draw_bounds_min;
	    vect_t draw_bounds_max;
	    if (ged_draw_bounds(gedp, &draw_bounds_min, &draw_bounds_max, 0) ||
		    draw_bounds_max[0] < 11.9)
		FAIL("GED draw bounds should read database-source bounds from the owned Obol controller");
	    vect_t context_bounds_min;
	    vect_t context_bounds_max;
	    if (ged_draw_scene_context_subtree_bounds(box_ctx,
		    &context_bounds_min, &context_bounds_max, 0) ||
		    context_bounds_max[0] < 11.9)
		FAIL("GED source contexts should read owned Obol subtree bounds");
	    if (owned_scene->moveDatabaseSourceToGroup("box.s",
		    "group_only.s") != 1)
		FAIL("owned Obol source should move to the sentinel group for context bounds");
	    if (owned_scene->setGroupDrawIntent("group_only.s",
		    "ged-draw-group:group_only.s", BRLOBOL_LOD_DRAW_WIRE,
		    BRLOBOL_LOD_DRAW_WIRE, TRUE, 0) < 0 ||
		    !ged_draw_group_context_is_overlay(group_only_ctx))
		FAIL("owned Obol group overlay state should remain authoritative before bounds");
	    if (!ged_draw_scene_context_subtree_bounds(group_only_ctx,
		    &context_bounds_min, &context_bounds_max, 0))
		FAIL("GED group context bounds should skip owned Obol overlay groups");
	    if (owned_scene->setGroupDrawIntent("group_only.s",
		    "ged-draw-group:group_only.s", BRLOBOL_LOD_DRAW_WIRE,
		    BRLOBOL_LOD_DRAW_WIRE, FALSE, 0) < 0 ||
		    ged_draw_group_context_is_overlay(group_only_ctx))
		FAIL("owned Obol group overlay state should clear before bounds");
	    if (ged_draw_scene_context_subtree_bounds(group_only_ctx,
		    &context_bounds_min, &context_bounds_max, 0) ||
		    context_bounds_max[0] < 11.9)
		FAIL("GED group contexts should read owned Obol subtree bounds");
	    if (owned_scene->moveDatabaseSourceToGroup("box.s", "/") != 1)
		FAIL("owned Obol source should move back to the root group after context bounds");
	    vect_t box_xlate;
	    VSET(box_xlate, 5.0, 0.0, 0.0);
	    if (!ged_draw_shape_ref_translate_geometry(gedp, box_record.ref,
		    box_xlate))
		FAIL("GED shape geometry translation should succeed");
	    if (!ged_draw_shape_ref_line_point_at(gedp, box_record.ref, 1,
		    box_line_point) ||
		    fabs(box_line_point[0] - 17.0) > 0.001 ||
		    fabs(box_line_point[1]) > 0.001 ||
		    fabs(box_line_point[2]) > 0.001)
		FAIL("GED shape translation should mutate owned Obol VLIST points");
	    if (ged_draw_bounds(gedp, &draw_bounds_min, &draw_bounds_max, 0) ||
		    draw_bounds_max[0] < 16.9)
		FAIL("GED draw bounds should reflect translated owned Obol VLIST points");
	    point_t published_points[3] = {
		{21.0, 0.0, 0.0},
		{22.0, 1.0, 0.0},
		{23.0, 0.0, 0.0}
	    };
	    int published_commands[3] = {
		GED_DRAW_VIEW_LINE_MOVE,
		GED_DRAW_VIEW_LINE_DRAW,
		GED_DRAW_VIEW_LINE_DRAW
	    };
	    if (!ged_draw_shape_ref_publish_line_set(gedp, box_record.ref,
		    (const point_t *)published_points, published_commands, 3))
		FAIL("GED shape line-set publish should succeed");
	    if (!ged_draw_shape_ref_line_summary(gedp, box_record.ref,
		    &box_line) ||
		    !box_line.valid ||
		    box_line.point_count != 3)
		FAIL("GED shape line-set publish should update owned Obol VLIST count");
	    if (!ged_draw_shape_ref_line_point_at(gedp, box_record.ref, 2,
		    box_line_point) ||
		    fabs(box_line_point[0] - 23.0) > 0.001 ||
		    fabs(box_line_point[1]) > 0.001 ||
		    fabs(box_line_point[2]) > 0.001)
		FAIL("GED shape line-set publish should update owned Obol VLIST points");
	    if (ged_draw_bounds(gedp, &draw_bounds_min, &draw_bounds_max, 0) ||
		    draw_bounds_max[0] < 22.9)
		FAIL("GED draw bounds should reflect published owned Obol VLIST points");
	    point_t explicit_center;
	    VSET(explicit_center, 30.0, 31.0, 32.0);
	    if (!ged_draw_shape_ref_set_center(gedp, box_record.ref,
		    explicit_center))
		FAIL("GED shape center setter should succeed");
	    SbVec3f obol_center = box_shape->drawCenter.getValue();
	    if (!box_shape->drawCenterValid.getValue() ||
		    fabs(obol_center[0] - 30.0f) > 0.001f ||
		    fabs(obol_center[1] - 31.0f) > 0.001f ||
		    fabs(obol_center[2] - 32.0f) > 0.001f)
		FAIL("GED shape center setter should mutate owned Obol VLIST center");
	    BRLObolDatabaseSourceSummary box_placement_summary;
	    if (!box_source->getSummary(box_placement_summary) ||
		    !box_placement_summary.valid ||
		    !box_placement_summary.drawCenterValid ||
		    fabs(box_placement_summary.drawCenter[0] - 30.0f) >
			0.001f ||
		    fabs(box_placement_summary.drawCenter[1] - 31.0f) >
			0.001f ||
		    fabs(box_placement_summary.drawCenter[2] - 32.0f) >
			0.001f)
		FAIL("GED shape center setter should update owned Obol source placement");
	    int bounds_bad_cmd = -1;
	    if (!ged_draw_shape_ref_update_bounds_from_geometry(gedp,
		    box_record.ref, &bounds_bad_cmd) ||
		    bounds_bad_cmd != 0)
		FAIL("GED shape bounds update should succeed for owned Obol VLIST");
	    obol_center = box_shape->drawCenter.getValue();
	    if (!box_shape->drawCenterValid.getValue() ||
		    fabs(obol_center[0] - 22.0f) > 0.001f ||
		    fabs(obol_center[1] - 0.5f) > 0.001f ||
		    fabs(obol_center[2]) > 0.001f ||
		    !box_shape->drawSizeValid.getValue() ||
		    fabs(box_shape->drawSize.getValue() - 2.0f) > 0.001f)
		FAIL("GED shape bounds update should mutate owned Obol VLIST bounds metadata");
	    if (!box_source->getSummary(box_placement_summary) ||
		    !box_placement_summary.valid ||
		    !box_placement_summary.drawCenterValid ||
		    fabs(box_placement_summary.drawCenter[0] - 22.0f) >
			0.001f ||
		    fabs(box_placement_summary.drawCenter[1] - 0.5f) >
			0.001f ||
		    fabs(box_placement_summary.drawCenter[2]) > 0.001f ||
		    !box_placement_summary.drawSizeValid ||
		    fabs(box_placement_summary.drawSize - 2.0f) > 0.001f)
		FAIL("GED shape bounds update should update owned Obol source placement");
	    if (!ged_draw_shape_ref_geometry_clear(gedp, box_record.ref))
		FAIL("GED shape geometry clear should succeed");
	    if (!ged_draw_shape_ref_line_summary(gedp, box_record.ref,
		    &box_line) ||
		    !box_line.valid ||
		    box_line.point_count != 0)
		FAIL("GED shape geometry clear should clear owned Obol VLIST points");
	    if (!ged_draw_shape_ref_geometry_summary(gedp, box_record.ref,
		    &box_geometry) ||
		    !box_geometry.valid ||
		    box_geometry.point_count != 0)
		FAIL("GED shape geometry clear should update owned Obol geometry summary");
		    if (ged_draw_bounds(gedp, &draw_bounds_min, &draw_bounds_max, 0) ||
			    draw_bounds_max[0] > 20.0)
			FAIL("GED draw bounds should reflect cleared owned Obol VLIST points");
		    point_t aux_points[2] = {
			{26.0, 0.0, 0.0},
			{27.0, 1.0, 0.0}
		    };
		    int aux_commands[2] = {
			GED_DRAW_VIEW_LINE_MOVE,
			GED_DRAW_VIEW_LINE_DRAW
		    };
		    if (!ged_draw_obol_database_source_publish_auxiliary_line_set_for_path(
			    gedp, "box.s", "obol_aux_clear_sentinel",
			    (const point_t *)aux_points, aux_commands, 2))
			FAIL("GED Obol auxiliary line-set bridge should publish to owned source");
		    SoBRLVListShape *aux_shape =
			box_source->findAuxiliaryVListShape(
				"obol_aux_clear_sentinel");
		    if (!aux_shape ||
			    aux_shape->point.getNum() != 2 ||
			    aux_shape->command.getNum() != 2 ||
			    strcmp(aux_shape->recordRole.getValue().getString(),
				"auxiliary") != 0)
			FAIL("GED Obol auxiliary line-set bridge should create an auxiliary VLIST");
		    if (!ged_draw_shape_ref_geometry_clear(gedp, box_record.ref) ||
			    box_source->findAuxiliaryVListShape(
				"obol_aux_clear_sentinel"))
			FAIL("GED shape geometry clear should clear owned Obol auxiliary VLISTs");
		    point_t point_set_points[2] = {
			{24.0, 0.0, 0.0},
			{25.0, 1.0, 0.0}
	    };
	    if (!ged_draw_shape_ref_publish_point_set(gedp, box_record.ref,
		    (const point_t *)point_set_points, 2))
		FAIL("GED shape point-set publish should succeed");
	    if (!box_shape ||
		    box_shape->point.getNum() != 2 ||
		    box_shape->command.getNum() != 2 ||
		    box_shape->command[0] != SoBRLVListShape::POINT ||
		    box_shape->command[1] != SoBRLVListShape::POINT)
		FAIL("GED shape point-set publish should mutate owned Obol point commands");
	    int point_set_command = -1;
	    if (!ged_draw_shape_ref_line_command_at(gedp, box_record.ref, 0,
		    &point_set_command) ||
		    point_set_command != GED_DRAW_VIEW_LINE_POINT_DRAW)
		FAIL("GED shape point-set command readback should use owned Obol POINT commands");
	    if (!ged_draw_shape_ref_geometry_summary(gedp, box_record.ref,
		    &box_geometry) ||
		    !box_geometry.valid ||
		    !box_geometry.geometry_name ||
		    !BU_STR_EQUAL(box_geometry.geometry_name, "point-set") ||
		    box_geometry.point_count != 2 ||
		    box_geometry.index_count != 0)
		FAIL("GED shape point-set geometry summary should read owned Obol point-set publication");
	    if (!ged_draw_shape_ref_geometry_clear(gedp, box_record.ref))
		FAIL("GED shape point-set geometry clear should succeed");
	    point_t mesh_points[4] = {
		{31.0, 0.0, 0.0},
		{32.0, 0.0, 0.0},
		{31.0, 1.0, 0.0},
		{32.0, 1.0, 0.0}
	    };
	    vect_t mesh_normals[4] = {
		{0.0, 0.0, 1.0},
		{0.0, 0.0, 1.0},
		{0.0, 0.0, 1.0},
		{0.0, 0.0, 1.0}
	    };
	    int mesh_indices[5] = {0, 1, 3, 2, -1};
	    if (!ged_draw_shape_ref_publish_indexed_face_set(gedp,
		    box_record.ref, (const point_t *)mesh_points, 4,
		    (const vect_t *)mesh_normals, 4, mesh_indices, 5))
		FAIL("GED shape indexed-face publish should succeed");
	    SoBRLMeshShape *box_mesh = box_source->getRealizedMesh();
	    if (!box_mesh ||
		    box_mesh->point.getNum() != 4 ||
		    box_mesh->coordIndex.getNum() != 6 ||
		    fabs(box_mesh->point[3][0] - 32.0f) > 0.001f ||
		    box_mesh->coordIndex[4] != 3)
		FAIL("GED shape indexed-face publish should mutate owned Obol mesh fields");
	    if (!ged_draw_shape_ref_geometry_summary(gedp, box_record.ref,
		    &box_geometry) ||
		    !box_geometry.valid ||
		    !box_geometry.geometry_name ||
		    !BU_STR_EQUAL(box_geometry.geometry_name,
			"indexed-face-set") ||
		    box_geometry.point_count != 4 ||
		    box_geometry.index_count != 6)
		FAIL("GED shape geometry summary should read owned Obol mesh publication");
	    if (ged_draw_bounds(gedp, &draw_bounds_min, &draw_bounds_max, 0) ||
		    draw_bounds_max[0] < 31.9)
		FAIL("GED draw bounds should reflect published owned Obol mesh points");
	    if (!ged_draw_shape_ref_geometry_clear(gedp, box_record.ref))
		FAIL("GED shape geometry clear should clear owned Obol mesh");
	    if (box_mesh->point.getNum() != 0 ||
		    box_mesh->coordIndex.getNum() != 0)
		FAIL("GED shape geometry clear should empty owned Obol mesh fields");
	    if (!ged_draw_shape_ref_geometry_summary(gedp, box_record.ref,
		    &box_geometry) ||
		    !box_geometry.valid ||
		    box_geometry.point_count != 0 ||
		    box_geometry.index_count != 0)
		FAIL("GED shape geometry summary should reflect cleared owned Obol mesh");
	    struct directory *box_dp = db_lookup(gedp->dbip, "box.s",
		    LOOKUP_QUIET);
	    if (!box_dp)
		FAIL("box.s should be available for primitive wireframe publication");
	    struct rt_db_internal box_intern;
	    RT_DB_INTERNAL_INIT(&box_intern);
	    mat_t identity;
	    MAT_IDN(identity);
	    if (rt_db_get_internal(&box_intern, box_dp, gedp->dbip,
		    identity) < 0)
		FAIL("box.s internal lookup should succeed");
	    if (ged_draw_shape_ref_publish_primitive_wireframe(gedp,
		    box_record.ref, &box_intern, NULL, NULL, NULL, 0) < 0)
		FAIL("GED primitive wireframe publication should succeed");
	    rt_db_free_internal(&box_intern);
	    box_record.found = 0;
	    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
		    &box_record);
	    if (!box_record.found)
		FAIL("GED shape record should refresh after primitive wireframe publication");
	    if (!ged_draw_shape_ref_line_summary(gedp, box_record.ref,
		    &box_line) ||
		    !box_line.valid ||
		    box_line.point_count == 0)
		FAIL("GED primitive wireframe publication should update owned Obol VLIST");
	    if (!ged_draw_shape_ref_geometry_summary(gedp, box_record.ref,
		    &box_geometry) ||
		    !box_geometry.valid ||
		    !box_geometry.geometry_name ||
		    !BU_STR_EQUAL(box_geometry.geometry_name, "line-set") ||
		    box_geometry.point_count == 0)
		FAIL("GED primitive wireframe publication should publish owned Obol line geometry");
	    if (!ged_draw_shape_ref_geometry_clear(gedp, box_record.ref))
		FAIL("GED shape geometry clear should clear primitive wireframe publication");
	    const char *draw_box_shaded[3] = {"draw", "-m2", "box.s"};
	    if (ged_exec_draw(gedp, 3, draw_box_shaded) != BRLCAD_OK)
		FAIL("GED shaded draw should succeed for source realization test");
	    box_record.found = 0;
	    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
		    &box_record);
	    if (!box_record.found)
		FAIL("GED shape record should refresh after shaded draw source realization");
	    box_source = owned_scene->findDatabaseSource("/box.s");
	    if (!box_source)
		box_source = source_for_path(owned_scene, "box.s");
	    if (!box_source)
		FAIL("owned Obol source should remain available after shaded draw source realization");
	    box_mesh = box_source->getRealizedMesh();
	    if (!box_mesh ||
		    box_mesh->point.getNum() == 0 ||
		    box_mesh->coordIndex.getNum() == 0)
		FAIL("GED shaded source realization should publish owned Obol mesh fields");
	    BRLObolDatabaseSourceSummary box_realized_summary;
	    if (!box_source->getSummary(box_realized_summary) ||
		    box_realized_summary.stale ||
		    box_realized_summary.staleReason !=
		    SoBRLDatabaseSource::STALE_NONE ||
		    box_realized_summary.realizationStatus !=
		    SoBRLDatabaseSource::REALIZED ||
		    box_realized_summary.realizedSourceRevision !=
		    box_realized_summary.sourceRevision ||
		    box_realized_summary.realizedInputsRevision !=
		    box_realized_summary.inputsRevision)
		FAIL("GED shaded source realization should update owned Obol realization status");
	    if (!ged_draw_shape_ref_geometry_summary(gedp, box_record.ref,
		    &box_geometry) ||
		    !box_geometry.valid ||
		    !box_geometry.geometry_name ||
		    !BU_STR_EQUAL(box_geometry.geometry_name,
			"indexed-face-set") ||
		    box_geometry.point_count == 0 ||
		    box_geometry.index_count == 0)
		FAIL("GED shaded source realization should publish owned Obol mesh summary");
	    if (!ged_draw_group_ref_set_mode(gedp, box_record.group,
		    GED_DRAW_MODE_WIRE))
		FAIL("GED group wire mode restore should succeed after source realization test");
	    mat_t source_placement_mat;
	    MAT_IDN(source_placement_mat);
	    source_placement_mat[MDX] = 24.0;
	    source_placement_mat[MDY] = 1.0;
	    source_placement_mat[MDZ] = 2.0;
	    point_t source_placement_center = {3.0, 4.0, 5.0};
	    SbMatrix source_placement_sb = SbMatrix::identity();
	    source_placement_sb.setTranslate(SbVec3f(24.0f, 1.0f, 2.0f));
	    if (!ged_draw_obol_database_source_set_placement_for_path(gedp,
		    "box.s", 1, source_placement_mat, 1,
		    source_placement_center, 1, 12.5))
		FAIL("GED Obol source placement bridge should set owned source state");
	    if (!box_source->getSummary(box_realized_summary) ||
		    !box_realized_summary.drawMatrixValid ||
		    !box_realized_summary.drawMatrix.equals(source_placement_sb,
			0.0001f) ||
		    !box_realized_summary.drawCenterValid ||
		    fabs(box_realized_summary.drawCenter[0] - 3.0f) > 0.001f ||
		    fabs(box_realized_summary.drawCenter[1] - 4.0f) > 0.001f ||
		    fabs(box_realized_summary.drawCenter[2] - 5.0f) > 0.001f ||
		    !box_realized_summary.drawSizeValid ||
		    fabs(box_realized_summary.drawSize - 12.5f) > 0.001f)
		FAIL("GED Obol source placement bridge should update source summary");
	    box_mesh = box_source->getRealizedMesh();
	    if (!box_mesh ||
		    !box_mesh->drawMatrixValid.getValue() ||
		    !box_mesh->drawMatrix.getValue().equals(
			source_placement_sb, 0.0001f) ||
		    !box_mesh->drawCenterValid.getValue() ||
		    fabs(box_mesh->drawCenter.getValue()[0] - 3.0f) >
			0.001f ||
		    !box_mesh->drawSizeValid.getValue() ||
		    fabs(box_mesh->drawSize.getValue() - 12.5f) > 0.001f)
		FAIL("GED Obol source placement bridge should sync realized shape placement");
	    std::string box_source_path =
		box_source->path.getValue().getString();
	    if (owned_scene->setDatabaseSourcePlacementState(
		    box_source_path.c_str(), FALSE, SbMatrix::identity(),
		    FALSE, SbVec3f(0.0f, 0.0f, 0.0f), FALSE, 0.0f) < 0)
		FAIL("owned Obol source placement reset should succeed");
	    SbMatrix obol_draw_matrix = SbMatrix::identity();
	    obol_draw_matrix.setTranslate(SbVec3f(40.0f, 0.0f, 0.0f));
	    struct ged_draw_obol_database_source_record obol_draw_record;
	    if (!ged_draw_obol_database_source_record_for_path(gedp, "box.s",
		    &obol_draw_record))
		FAIL("GED Obol draw-state bridge should read source record before redraw");
	    obol_draw_record.draw_mode = GED_DRAW_MODE_WIRE;
	    if (!ged_draw_obol_database_source_apply_record_for_path(gedp,
		    "box.s", &obol_draw_record))
		FAIL("GED Obol draw-state bridge should set source wire mode before redraw");
	    mat_t obol_draw_mat;
	    MAT_IDN(obol_draw_mat);
	    obol_draw_mat[MDX] = 40.0;
	    if (!ged_draw_obol_database_source_set_placement_for_path(gedp,
		    "box.s", 1, obol_draw_mat, 0, NULL, 0, 0.0))
		FAIL("GED Obol draw-state bridge should set owned source draw matrix");
	    box_source->lineStyle = 5;
	    struct ged_draw_obol_draw_state_summary obol_draw_state;
	    if (!ged_draw_obol_database_source_draw_state_for_path(gedp,
		    "box.s", &obol_draw_state) ||
		    !obol_draw_state.valid ||
		    !obol_draw_state.draw_mode_valid ||
		    obol_draw_state.draw_mode != GED_DRAW_MODE_WIRE ||
		    obol_draw_state.line_style != 5 ||
		    !obol_draw_state.draw_mat_valid ||
		    fabs(obol_draw_state.draw_mat[MDX] - 40.0) > 0.001 ||
		    fabs(obol_draw_state.draw_mat[MDY]) > 0.001 ||
		    fabs(obol_draw_state.draw_mat[MDZ]) > 0.001)
		FAIL("GED Obol draw-state bridge should read owned line style and draw matrix");
	    if (ged_draw_shape_ref_redraw_wireframe(gedp, box_record.ref,
		    gedp->dbip, NULL, NULL, NULL, 0) < 0)
		FAIL("GED wire redraw should succeed with owned Obol draw matrix");
	    box_record.found = 0;
	    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
		    &box_record);
	    if (!box_record.found)
		FAIL("GED shape record should refresh after owned Obol draw matrix redraw");
	    box_source = owned_scene->findDatabaseSource("/box.s");
	    if (!box_source)
		box_source = source_for_path(owned_scene, "box.s");
	    if (!box_source)
		FAIL("owned Obol source should remain available after draw matrix redraw");
	    box_mesh = box_source->getRealizedMesh();
	    if (box_mesh) {
		box_mesh->drawMatrixValid = FALSE;
		box_mesh->drawMatrix = SbMatrix::identity();
	    }
	    if (!ged_draw_shape_ref_line_summary(gedp, box_record.ref,
		    &box_line) ||
		    !box_line.valid ||
		    box_line.point_count == 0)
		FAIL("GED redraw should publish wire geometry after owned Obol draw matrix readback");
	    if (ged_draw_bounds(gedp, &draw_bounds_min, &draw_bounds_max, 0) ||
		    draw_bounds_max[0] < 40.9)
		FAIL("GED redraw should use owned Obol draw matrix");
	    if (owned_scene->setDatabaseSourcePlacementState(
		    box_source_path.c_str(), FALSE, SbMatrix::identity(),
		    FALSE, SbVec3f(0.0f, 0.0f, 0.0f), FALSE, 0.0f) < 0)
		FAIL("owned Obol source draw matrix reset should succeed");
	    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK)
		FAIL("GED wire redraw should succeed for LoD policy test");
	    box_record.found = 0;
	    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
		    &box_record);
	    if (!box_record.found)
		FAIL("GED shape record should refresh before LoD policy test");
	    box_source = owned_scene->findDatabaseSource("/box.s");
	    if (!box_source)
		box_source = source_for_path(owned_scene, "box.s");
	    if (!box_source)
		FAIL("owned Obol source should remain available before LoD policy test");
	    void *lod_view_ctx = ged_draw_shape_ref_view_context(gedp,
		    box_record.ref);
	    if (!lod_view_ctx)
		FAIL("GED LoD view context should be available");
	    if (!ged_view_context_scale_set(lod_view_ctx, 7.0))
		FAIL("GED LoD view context scale should be settable");
	    struct rt_view_lod_policy lod_policy = RT_VIEW_LOD_POLICY_INIT;
	    lod_policy.csg_enabled = 1;
	    lod_policy.mesh_enabled = 0;
	    lod_policy.bot_threshold = 77;
	    lod_policy.curve_scale = 2.25;
	    lod_policy.point_scale = 3.25;
	    if (!ged_view_context_lod_policy_apply(lod_view_ctx, &lod_policy))
		FAIL("GED LoD view policy should be settable");
	    void *lod_view_ctxs[1] = {lod_view_ctx};
	    if (!ged_draw_shape_ref_lod_ensure(gedp, box_record.ref,
		    lod_view_ctx, lod_view_ctxs, 1))
		FAIL("GED LoD ensure should succeed for Obol source policy test");
	    if (!box_source->getSummary(box_realized_summary) ||
		    !(box_realized_summary.realizationRoleFlags &
			SoBRLDatabaseSource::REALIZATION_ROLE_CSG) ||
		    !box_realized_summary.realizationViewDependent ||
		    fabs(box_realized_summary.realizationViewScale - 7.0f) >
			0.001 ||
		    box_realized_summary.realizationBotThreshold != 77 ||
		    fabs(box_realized_summary.realizationCurveScale - 2.25f) >
			0.001 ||
		    fabs(box_realized_summary.realizationPointScale - 3.25f) >
			0.001) {
		fprintf(stderr,
			"LoD Obol summary role=%d view=%d scale=%g bot=%u curve=%g point=%g\n",
			box_realized_summary.realizationRoleFlags,
			box_realized_summary.realizationViewDependent ? 1 : 0,
			(double)box_realized_summary.realizationViewScale,
			(unsigned)box_realized_summary.realizationBotThreshold,
			(double)box_realized_summary.realizationCurveScale,
			(double)box_realized_summary.realizationPointScale);
		FAIL("GED LoD realization policy should update owned Obol source state");
	    }
	    box_record.found = 0;
	    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
		    &box_record);
	    if (!box_record.found)
		FAIL("GED shape record should refresh after LoD policy test");
	    if (!ged_draw_shape_ref_geometry_clear(gedp, box_record.ref))
		FAIL("GED shape geometry clear should clear shaded source realization mesh");
	    if (!ged_draw_shape_ref_set_visible(gedp, box_record.ref, 1))
		FAIL("GED visible setter should succeed");
    if (!ged_draw_shape_set_highlighted(gedp, box_record.ref, 0))
	FAIL("GED highlighted setter should succeed");
    const unsigned char override_color[3] = {10, 20, 30};
    if (!ged_draw_shape_ref_set_color(gedp, box_record.ref, override_color))
	FAIL("GED color setter should succeed");
    box_source = source_for_path(owned_scene, "box.s");
    BRLObolDatabaseSourceSummary box_source_summary;
    if (!box_source ||
	    !box_source->getSummary(box_source_summary) ||
	    !box_source_summary.visible ||
	    box_source_summary.highlighted ||
	    !box_source_summary.colorOverride ||
	    fabsf(box_source_summary.color[0] -
		(10.0f / 255.0f)) > 1.0e-6f ||
	    fabsf(box_source_summary.color[1] -
		(20.0f / 255.0f)) > 1.0e-6f ||
	    fabsf(box_source_summary.color[2] -
		(30.0f / 255.0f)) > 1.0e-6f)
	FAIL("GED visible/highlight/color setters should mutate the owned Obol source");
    int display_synced = 0;
    SoBRLVListShape *display_vlist = box_source->getRealizedShape();
    if (display_vlist) {
	SbColor display_color = display_vlist->color.getValue();
	display_synced =
	    display_vlist->visible.getValue() &&
	    !display_vlist->highlighted.getValue() &&
	    display_vlist->colorOverride.getValue() &&
	    fabsf(display_color[0] - (10.0f / 255.0f)) < 1.0e-6f &&
	    fabsf(display_color[1] - (20.0f / 255.0f)) < 1.0e-6f &&
	    fabsf(display_color[2] - (30.0f / 255.0f)) < 1.0e-6f;
    }
    SoBRLMeshShape *display_mesh = box_source->getRealizedMesh();
    if (!display_synced && display_mesh) {
	SbColor display_color = display_mesh->color.getValue();
	display_synced =
	    display_mesh->visible.getValue() &&
	    !display_mesh->highlighted.getValue() &&
	    display_mesh->colorOverride.getValue() &&
	    fabsf(display_color[0] - (10.0f / 255.0f)) < 1.0e-6f &&
	    fabsf(display_color[1] - (20.0f / 255.0f)) < 1.0e-6f &&
	    fabsf(display_color[2] - (30.0f / 255.0f)) < 1.0e-6f;
    }
    if (!display_synced)
	FAIL("GED visible/highlight/color setters should sync realized Obol shape display state");
    memset(&box_record, 0, sizeof(box_record));
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &box_record);
    if (!box_record.found)
	FAIL("GED shape record should still be available after display setters");
    uint32_t previous_material_revision =
	box_source_summary.materialRevision;
    const unsigned char material_color[3] = {200, 100, 50};
    if (!ged_draw_shape_ref_set_material_color(gedp, box_record.ref,
	    material_color))
	FAIL("GED material color setter should succeed");
    if (!box_source->getSummary(box_source_summary) ||
	    !box_source_summary.materialColorValid ||
	    box_source_summary.materialRevision <= previous_material_revision ||
	    fabsf(box_source_summary.materialColor[0] -
		(200.0f / 255.0f)) > 1.0e-6f ||
	    fabsf(box_source_summary.materialColor[1] -
		(100.0f / 255.0f)) > 1.0e-6f ||
	    fabsf(box_source_summary.materialColor[2] -
		(50.0f / 255.0f)) > 1.0e-6f)
	FAIL("GED material color setter should mutate the owned Obol source");
    const uint32_t refresh_material_revision = 4321;
    if (!ged_draw_shape_ref_refresh_material_color(gedp, box_record.ref,
	    gedp->dbip, refresh_material_revision))
	FAIL("GED material color refresh should succeed");
    struct ged_draw_shape_material_summary refreshed_material;
    memset(&refreshed_material, 0, sizeof(refreshed_material));
    if (!box_source->getSummary(box_source_summary) ||
	    box_source_summary.materialRevision !=
		refresh_material_revision ||
	    !ged_draw_shape_ref_material_summary(gedp, box_record.ref,
		&refreshed_material) ||
	    !refreshed_material.valid ||
	    refreshed_material.material_revision !=
		refresh_material_revision)
	FAIL("GED material color refresh should stamp the owned Obol source revision");
    box_source->tessellationAbsTol = 0.125f;
    box_source->tessellationRelTol = 0.25f;
    box_source->tessellationNormTol = 0.5f;
    struct ged_draw_shape_source_snapshot obol_source_snapshot;
    memset(&obol_source_snapshot, 0, sizeof(obol_source_snapshot));
    if (!ged_draw_shape_ref_source_snapshot(gedp, box_record.ref,
	    &obol_source_snapshot) ||
	    obol_source_snapshot.dbip != gedp->dbip ||
	    !obol_source_snapshot.fullpath ||
	    !DB_FULL_PATH_CUR_DIR(obol_source_snapshot.fullpath) ||
	    !path_equal(DB_FULL_PATH_CUR_DIR(obol_source_snapshot.fullpath)->d_namep,
		"box.s") ||
	    !obol_source_snapshot.tol ||
	    obol_source_snapshot.tol->magic != BN_TOL_MAGIC ||
	    !obol_source_snapshot.ttol ||
	    obol_source_snapshot.ttol->magic != BG_TESS_TOL_MAGIC ||
	    fabs(obol_source_snapshot.ttol->abs - 0.125) > 0.001 ||
	    fabs(obol_source_snapshot.ttol->rel - 0.25) > 0.001 ||
	    fabs(obol_source_snapshot.ttol->norm - 0.5) > 0.001)
	FAIL("GED source snapshots should read database and tessellation state from owned Obol sources");
    if (box_source->setRealizationState(SoBRLDatabaseSource::REALIZED,
	    box_source->sourceRevision.getValue(),
	    box_source->inputsRevision.getValue(),
	    SoBRLDatabaseSource::STALE_NONE) < 0)
	FAIL("owned Obol source should restore realization state after snapshot sentinel");
    ged_draw_shape_ref stale_box_ref = box_record.ref;
    uint64_t stale_box_revision = ged_draw_scene_revision(gedp);
    const char *erase_ball_for_stale_ref[2] = {"erase", "ball.s"};
    if (ged_exec_erase(gedp, 2, erase_ball_for_stale_ref) != BRLCAD_OK ||
	    ged_draw_scene_revision(gedp) <= stale_box_revision)
	FAIL("GED stale-ref sentinel should advance the draw scene revision");
    void *stale_lod_view_ctx = ged_draw_shape_ref_view_context(gedp,
	    stale_box_ref);
    if (!stale_lod_view_ctx)
	FAIL("GED stale shape-ref view context should recover cached source state");
    if (!ged_view_context_scale_set(stale_lod_view_ctx, 9.0))
	FAIL("GED stale shape-ref LoD context scale should be settable");
    struct rt_view_lod_policy stale_lod_policy = RT_VIEW_LOD_POLICY_INIT;
    stale_lod_policy.csg_enabled = 1;
    stale_lod_policy.mesh_enabled = 0;
    stale_lod_policy.bot_threshold = 91;
    stale_lod_policy.curve_scale = 6.25;
    stale_lod_policy.point_scale = 7.25;
    if (!ged_view_context_lod_policy_apply(stale_lod_view_ctx,
	    &stale_lod_policy))
	FAIL("GED stale shape-ref LoD policy should be settable");
    void *stale_lod_view_ctxs[1] = {stale_lod_view_ctx};
    if (!ged_draw_shape_ref_lod_ensure(gedp, stale_box_ref,
	    stale_lod_view_ctx, stale_lod_view_ctxs, 1))
	FAIL("GED stale shape-ref LoD ensure should recover cached Obol source runtime");
    if (!box_source->getSummary(box_realized_summary) ||
	    !(box_realized_summary.realizationRoleFlags &
		SoBRLDatabaseSource::REALIZATION_ROLE_CSG) ||
	    !box_realized_summary.realizationViewDependent ||
	    fabs(box_realized_summary.realizationViewScale - 9.0f) >
		0.001 ||
	    box_realized_summary.realizationBotThreshold != 91 ||
	    fabs(box_realized_summary.realizationCurveScale - 6.25f) >
		0.001 ||
	    fabs(box_realized_summary.realizationPointScale - 7.25f) >
		0.001)
	FAIL("GED stale shape-ref LoD ensure should update owned Obol source policy");
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("GED stale-ref sentinel should restore ball draw state");
    memset(&box_record, 0, sizeof(box_record));
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &box_record);
    if (!box_record.found)
	FAIL("GED shape record should refresh after stale-ref LoD sentinel");
    box_source = source_for_path(owned_scene, "box.s");
    if (!box_source)
	FAIL("owned Obol source should remain available after stale-ref LoD sentinel");
    group_source_state group_state = {0, GED_DRAW_GROUP_REF_NULL_INIT, NULL,
	"box.s"};
    ged_draw_foreach_group_record(gedp, group_source_state_cb,
	    &group_state);
    if (!group_state.found)
	FAIL("GED group record should be available");
    std::string group_path(group_state.path);
    if (!ged_draw_group_ref_set_visible(gedp, group_state.ref, 0))
	FAIL("GED group visible setter should succeed");
    SoGroup *box_group = owned_scene->findGroup(group_path.c_str());
    if (!box_group ||
	    !box_group->isOfType(SoBRLSceneGroup::getClassTypeId()) ||
	    static_cast<SoBRLSceneGroup *>(box_group)->visible.getValue())
	FAIL("GED group visible setter should mutate the owned Obol group");
    group_state.found = 0;
    group_state.ref = GED_DRAW_GROUP_REF_NULL;
    group_state.path = NULL;
    group_state.matchPath = group_path.c_str();
    ged_draw_foreach_group_record(gedp, group_source_state_cb,
	    &group_state);
    if (!group_state.found)
	FAIL("GED group record should still be available after visibility mutation");
    if (!ged_draw_group_ref_set_visible(gedp, group_state.ref, 1))
	FAIL("GED group visible restore should succeed");
    group_state.found = 0;
    group_state.ref = GED_DRAW_GROUP_REF_NULL;
    group_state.path = NULL;
    group_state.matchPath = group_path.c_str();
    ged_draw_foreach_group_record(gedp, group_source_state_cb,
	    &group_state);
    if (!group_state.found)
	FAIL("GED group record should still be available after visibility restore");
    if (!ged_draw_group_ref_set_mode(gedp, group_state.ref,
	    GED_DRAW_MODE_SHADED))
	FAIL("GED group mode setter should succeed");
    box_group = owned_scene->findGroup(group_path.c_str());
    if (!box_group ||
	    !box_group->isOfType(SoBRLSceneGroup::getClassTypeId()) ||
	    !static_cast<SoBRLSceneGroup *>(box_group)->
		drawIntentValid.getValue() ||
	    static_cast<SoBRLSceneGroup *>(box_group)->drawMode.getValue() !=
		BRLOBOL_LOD_DRAW_SHADED)
	FAIL("GED group mode setter should mutate the owned Obol group draw intent");
    if (!ged_draw_group_ref_set_mode(gedp, group_state.ref,
	    GED_DRAW_MODE_WIRE))
	FAIL("GED group mode restore should succeed");
    struct ged_draw_appearance_settings group_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    group_appearance.transparency = 0.45;
    group_appearance.color_override = 1;
    group_appearance.color[0] = 90;
    group_appearance.color[1] = 100;
    group_appearance.color[2] = 110;
    group_appearance.s_line_width = 6;
    if (!ged_draw_group_ref_set_appearance_settings(gedp, group_state.ref,
	    &group_appearance))
	FAIL("GED group appearance setter should succeed");
    box_group = owned_scene->findGroup(group_path.c_str());
    if (!box_group ||
	    !box_group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	FAIL("owned Obol group should remain available after appearance mutation");
    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(box_group);
    SbColor group_color = scene_group->color.getValue();
    if (scene_group->lineWidth.getValue() != 6 ||
	    fabs(scene_group->transparency.getValue() - 0.45) > 0.001 ||
	    !scene_group->colorOverride.getValue() ||
	    fabs(group_color[0] - (90.0f / 255.0f)) > 1.0e-6f ||
	    fabs(group_color[1] - (100.0f / 255.0f)) > 1.0e-6f ||
	    fabs(group_color[2] - (110.0f / 255.0f)) > 1.0e-6f)
	FAIL("GED group appearance setter should mutate the owned Obol group");
    if (owned_scene->setGroupDisplayState(group_path.c_str(),
	    scene_group->visible.getValue(),
	    scene_group->selected.getValue(),
	    scene_group->highlighted.getValue(),
	    scene_group->lineStyle.getValue(),
	    9,
	    0.23f,
	    TRUE,
	    SbColor(12.0f / 255.0f, 34.0f / 255.0f, 56.0f / 255.0f),
	    scene_group->materialColorValid.getValue(),
	    scene_group->materialColor.getValue(),
	    scene_group->materialRevision.getValue()) <= 0)
	FAIL("owned Obol group appearance sentinel update should succeed");
    struct ged_draw_appearance_settings group_appearance_readback =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    if (!ged_draw_group_ref_appearance_settings(gedp, group_state.ref,
	    &group_appearance_readback) ||
	    group_appearance_readback.s_line_width != 9 ||
	    fabs(group_appearance_readback.transparency - 0.23) > 0.001 ||
	    !group_appearance_readback.color_override ||
	    group_appearance_readback.color[0] != 12 ||
	    group_appearance_readback.color[1] != 34 ||
	    group_appearance_readback.color[2] != 56)
	FAIL("GED group appearance readback should prefer owned Obol group state");
    if (owned_scene->setGroupDisplayState(group_path.c_str(),
	    scene_group->visible.getValue(),
	    scene_group->selected.getValue(),
	    scene_group->highlighted.getValue(),
	    scene_group->lineStyle.getValue(),
	    6,
	    0.45f,
	    TRUE,
	    SbColor(90.0f / 255.0f, 100.0f / 255.0f, 110.0f / 255.0f),
	    scene_group->materialColorValid.getValue(),
	    scene_group->materialColor.getValue(),
	    scene_group->materialRevision.getValue()) <= 0)
	FAIL("owned Obol group appearance sentinel restore should succeed");
    struct ged_draw_group_record group_record;
    memset(&group_record, 0, sizeof(group_record));
    if (!ged_draw_group_record_get(gedp, group_state.ref, &group_record) ||
	    fabs(group_record.transparency - 0.45) > 0.001 ||
	    !group_record.visible ||
	    group_record.draw_mode != GED_DRAW_MODE_WIRE)
	FAIL("GED group records should read owned Obol group display state");
    const int original_group_shape_count = group_record.shape_count;
    if (original_group_shape_count <= 0)
	FAIL("GED group record should report database sources under the group");
    if (!ged_draw_obol_database_source_ensure_for_path(gedp,
	    "__obol_count_sentinel.s", gedp->dbip, GED_DRAW_MODE_WIRE,
	    1001))
	FAIL("GED Obol source count sentinel bridge should insert the source");
    if (!ged_draw_obol_database_source_move_to_group_for_path(gedp,
	    "__obol_count_sentinel.s", group_path.c_str()))
	FAIL("GED Obol source count sentinel bridge should move under the group");
    if (owned_scene->getGroupDatabaseSourceCount(group_path.c_str()) !=
	    original_group_shape_count + 1)
	FAIL("owned Obol group source count should include the sentinel");
    memset(&group_record, 0, sizeof(group_record));
    if (!ged_draw_group_record_get(gedp, group_state.ref, &group_record) ||
	    group_record.shape_count != original_group_shape_count + 1)
	FAIL("GED group record shape count should prefer owned Obol group sources");
    if (owned_scene->removeDatabaseSource("__obol_count_sentinel.s") <= 0)
	FAIL("owned Obol source count sentinel should be removable");
    struct ged_draw_obol_database_source_record source_record;
    memset(&source_record, 0, sizeof(source_record));
    if (!ged_draw_obol_database_source_record_for_path(gedp, "box.s",
	    &source_record) ||
	    !source_record.valid ||
	    source_record.realization_status !=
	    GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_CURRENT)
	FAIL("GED Obol source-record bridge should read owned source state");
    source_record.source_revision += 17;
    source_record.inputs_revision += 23;
    source_record.realization_status =
	GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_STALE;
    source_record.stale_reason = GED_DRAW_STALE_VIEW_INPUT_CHANGED;
    source_record.material_policy =
	GED_DRAW_OBOL_DATABASE_SOURCE_MATERIAL_INHERIT;
    source_record.realization_role_flags =
	SoBRLDatabaseSource::REALIZATION_ROLE_MESH;
    source_record.realization_view_dependent = 1;
    source_record.realization_view_scale = 11.0;
    source_record.realization_bot_threshold = 88;
    source_record.realization_curve_scale = 4.5;
    source_record.realization_point_scale = 5.5;
    if (!ged_draw_obol_database_source_apply_record_for_path(gedp, "box.s",
	    &source_record))
	FAIL("GED Obol source-record bridge should apply owned source state");
    if (!box_source->getSummary(box_source_summary) ||
	    !box_source_summary.stale ||
	    !(box_source_summary.staleReason &
		SoBRLDatabaseSource::STALE_INPUTS) ||
	    box_source_summary.sourceRevision !=
	    (uint32_t)source_record.source_revision ||
	    box_source_summary.inputsRevision !=
	    (uint32_t)source_record.inputs_revision ||
	    box_source_summary.materialPolicy !=
	    SoBRLDatabaseSource::MATERIAL_INHERIT ||
	    box_source_summary.realizationStatus !=
	    SoBRLDatabaseSource::UNREALIZED ||
	    box_source_summary.realizationRoleFlags !=
	    SoBRLDatabaseSource::REALIZATION_ROLE_MESH ||
	    !box_source_summary.realizationViewDependent ||
	    fabs(box_source_summary.realizationViewScale - 11.0f) > 0.001f ||
	    box_source_summary.realizationBotThreshold != 88 ||
	    fabs(box_source_summary.realizationCurveScale - 4.5f) > 0.001f ||
	    fabs(box_source_summary.realizationPointScale - 5.5f) > 0.001f)
	FAIL("GED Obol source-record bridge should mutate owned source metadata");
    (void)owned_scene->realizePending();
    if (!box_source->getSummary(box_source_summary) ||
	    box_source_summary.stale)
	FAIL("owned Obol source should be current before GED stale mutation check");
    if (ged_draw_mark_database_change(gedp, "box.s",
	    GED_DRAW_STALE_VIEW_INPUT_CHANGED) <= 0)
	FAIL("GED database-change marker should succeed");
    if (!box_source->getSummary(box_source_summary) ||
	    !box_source_summary.stale ||
	    !(box_source_summary.staleReason &
		SoBRLDatabaseSource::STALE_INPUTS))
	FAIL("GED database-change marker should mutate the owned Obol source stale state");

    const char *erase_ball[2] = {"erase", "ball.s"};
    if (ged_exec_erase(gedp, 2, erase_ball) != BRLCAD_OK)
	FAIL("owned-controller erase command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(owned_scene, "box.s") ||
	    source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror erase transactions");

    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller redraw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw transactions");
    if (owned_scene->replaceDatabaseSource("draft_move.s", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 6060) <= 0 ||
	    !source_for_path(owned_scene, "draft_move.s") ||
	    owned_scene->getDatabaseSourceCount() != 3)
	FAIL("owned Obol transaction canary source should be created");
    struct ged_draw_transaction display_txn =
	ged_draw_transaction_make_value(GED_DRAW_TXN_VISIBILITY,
		"box.s", 0.0);
    struct ged_draw_transaction_result display_result;
    ged_draw_transaction_result_init(&display_result);
    if (ged_draw_apply_transaction(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED visibility transaction should succeed");
    ged_draw_transaction_result_free(&display_result);
    box_source = source_for_path(owned_scene, "box.s");
    if (!box_source ||
	    !box_source->getSummary(box_source_summary) ||
	    box_source_summary.visible ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED visibility transaction should update owned Obol state without full-scene sync");
    display_txn = ged_draw_transaction_make_value(GED_DRAW_TXN_VISIBILITY,
	    "box.s", 1.0);
    ged_draw_transaction_result_init(&display_result);
    if (ged_draw_apply_transaction(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED visibility restore transaction should succeed");
    ged_draw_transaction_result_free(&display_result);
    display_txn = ged_draw_transaction_make_value(GED_DRAW_TXN_TRANSPARENCY,
	    "box.s", 0.125);
    ged_draw_transaction_result_init(&display_result);
    if (ged_draw_apply_transaction(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED transparency transaction should succeed");
    ged_draw_transaction_result_free(&display_result);
    if (!box_source->getSummary(box_source_summary) ||
	    fabs(box_source_summary.transparency - 0.125f) > 0.001f ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED transparency transaction should preserve Obol-only sources");
    display_txn = ged_draw_transaction_make_value(GED_DRAW_TXN_HIGHLIGHT,
	    "box.s", 1.0);
    ged_draw_transaction_result_init(&display_result);
    if (ged_draw_apply_transaction(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED highlight transaction should succeed");
    ged_draw_transaction_result_free(&display_result);
    if (!box_source->getSummary(box_source_summary) ||
	    !box_source_summary.highlighted ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED highlight transaction should preserve Obol-only sources");
    display_txn = ged_draw_transaction_make_value(GED_DRAW_TXN_HIGHLIGHT,
	    "box.s", 0.0);
    ged_draw_transaction_result_init(&display_result);
    if (ged_draw_apply_transaction(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED highlight restore transaction should succeed");
    ged_draw_transaction_result_free(&display_result);
    display_txn = ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE,
	    "box.s");
    display_txn.stale_reason = GED_DRAW_STALE_SETTINGS_CHANGED;
    ged_draw_transaction_result_init(&display_result);
    if (ged_draw_apply_transaction(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED stale-source transaction should succeed");
    ged_draw_transaction_result_free(&display_result);
    if (!box_source->getSummary(box_source_summary) ||
	    !box_source_summary.stale ||
	    !(box_source_summary.staleReason &
		SoBRLDatabaseSource::STALE_DRAW) ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED stale-source transaction should target owned Obol state without full-scene sync");
    display_txn = ged_draw_transaction_make(GED_DRAW_TXN_MATERIAL_CHANGED,
	    NULL);
    ged_draw_transaction_result_init(&display_result);
    if (ged_draw_apply_transaction(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED material-changed transaction should succeed");
    ged_draw_transaction_result_free(&display_result);
    if (!source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED material-changed transaction should preserve Obol-only sources");
    display_txn = ged_draw_transaction_make(GED_DRAW_TXN_REFRESH_MATERIAL_COLORS,
	    NULL);
    ged_draw_transaction_result_init(&display_result);
    if (ged_draw_apply_transaction(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED material-refresh transaction should succeed");
    ged_draw_transaction_result_free(&display_result);
    if (!source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED material-refresh transaction should preserve Obol-only sources");
    display_txn = ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    ged_draw_transaction_result_init(&display_result);
    if (ged_draw_apply_transaction(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED redraw transaction should succeed");
    ged_draw_transaction_result_free(&display_result);
    if (!source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s") ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED redraw transaction should refresh retained sources without clearing Obol-only sources");
    display_txn = ged_draw_transaction_make(GED_DRAW_TXN_ERASE_PREFIX,
	    "box.s");
    ged_draw_transaction_result_init(&display_result);
    if (ged_draw_apply_transaction(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED erase-prefix transaction should succeed");
    ged_draw_transaction_result_free(&display_result);
    if (source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s") ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED erase-prefix transaction should remove only matching owned Obol sources");
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK ||
	    !source_for_path(owned_scene, "box.s"))
	FAIL("GED erase-prefix canary redraw should restore box source");
    if (!source_for_path(owned_scene, "draft_move.s") ||
	    owned_scene->removeDatabaseSource("draft_move.s") <= 0 ||
	    owned_scene->getDatabaseSourceCount() != 2)
	FAIL("GED display/material transaction canary cleanup should restore baseline");
    record_source_state ball_record = {0, GED_DRAW_SHAPE_REF_NULL_INIT,
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "ball.s", 0, 0, 0, 0, 0.0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &ball_record);
    if (!ball_record.found)
	FAIL("GED ball shape record should be available before direct release");
    if (!ged_draw_shape_ref_release(gedp, ball_record.ref))
	FAIL("GED shape direct release should succeed");
    if (owned_scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(owned_scene, "box.s") ||
	    source_for_path(owned_scene, "ball.s"))
	FAIL("GED shape direct release should remove the owned Obol source");
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller direct-release redraw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after direct release");
    if (!ged_draw_source_erase_path_in_active_scope(gedp, "ball.s",
	    ged_draw_active_view_ctx(gedp), -1))
	FAIL("GED direct active-scope erase path should succeed");
    if (owned_scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(owned_scene, "box.s") ||
	    source_for_path(owned_scene, "ball.s"))
	FAIL("GED direct active-scope group erase should remove the owned Obol group subtree");
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller active-scope erase redraw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after direct active-scope erase");
    if (!ged_draw_scene_root_erase_path(gedp, "ball.s"))
	FAIL("GED direct root erase path should succeed");
    if (owned_scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(owned_scene, "box.s") ||
	    source_for_path(owned_scene, "ball.s"))
	FAIL("GED direct root erase should remove the owned Obol source");
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller root-erase redraw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after direct root erase");
    if (!ged_draw_scene_root_erase_groups_by_name(gedp, "ball.s"))
	FAIL("GED direct root group-name erase should succeed");
    if (owned_scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(owned_scene, "box.s") ||
	    source_for_path(owned_scene, "ball.s"))
	FAIL("GED direct root group-name erase should remove the owned Obol group");
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller group-name erase redraw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after direct group-name erase");
    const char *nested_leaf_source_path =
	"nested_parent.c/nested_child.c/nested_leaf.s";
    const char *nested_sibling_source_path =
	"nested_parent.c/nested_sibling.s";
    const char *draw_nested_leaf[2] = {"draw", nested_leaf_source_path};
    const char *draw_nested_sibling[2] = {
	"draw", nested_sibling_source_path
    };
    if (ged_exec_draw(gedp, 2, draw_nested_leaf) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 2, draw_nested_sibling) != BRLCAD_OK)
	FAIL("GED nested path draws should succeed");
    SoBRLDatabaseSource *nested_leaf_source =
	source_for_path(owned_scene, nested_leaf_source_path);
    SoBRLDatabaseSource *nested_sibling_source =
	source_for_path(owned_scene, nested_sibling_source_path);
    if (!nested_leaf_source || !nested_sibling_source ||
	    owned_scene->getDatabaseSourceCount() != 4)
	FAIL("GED nested path draws should create owned Obol child sources");
    BRLObolDatabaseSourceSummary nested_sibling_initial_summary;
    if (!nested_sibling_source->getSummary(nested_sibling_initial_summary))
	FAIL("GED nested sibling source should expose a state summary");
    if (owned_scene->setDatabaseSourceState(nested_sibling_source_path,
	    TRUE,
	    nested_sibling_initial_summary.sourceRevision,
	    nested_sibling_initial_summary.inputsRevision,
	    nested_sibling_initial_summary.visible,
	    nested_sibling_initial_summary.highlighted,
	    nested_sibling_initial_summary.lineStyle,
	    23,
	    nested_sibling_initial_summary.transparency,
	    nested_sibling_initial_summary.colorOverride,
	    nested_sibling_initial_summary.color,
	    nested_sibling_initial_summary.materialColorValid,
	    nested_sibling_initial_summary.materialColor,
	    nested_sibling_initial_summary.materialRevision) <= 0)
	FAIL("owned Obol nested sibling state sentinel update should succeed");
    struct ged_draw_transaction stale_nested_leaf =
	ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED,
		"nested_leaf.s");
    stale_nested_leaf.redraw = 0;
    struct ged_draw_transaction_result nested_result;
    ged_draw_transaction_result_init(&nested_result);
    if (ged_draw_apply_transaction(gedp, &stale_nested_leaf,
	    &nested_result) <= 0)
	FAIL("GED nested leaf stale transaction should succeed");
    ged_draw_transaction_result_free(&nested_result);
    BRLObolDatabaseSourceSummary nested_leaf_summary;
    BRLObolDatabaseSourceSummary nested_sibling_summary;
    nested_leaf_source = source_for_path(owned_scene,
	    nested_leaf_source_path);
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (!nested_leaf_source ||
	    !nested_leaf_source->getSummary(nested_leaf_summary))
	FAIL("GED nested leaf stale transaction should preserve the owned Obol leaf source");
    if (!nested_leaf_summary.stale)
	FAIL("GED nested leaf stale transaction should mark the owned Obol leaf source stale");
    if (!(nested_leaf_summary.staleReason &
	    SoBRLDatabaseSource::STALE_SOURCE))
	FAIL("GED nested leaf stale transaction should use source-stale metadata");
    if (!nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary))
	FAIL("GED nested leaf stale transaction should preserve the sibling source");
    if (nested_sibling_summary.stale)
	FAIL("GED nested leaf stale transaction should not stale the sibling source");
    if (nested_sibling_summary.lineWidth != 23)
	FAIL("GED nested leaf stale transaction should preserve sibling source state");
    struct ged_draw_transaction redraw_nested_leaf =
	ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED,
		"nested_leaf.s");
    redraw_nested_leaf.redraw = 1;
    ged_draw_transaction_result_init(&nested_result);
    if (ged_draw_apply_transaction(gedp, &redraw_nested_leaf,
	    &nested_result) <= 0)
	FAIL("GED nested leaf redraw transaction should succeed");
    ged_draw_transaction_result_free(&nested_result);
    nested_leaf_source = source_for_path(owned_scene,
	    nested_leaf_source_path);
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (!nested_leaf_source ||
	    !nested_leaf_source->getSummary(nested_leaf_summary) ||
	    nested_leaf_summary.stale ||
	    !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED nested leaf redraw transaction should preserve unrelated owned Obol source state");
    if (!ged_draw_source_erase_component_name_in_active_scope(gedp,
	    "nested_leaf.s", ged_draw_active_view_ctx(gedp), -1, 1))
	FAIL("GED scoped component erase should succeed through active-scope source adapter");
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (source_for_path(owned_scene, nested_leaf_source_path) ||
	    !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED scoped component erase should remove matching owned Obol sources only");
    if (ged_exec_draw(gedp, 2, draw_nested_leaf) != BRLCAD_OK)
	FAIL("GED scoped component erase redraw restore should succeed");
    nested_leaf_source = source_for_path(owned_scene,
	    nested_leaf_source_path);
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (!nested_leaf_source || !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED scoped component erase redraw restore should preserve sibling owned Obol source state");
    struct ged_draw_transaction remove_nested_leaf_ref =
	ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED,
		"nested_leaf.s");
    ged_draw_transaction_result_init(&nested_result);
    if (ged_draw_apply_transaction(gedp, &remove_nested_leaf_ref,
	    &nested_result) <= 0)
	FAIL("GED nested leaf reference removal transaction should succeed");
    ged_draw_transaction_result_free(&nested_result);
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (source_for_path(owned_scene, nested_leaf_source_path) ||
	    !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED nested leaf reference removal should remove only non-root owned Obol sources");
    struct ged_draw_transaction remove_nested_sibling =
	ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED,
		"nested_sibling.s");
    remove_nested_sibling.removed = 1;
    ged_draw_transaction_result_init(&nested_result);
    if (ged_draw_apply_transaction(gedp, &remove_nested_sibling,
	    &nested_result) <= 0)
	FAIL("GED nested sibling source removal transaction should succeed");
    ged_draw_transaction_result_free(&nested_result);
    if (source_for_path(owned_scene, nested_leaf_source_path) ||
	    source_for_path(owned_scene, nested_sibling_source_path) ||
	    owned_scene->getDatabaseSourceCount() != 2)
	FAIL("GED source removal transaction should remove matching owned Obol component sources");
    struct db_full_path nested_child_path;
    db_full_path_init(&nested_child_path);
    if (db_string_to_path(&nested_child_path, gedp->dbip,
	    "nested_parent.c/nested_child.c") != 0)
	FAIL("GED nested erase sentinel path should resolve");
    ged_draw_group_ref nested_child_group =
	ged_draw_group_ref_lookup_or_create(gedp, &nested_child_path);
    db_free_full_path(&nested_child_path);
    if (ged_draw_group_ref_is_null(nested_child_group))
	FAIL("GED nested erase sentinel group should be created");
    if (!owned_scene->findGroup("nested_parent.c/nested_child.c"))
	FAIL("owned Obol nested erase sentinel child group should exist after creation");
    if (!ged_draw_source_erase_path_in_active_scope(gedp,
	    "nested_parent.c/nested_child.c", ged_draw_active_view_ctx(gedp),
	    -1))
	FAIL("GED direct active-scope nested group path erase should succeed");
    if (!owned_scene->findGroup("nested_parent.c") ||
	    owned_scene->findGroup("nested_parent.c/nested_child.c"))
	FAIL("GED direct active-scope nested group path erase should remove the owned Obol nested group");
    if (!ged_draw_source_erase_path_in_active_scope(gedp, "nested_parent.c",
	    ged_draw_active_view_ctx(gedp), -1) ||
	    owned_scene->findGroup("nested_parent.c") ||
	    owned_scene->getDatabaseSourceCount() != 2)
	FAIL("GED direct active-scope nested group path cleanup should restore owned Obol baseline");
    std::string clear_child_path = group_path + "/__obol_clear_child";
    if (!owned_scene->ensureGroup(clear_child_path.c_str()) ||
	    !owned_scene->findGroup(clear_child_path.c_str()))
	FAIL("owned Obol direct clear sentinel child group should be created");
    ged_draw_clear(gedp);
    if (owned_scene->getDatabaseSourceCount() != 0 ||
	    owned_scene->findGroup(clear_child_path.c_str()) ||
	    owned_scene->findGroup(group_path.c_str()))
	FAIL("GED direct draw clear should remove owned Obol group/source subtrees");
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller direct-clear redraw commands should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after direct clear");
    record_source_state draft_box_record = {0, GED_DRAW_SHAPE_REF_NULL_INIT,
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "box.s", 0, 0, 0, 0, 0.0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &draft_box_record);
    if (!draft_box_record.found)
	FAIL("GED Obol draft append sentinel target group should be available");
    struct ged_draw_group_record draft_group_record;
    memset(&draft_group_record, 0, sizeof(draft_group_record));
    if (!ged_draw_group_record_get(gedp, draft_box_record.group,
	    &draft_group_record) ||
	    !draft_group_record.path)
	FAIL("GED Obol draft append sentinel target group should have a path");
    int draft_group_source_count =
	owned_scene->getGroupDatabaseSourceCount(draft_group_record.path);
    if (draft_group_source_count < 0)
	FAIL("owned Obol draft append sentinel target group should exist");
    if (!ged_draw_obol_database_source_ensure_for_path(gedp,
	    "draft_move.s", gedp->dbip, GED_DRAW_MODE_WIRE, 5150))
	FAIL("GED Obol draft append sentinel source should be ensured");
    struct db_full_path draft_move_path;
    db_full_path_init(&draft_move_path);
    if (db_string_to_path(&draft_move_path, gedp->dbip,
	    "draft_move.s") != 0)
	FAIL("GED Obol draft append sentinel path should resolve");
    ged_draw_shape_draft *draft_move =
	ged_draw_shape_draft_create_context(gedp,
		ged_draw_active_view_ctx(gedp), 1);
    if (!draft_move)
	FAIL("GED Obol draft append sentinel draft should be created");
    struct bn_tol draft_tol;
    struct bg_tess_tol draft_ttol;
    BN_TOL_INIT_SET_TOL(&draft_tol);
    BG_TESS_TOL_INIT_SET_TOL(&draft_ttol);
    if (!ged_draw_shape_draft_apply_path_source_state(draft_move,
	    gedp->dbip, &draft_move_path, &draft_tol, &draft_ttol,
	    0, NULL, "draft_move.s"))
	FAIL("GED Obol draft append sentinel source state should apply");
    point_t draft_points[2] = {
	{0.0, 5.0, 0.0},
	{1.0, 5.0, 0.0}
    };
    int draft_cmds[2] = {
	GED_DRAW_VIEW_LINE_MOVE,
	GED_DRAW_VIEW_LINE_DRAW
    };
    if (!ged_draw_shape_draft_publish_line_set(draft_move,
	    (const point_t *)draft_points, draft_cmds, 2))
	FAIL("GED Obol draft append sentinel line set should publish");
    ged_draw_shape_ref draft_move_ref =
	ged_draw_shape_draft_commit_to_group(draft_move,
		draft_box_record.group);
    db_free_full_path(&draft_move_path);
    if (ged_draw_shape_ref_is_null(draft_move_ref))
	FAIL("GED Obol draft append sentinel should commit to the group");
    if (owned_scene->getGroupDatabaseSourceCount(draft_group_record.path) !=
	    draft_group_source_count + 1 ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED source-owner append should move the owned Obol source into the target group");
    if (!ged_draw_shape_ref_release(gedp, draft_move_ref) ||
	    source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED Obol draft append sentinel cleanup should remove the owned source");

    ged_draw_obol_scene_controller_detach(gedp);
    if (ged_draw_obol_scene_controller(gedp) ||
	    ged_draw_obol_scene_controller_owned(gedp))
	FAIL("owned Obol scene controller detach should clear ownership state");

    BRLObolViewController view_controller(root);
    if (!ged_draw_obol_controller_attach(gedp, &view_controller, 1))
	FAIL("GED Obol view-controller attachment should succeed");
    if (ged_draw_obol_controller(gedp) != &view_controller ||
	    ged_draw_obol_scene_controller(gedp) !=
		view_controller.getSceneController())
	FAIL("GED Obol view-controller attachment should expose its scene controller");
    SoBRLSceneController *view_scene = view_controller.getSceneController();
    if (view_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(view_scene, "box.s") ||
	    !source_for_path(view_scene, "ball.s"))
	FAIL("view-wrapper reattach full sync should rebuild current GED draw state");

    const char *erase_box[2] = {"erase", "box.s"};
    if (ged_exec_erase(gedp, 2, erase_box) != BRLCAD_OK)
	FAIL("real GED erase command should succeed");
    if (view_scene->getDatabaseSourceCount() != 1 ||
	    source_for_path(view_scene, "box.s") ||
	    !source_for_path(view_scene, "ball.s"))
	FAIL("attached Obol scene controller should mirror GED erase command");

    const char *zap_cmd[1] = {"zap"};
    if (ged_exec_zap(gedp, 1, zap_cmd) != BRLCAD_OK)
	FAIL("real GED zap command should succeed");
    if (view_scene->getDatabaseSourceCount() != 0)
	FAIL("attached Obol scene controller should mirror GED clear command");

    ged_draw_obol_controller_detach(gedp);
    ged_close(gedp);
    root->unref();
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
