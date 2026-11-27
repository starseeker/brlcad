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

#include "ged/display_obol_private.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BDrawCache.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BExportAction.h"
#include "BObol/BInit.h"
#include "BObol/BLodRealization.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshShape.h"
#include "BObol/BPickDetail.h"
#include "BObol/BSceneController.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BVListShape.h"
#include "BObol/BViewController.h"
#include "BObol/BViewLod.h"
#include "BObol/BViewStore.h"
#include "bg/line_layer.h"
#include "bg/plot3.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/process.h"
#include "bu/str.h"
#include "ged.h"
#include "ged/commands.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "../../ged_scene_record_api_private.h"
#include "ged/selection.h"
#include "ged/view.h"
#include "opennurbs_sphere.h"
#include "rt/db_internal.h"
#include "rt/edit.h"
#include "rt/view.h"
#include "view_test_util.h"
#include "wdb.h"

#include "../../ged_private.h"
#include "../../ged_bobol_private.hpp"
#include "../../ged_draw_private.h"
#include "../../../libBObol/tests/transaction_fault_test_private.h"

#include <algorithm>
#include <exception>
#include <array>
#include <Inventor/sensors/SoFieldSensor.h>
#include <Inventor/sensors/SoNodeSensor.h>

#include <Inventor/SbMatrix.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/tools/SbModernUtils.h>

#include <Obol/cad/SoCADAssembly.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <chrono>
#include <limits>
#include <set>
#include <string>
#include <thread>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

struct occurrence_visit_state {
    size_t count;
    int valid;
};

static BObolPresentationTimingContext
test_capacity_cad_timing(void)
{
    return BObolPresentationTimingContext(
	BObolLodCapacityRelevance::RELEVANT,
	BObolLodPlanningRelevance::RELEVANT,
	BObolCadPresentationExecution::EXECUTED,
	BOBOL_CAD_PREPARATION_NONE,
	BObolCadPresentationCompleteness::EXACT);
}

static int
occurrence_visit_cb(const struct ged_scene_occurrence_info *occurrence,
	void *client_data)
{
    struct occurrence_visit_state *state =
	static_cast<struct occurrence_visit_state *>(client_data);
    if (!state || !occurrence)
	return 0;
    state->count++;
    if (ged_scene_occurrence_ref_is_null(occurrence->ref) ||
	!occurrence->path || !occurrence->path[0])
	state->valid = 0;
    return 1;
}

static struct ged_view_feature_batch *
test_feature_batch(struct ged_view_context *view_ctx, int local,
	int overlay_class, int lifecycle, int order)
{
    struct ged_view_feature_batch_desc desc = ged_view_feature_batch_desc_default();
    desc.owner_id = "obol-sync-test";
    desc.owner_role = "test-feature";
    desc.local = local;
    desc.overlay_class = overlay_class;
    desc.lifecycle = lifecycle;
    desc.overlay_order = order;
    return ged_view_feature_batch_begin(view_ctx, &desc);
}

static int
test_line_replace(struct ged_view_context *view_ctx, const char *name,
	int local, const point_t *points, const int *commands, size_t count,
	const struct ged_view_feature_style *style, int overlay_class,
	int lifecycle, int order)
{
    struct ged_view_feature_batch *batch = test_feature_batch(view_ctx,
	local, overlay_class, lifecycle, order);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_line_set_replace(batch, name, points,
	    commands, count, style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
test_labels_replace(struct ged_view_context *view_ctx, const char *name,
	int local, const struct ged_view_feature_label *labels, size_t count,
	const struct ged_view_feature_style *style)
{
    struct ged_view_feature_batch *batch = test_feature_batch(view_ctx,
	local, GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION,
	GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT,
	GED_VIEW_FEATURE_OVERLAY_ORDER_MODEL);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_labels_replace(batch, name, labels, count,
	    style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
test_arrows_replace(struct ged_view_context *view_ctx, const char *name,
	const point_t *points, size_t count,
	const struct ged_view_feature_style *style)
{
    struct ged_view_feature_batch *batch = test_feature_batch(view_ctx, 1,
	GED_VIEW_FEATURE_OVERLAY_CLASS_TCL_OVERLAY,
	GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND,
	GED_VIEW_FEATURE_OVERLAY_ORDER_POST_TRANSPARENT);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_arrow_replace(batch, name, points, count,
	    style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
test_axes_replace(struct ged_view_context *view_ctx, const char *name,
	int local, const point_t *centers, size_t count, fastf_t half_size,
	const struct ged_view_feature_style *style, int overlay_class,
	int lifecycle, int order)
{
    struct ged_view_feature_batch *batch = test_feature_batch(view_ctx,
	local, overlay_class, lifecycle, order);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_axes_replace(batch, name, centers, count,
	    half_size, style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
test_mesh_replace(struct ged_view_context *view_ctx, const char *name,
	const point_t *points, size_t point_count, const int *indices,
	size_t index_count, const struct ged_view_feature_style *style)
{
    struct ged_view_feature_batch *batch = test_feature_batch(view_ctx, 0,
	GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION,
	GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT,
	GED_VIEW_FEATURE_OVERLAY_ORDER_MODEL);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_indexed_face_set_replace(batch, name, points,
	    point_count, NULL, 0, indices, index_count, style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
test_line_builder_replace(struct ged_view_context *view_ctx, const char *name,
	const struct bg_line_layer_builder *builder)
{
    struct ged_view_feature_batch *batch = test_feature_batch(view_ctx, 0,
	GED_VIEW_FEATURE_OVERLAY_CLASS_DIAGNOSTIC,
	GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND,
	GED_VIEW_FEATURE_OVERLAY_ORDER_POST_TRANSPARENT);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_line_layer_builder_replace(batch, name,
	    builder, NULL)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
exercise_ged_value_handle_lifetimes(struct ged *gedp)
{
    struct ged_view_context *view_ctx = ged_view_context_create();
    if (!view_ctx || !ged_view_context_host_attach(gedp, view_ctx))
	FAIL("value-handle lifetime test should create a hosted view");

    bobol_display_endpoint_t *first_endpoint =
	bobol_display_endpoint_create(NULL, 0);
    if (!first_endpoint ||
	!ged_view_context_obol_endpoint_set(view_ctx, first_endpoint, 1)) {
	if (first_endpoint)
	    bobol_display_endpoint_destroy(first_endpoint);
	ged_view_context_free(view_ctx);
	FAIL("value-handle lifetime test should attach its first endpoint");
    }

    ged_view_edit_ref first_feature =
	ged_view_edit_overlay_ensure(view_ctx,
	    "handle-lifetime::feature", "handle-lifetime::source.s");
    point_t origin = {0.0, 0.0, 0.0};
    ged_view_polygon_ref first_polygon =
	ged_view_polygon_create(view_ctx,
	    "handle-lifetime::polygon", 1, GED_VIEW_POLYGON_SQUARE,
	    origin);
    if (ged_view_edit_ref_is_null(first_feature) ||
	ged_view_polygon_ref_is_null(first_polygon)) {
	ged_view_context_free(view_ctx);
	FAIL("value-handle lifetime test should create live references");
    }

    struct ged_view_context *foreign_view = ged_view_context_create();
    struct ged_view_polygon_record foreign_record;
    if (!foreign_view ||
	!ged_view_context_host_attach(gedp, foreign_view) ||
	ged_view_edit_visible_set(foreign_view, first_feature, 1) ||
	ged_view_polygon_record_get(foreign_view, first_polygon,
	    &foreign_record)) {
	if (foreign_view)
	    ged_view_context_free(foreign_view);
	ged_view_context_free(view_ctx);
	FAIL("a polygon reference must fail against a foreign view owner");
    }
    ged_view_context_free(foreign_view);

    bobol_display_endpoint_t *second_endpoint =
	bobol_display_endpoint_create(NULL, 0);
    if (!second_endpoint ||
	!ged_view_context_obol_endpoint_set(view_ctx, second_endpoint, 1)) {
	if (second_endpoint)
	    bobol_display_endpoint_destroy(second_endpoint);
	ged_view_context_free(view_ctx);
	FAIL("value-handle lifetime test should replace its endpoint");
    }

    struct ged_view_polygon_record polygon_record;
    if (ged_view_edit_visible_set(view_ctx, first_feature, 1) ||
	ged_view_polygon_record_get(view_ctx, first_polygon, &polygon_record)) {
	ged_view_context_free(view_ctx);
	FAIL("references from a replaced controller must be stale");
    }

    ged_view_edit_ref second_feature =
	ged_view_edit_overlay_ensure(view_ctx,
	    "handle-lifetime::feature", "handle-lifetime::source.s");
    ged_view_polygon_ref second_polygon =
	ged_view_polygon_create(view_ctx,
	    "handle-lifetime::polygon", 1, GED_VIEW_POLYGON_SQUARE,
	    origin);
    if (ged_view_edit_ref_is_null(second_feature) ||
	ged_view_polygon_ref_is_null(second_polygon) ||
	first_feature.generation == second_feature.generation ||
	first_polygon.generation == second_polygon.generation) {
	ged_view_context_free(view_ctx);
	FAIL("replacement stores must issue a new reference generation");
    }

    if (!ged_view_feature_remove(view_ctx,
	    "handle-lifetime::feature") ||
	ged_view_edit_visible_set(view_ctx, second_feature, 1) ||
	!ged_view_polygon_remove(view_ctx, second_polygon) ||
	ged_view_polygon_record_get(view_ctx, second_polygon, &polygon_record)) {
	ged_view_context_free(view_ctx);
	FAIL("removed feature and polygon references must be stale");
    }

    ged_view_edit_ref recreated_feature =
	ged_view_edit_overlay_ensure(view_ctx,
	    "handle-lifetime::feature", "handle-lifetime::source.s");
    ged_view_polygon_ref recreated_polygon =
	ged_view_polygon_create(view_ctx,
	    "handle-lifetime::polygon", 1, GED_VIEW_POLYGON_SQUARE,
	    origin);
    if (ged_view_edit_ref_is_null(recreated_feature) ||
	ged_view_polygon_ref_is_null(recreated_polygon) ||
	recreated_feature.id == second_feature.id ||
	recreated_polygon.id == second_polygon.id ||
	!ged_view_edit_visible_set(view_ctx, recreated_feature, 1) ||
	!ged_view_polygon_record_get(view_ctx, recreated_polygon,
	    &polygon_record)) {
	ged_view_context_free(view_ctx);
	FAIL("recreated objects must have distinct, live references");
    }

    ged_view_context_free(view_ctx);

    /* Keep this loop in the regular test so ASan/LSan configurations exercise
     * registry, endpoint, store, and reference teardown repeatedly. */
    for (int i = 0; i < 64; i++) {
	struct ged_view_context *cycle_view = ged_view_context_create();
	bobol_display_endpoint_t *cycle_endpoint =
	    bobol_display_endpoint_create(NULL, 0);
	if (!cycle_view || !cycle_endpoint ||
	    !ged_view_context_host_attach(gedp, cycle_view) ||
	    !ged_view_context_obol_endpoint_set(cycle_view,
		cycle_endpoint, 1)) {
	    if (cycle_endpoint && (!cycle_view ||
		ged_view_context_obol_endpoint_get(cycle_view) !=
		    cycle_endpoint))
		bobol_display_endpoint_destroy(cycle_endpoint);
	    if (cycle_view)
		ged_view_context_free(cycle_view);
	    FAIL("repeated handle lifetime cycle should create its endpoint");
	}
	ged_view_edit_ref cycle_ref =
	    ged_view_edit_overlay_ensure(cycle_view,
		"handle-lifetime::cycle", NULL);
	if (ged_view_edit_ref_is_null(cycle_ref)) {
	    ged_view_context_free(cycle_view);
	    FAIL("repeated handle lifetime cycle should issue a reference");
	}
	ged_view_context_free(cycle_view);
    }

    return 0;
}

static int
make_obol_sync_brep_sphere(struct rt_wdb *wdbp, const char *name)
{
    ON_Sphere sphere(ON_3dPoint(4.0, 4.0, 0.0), 1.5);
    ON_Brep *brep = ON_BrepSphere(sphere);
    int ret = (brep && mk_brep(wdbp, name, (void *)brep) == 0);
    delete brep;
    return ret;
}

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
    point_t reuse_min = {-1.0, -1.0, -1.0};
    point_t reuse_max = { 1.0,  1.0,  1.0};
    point_t duplicate_min = {70.0, -1.0, -1.0};
    point_t duplicate_max = {72.0,  1.0,  1.0};
    const int mesh_grid = 12;
    const int mesh_vertex_count = (mesh_grid + 1) * (mesh_grid + 1);
    const int mesh_face_count = mesh_grid * mesh_grid * 2;
    std::vector<fastf_t> mesh_owner_vertices(mesh_vertex_count * 3, 0.0);
    std::vector<int> mesh_owner_faces(mesh_face_count * 3, 0);
    for (int y = 0; y <= mesh_grid; y++) {
	for (int x = 0; x <= mesh_grid; x++) {
	    const int vertex = y * (mesh_grid + 1) + x;
	    mesh_owner_vertices[3 * vertex + X] = (fastf_t)x;
	    mesh_owner_vertices[3 * vertex + Y] = (fastf_t)y;
	    mesh_owner_vertices[3 * vertex + Z] =
		(fastf_t)((x + y) % 3) * 0.05;
	}
    }
    for (int y = 0; y < mesh_grid; y++) {
	for (int x = 0; x < mesh_grid; x++) {
	    const int cell = y * mesh_grid + x;
	    const int v0 = y * (mesh_grid + 1) + x;
	    const int v1 = v0 + 1;
	    const int v2 = v0 + mesh_grid + 1;
	    const int v3 = v2 + 1;
	    mesh_owner_faces[6 * cell + 0] = v0;
	    mesh_owner_faces[6 * cell + 1] = v1;
	    mesh_owner_faces[6 * cell + 2] = v3;
	    mesh_owner_faces[6 * cell + 3] = v0;
	    mesh_owner_faces[6 * cell + 4] = v3;
	    mesh_owner_faces[6 * cell + 5] = v2;
	}
    }
    char binunif_payload[4] = {1, 2, 3, 4};

    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0 &&
	mk_sph(wdbp, "ball.s", center, 1.0) == 0 &&
	mk_sph(wdbp, "group_only.s", group_center, 1.0) == 0 &&
	mk_sph(wdbp, "draft_move.s", draft_center, 1.0) == 0 &&
	mk_sph(wdbp, "nested_leaf.s", nested_leaf_center, 1.0) == 0 &&
	mk_sph(wdbp, "nested_sibling.s", nested_sibling_center, 1.0) == 0 &&
	mk_sph(wdbp, "rename_source.s", rename_center, 1.0) == 0 &&
	mk_rpp(wdbp, "reuse_leaf.s", reuse_min, reuse_max) == 0 &&
	mk_rpp(wdbp, "dup_leaf.s", duplicate_min, duplicate_max) == 0 &&
	mk_bot(wdbp, "mesh_owner.bot", RT_BOT_SURFACE, RT_BOT_CCW, 0,
		mesh_vertex_count, mesh_face_count, mesh_owner_vertices.data(),
		mesh_owner_faces.data(), NULL, NULL) == 0 &&
	make_obol_sync_brep_sphere(wdbp, "brep_owner.brep") &&
	mk_binunif(wdbp, "payload.binunif", binunif_payload,
		WDB_BINUNIF_INT8, 4) == 0 &&
	mk_submodel(wdbp, "submodel_owner.s", NULL, "box.s", 0) == 0 &&
	mk_submodel(wdbp, "submodel_temp_owner.s", NULL,
		"nested_leaf.s", 0) == 0;
    if (ret) {
	struct wmember child_wm;
	struct wmember renamed_child_wm;
	struct wmember parent_wm;
	BU_LIST_INIT(&child_wm.l);
	BU_LIST_INIT(&renamed_child_wm.l);
	BU_LIST_INIT(&parent_wm.l);
	ret = mk_addmember("nested_leaf.s", &child_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_comb(wdbp, "nested_child.c", &child_wm.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) == 0 &&
	    mk_addmember("nested_leaf.s", &renamed_child_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_comb(wdbp, "nested_child_renamed.c", &renamed_child_wm.l, 0,
		NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0) == 0 &&
	    mk_addmember("nested_child.c", &parent_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_addmember("nested_child_renamed.c", &parent_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_addmember("nested_sibling.s", &parent_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_comb(wdbp, "nested_parent.c", &parent_wm.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) == 0;
    }
    if (ret) {
	struct wmember shared_wm;
	struct wmember inst_a_wm;
	struct wmember inst_b_wm;
	struct wmember root_wm;
	mat_t shared_leaf_mat;
	mat_t inst_a_mat;
	mat_t inst_b_mat;
	MAT_IDN(shared_leaf_mat);
	MAT_IDN(inst_a_mat);
	MAT_IDN(inst_b_mat);
	shared_leaf_mat[MDY] = 3.0;
	inst_a_mat[MDX] = -12.0;
	inst_b_mat[MDX] = 18.0;
	BU_LIST_INIT(&shared_wm.l);
	BU_LIST_INIT(&inst_a_wm.l);
	BU_LIST_INIT(&inst_b_wm.l);
	BU_LIST_INIT(&root_wm.l);
	ret = mk_addmember("reuse_leaf.s", &shared_wm.l, shared_leaf_mat,
		WMOP_UNION) != NULL &&
	    mk_comb(wdbp, "reuse_shared.c", &shared_wm.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) == 0 &&
	    mk_addmember("reuse_shared.c", &inst_a_wm.l, inst_a_mat,
		WMOP_UNION) != NULL &&
	    mk_comb(wdbp, "reuse_inst_a.c", &inst_a_wm.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) == 0 &&
	    mk_addmember("reuse_shared.c", &inst_b_wm.l, inst_b_mat,
		WMOP_UNION) != NULL &&
	    mk_comb(wdbp, "reuse_inst_b.c", &inst_b_wm.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) == 0 &&
	    mk_addmember("reuse_inst_a.c", &root_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_addmember("reuse_inst_b.c", &root_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_comb(wdbp, "reuse_root.c", &root_wm.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) == 0;
    }
    if (ret) {
	struct wmember duplicate_wm;
	BU_LIST_INIT(&duplicate_wm.l);
	ret = mk_addmember("dup_leaf.s", &duplicate_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_addmember("dup_leaf.s", &duplicate_wm.l, NULL,
		WMOP_UNION) != NULL &&
	    mk_comb(wdbp, "dup_twice.c", &duplicate_wm.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) == 0;
    }
    if (ret) {
	struct wmember progressive_wm;
	mat_t first_mat;
	mat_t duplicate_mat;
	unsigned char progressive_rgb[3] = {42, 84, 126};
	MAT_IDN(first_mat);
	MAT_IDN(duplicate_mat);
	MAT_DELTAS(first_mat, 11.0, 12.0, 13.0);
	MAT_DELTAS(duplicate_mat, 21.0, 22.0, 23.0);
	BU_LIST_INIT(&progressive_wm.l);
	ret = mk_addmember("dup_leaf.s", &progressive_wm.l, first_mat,
		WMOP_UNION) != NULL &&
	    mk_addmember("dup_leaf.s", &progressive_wm.l, duplicate_mat,
		WMOP_UNION) != NULL &&
	    mk_addmember("box.s", &progressive_wm.l, NULL,
		WMOP_SUBTRACT) != NULL &&
	    mk_addmember("ball.s", &progressive_wm.l, NULL,
		WMOP_INTERSECT) != NULL &&
	    mk_comb(wdbp, "progressive_root.c", &progressive_wm.l, 0,
		NULL, NULL, progressive_rgb, 0, 0, 0, 0, 0, 0, 0) == 0;
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

struct command_result_callback_state {
    int callback_count;
    int accepted_count;
    int updated_count;
    int removed_count;
    int failed_count;
    int saw_line_layers_update;
    int saw_metadata_update;
    int saw_primitive_metadata_update;
    int saw_remove_prefix;
    int saw_stale_failure;
    int saw_commit_failure;
    uint64_t line_layers_feature_id;
    uint64_t metadata_feature_id;
    uint64_t primitive_metadata_feature_id;
};

static void
command_event_cb(const struct ged_view_feature_batch_event *result, void *data)
{
    struct command_result_callback_state *ctx =
	(struct command_result_callback_state *)data;
    if (!ctx || !result)
	return;

    ctx->callback_count++;
    switch (result->status) {
	case GED_VIEW_FEATURE_BATCH_ACCEPTED:
	    ctx->accepted_count++;
	    break;
	case GED_VIEW_FEATURE_BATCH_UPDATED:
	    ctx->updated_count++;
	    break;
	case GED_VIEW_FEATURE_BATCH_REMOVED:
	    ctx->removed_count++;
	    break;
	case GED_VIEW_FEATURE_BATCH_FAILED:
	    ctx->failed_count++;
	    break;
	default:
	    break;
    }

    const char *name = result->feature_name ? result->feature_name : "";
    const char *command = result->command ? result->command : "";
    if (result->status == GED_VIEW_FEATURE_BATCH_UPDATED &&
	    BU_STR_EQUAL(name, "rtcheck::overlaps") &&
	    BU_STR_EQUAL(command, "lineLayersReplace")) {
	ctx->saw_line_layers_update = 1;
	ctx->line_layers_feature_id = result->feature_id;
    }
    if (result->status == GED_VIEW_FEATURE_BATCH_UPDATED &&
	    BU_STR_EQUAL(name, "rtcheck::overlaps") &&
	    BU_STR_EQUAL(command, "metadataReplace")) {
	ctx->saw_metadata_update = 1;
	ctx->metadata_feature_id = result->feature_id;
    }
    if (result->status == GED_VIEW_FEATURE_BATCH_UPDATED &&
	    BU_STR_EQUAL(name, "rtcheck::overlaps") &&
	    BU_STR_EQUAL(command, "primitiveMetadataReplace")) {
	ctx->saw_primitive_metadata_update = 1;
	ctx->primitive_metadata_feature_id = result->feature_id;
    }
    if (result->status == GED_VIEW_FEATURE_BATCH_REMOVED &&
	    BU_STR_EQUAL(command, "removePrefix") &&
	    bu_strncmp(name, "rtcheck::", strlen("rtcheck::")) == 0)
	ctx->saw_remove_prefix = 1;
    if (result->status == GED_VIEW_FEATURE_BATCH_FAILED &&
	    BU_STR_EQUAL(name, "rtcheck::generation") &&
	    BU_STR_EQUAL(command, "lineLayersReplace"))
	ctx->saw_stale_failure = 1;
    if (result->status == GED_VIEW_FEATURE_BATCH_FAILED &&
	    BU_STR_EQUAL(command, "commit"))
	ctx->saw_commit_failure = 1;
}

static int
feature_overlay_matches(BObolViewController *controller,
	const char *name,
	BObolOverlayClass overlay_class,
	BObolOverlayLifecycle lifecycle,
	BObolOverlayOrder order)
{
    if (!controller || !name)
	return 0;

    BObolFeatureHandle handle = controller->features().find(name);
    BObolOverlayInfo overlay;
    if (!handle.isValid() ||
	    !controller->features().overlayInfo(handle, overlay))
	return 0;

    return overlay.isOverlay &&
	overlay.role == BObolOverlayRole::Model &&
	overlay.overlayClass == overlay_class &&
	overlay.lifecycle == lifecycle &&
	overlay.order == order &&
	BU_STR_EQUAL(overlay.sourcePath.getString(), name);
}

static SoBRLVListShape *
first_feature_vlist(SoNode *node)
{
    if (!node)
	return NULL;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<SoBRLVListShape *>(node);
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoBRLVListShape *found = first_feature_vlist(group->getChild(i));
	if (found)
	    return found;
    }
    return NULL;
}

static int
scene_has_source_geometry(SoNode *node, const char *instance_key,
	int mesh)
{
    if (!node || !instance_key || !instance_key[0])
	return 0;
    if (!mesh && node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	return BU_STR_EQUAL(shape->ownerSourceInstanceKey.getValue().getString(),
		instance_key);
    }
    if (mesh && node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	return BU_STR_EQUAL(shape->ownerSourceInstanceKey.getValue().getString(),
		instance_key);
    }
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return 0;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (scene_has_source_geometry(group->getChild(i), instance_key, mesh))
	    return 1;
    }
    return 0;
}

static SoBRLDatabaseSource *
source_for_path(BObolSceneController *controller, const char *path)
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

static int
seed_view_lod_probe_payload(BObolViewController *controller,
			    const char *path,
			    const char *name)
{
    if (!controller || !controller->getViewLodState() || !path || !name)
	return 0;

    SoBRLMeshShape *mesh = new SoBRLMeshShape;
    mesh->ref();
    mesh->sourcePath = path;
    mesh->sourceName = name;

    BObolLodRequest request;
    request.databaseId = "ged-obol-sync";
    request.sourceRevision = 1;
    request.sourceContentHash = 1;
    request.objectPath = path;
    request.objectName = name;
    request.viewRevision = 1;
    request.policyRevision = 1;
    request.drawMode = BOBOL_LOD_DRAW_WIRE;
    request.providerId = "ged-obol-sync-probe";
    request.providerVersion = "1";
    request.qualityTier = BOBOL_LOD_QUALITY_PROXY;
    request.bounds = SbBox3f(SbVec3f(-1.0f, -1.0f, -1.0f),
	    SbVec3f(1.0f, 1.0f, 1.0f));

    BObolLodCounts counts;
    counts.faceCount = 1;
    counts.pointCount = 2;
    BObolLodResult result =
	bobol_lod_aabb_result(request, request.bounds, &counts);

    const int seeded = controller->getViewLodState()->applyProxyResult(
	    mesh, result) &&
	controller->getViewLodState()->payloadCount() > 0;
    mesh->unref();
    return seeded;
}

static int
apply_attached_view_lod_invalidation_probe(struct ged *gedp,
	BObolViewController *controller,
	struct ged_scene_reducer_request *txn,
	const char *label)
{
    if (!gedp || !controller || !txn || !label)
	return 1;

    if (!seed_view_lod_probe_payload(controller, "box.s", "box.s")) {
	fprintf(stderr,
		"FAIL: attached Obol view-controller LoD invalidation probe should seed %s payload\n",
		label);
	return 1;
    }

    struct ged_scene_reducer_result result;
    ged_scene_reducer_result_init(&result);
    int ret = ged_scene_reduce(gedp, txn, &result);
    ged_scene_reducer_result_free(&result);
    if (ret <= 0) {
	fprintf(stderr,
		"FAIL: attached Obol %s transaction should succeed\n",
		label);
	return 1;
    }
    if (controller->getViewLodState()->payloadCount() != 0) {
	fprintf(stderr,
		"FAIL: attached Obol %s transaction should clear view-local LoD state\n",
		label);
	return 1;
    }

    return 0;
}

static int
source_instance_is_view_scoped(SoBRLDatabaseSource *source, const char *view_name)
{
    if (!source || !view_name || !view_name[0])
	return 0;
    std::string prefix("ged-view:");
    prefix += view_name;
    prefix += ":";
    const char *instance_key = source->instanceKey.getValue().getString();
    return instance_key && bu_strncmp(instance_key, prefix.c_str(),
	    prefix.length()) == 0;
}

static int
source_instance_is_any_view_scoped(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;
    const char *instance_key = source->instanceKey.getValue().getString();
    return instance_key && bu_strncmp(instance_key, "ged-view:", 9) == 0;
}

static SoBRLDatabaseSource *
source_for_view_path(BObolSceneController *controller,
	const char *view_name,
	const char *path)
{
    if (!controller || !view_name || !view_name[0] || !path)
	return NULL;
    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (source && path_equal(source->path.getValue().getString(), path) &&
		source_instance_is_view_scoped(source, view_name))
	    return source;
    }
    return NULL;
}

static SoBRLDatabaseSource *
source_for_shared_path(BObolSceneController *controller, const char *path)
{
    if (!controller || !path)
	return NULL;
    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (source && path_equal(source->path.getValue().getString(), path) &&
		!source_instance_is_any_view_scoped(source))
	    return source;
    }
    return NULL;
}

static SoBRLDatabaseSource *
source_for_representation(BObolSceneController *controller,
	const char *path,
	int representation_mode)
{
    if (!controller || !path)
	return NULL;

    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (!source)
	    continue;
	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    continue;
	if (path_equal(summary.path.getString(), path) &&
		summary.representationMode == representation_mode)
	    return source;
    }

    return NULL;
}

static int
source_representation_count(BObolSceneController *controller,
	const char *path,
	int representation_mode)
{
    if (!controller || !path)
	return 0;

    int count = 0;
    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (!source)
	    continue;
	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    continue;
	if (path_equal(summary.path.getString(), path) &&
		summary.representationMode == representation_mode)
	    count++;
    }

    return count;
}

static int
source_path_count(BObolSceneController *controller,
	const char *path)
{
    if (!controller || !path)
	return 0;

    int count = 0;
    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (source && path_equal(source->path.getValue().getString(), path))
	    count++;
    }

    return count;
}

static int
verify_mode_source(BObolSceneController *controller,
	const char *path,
	int representation_mode,
	int expect_vlist,
	int expect_mesh,
	int expect_visible,
	int expect_stale,
	const char *label)
{
    SoBRLDatabaseSource *source =
	source_for_representation(controller, path, representation_mode);
    if (!source)
	FAIL("mode-specific source should exist");

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	FAIL("mode-specific source summary should be readable");

    if (summary.representationMode != representation_mode)
	FAIL("mode-specific source should preserve exact representation mode");
    if ((summary.visible ? 1 : 0) != expect_visible)
	FAIL("mode-specific source visibility should match expected state");
    if ((summary.stale ? 1 : 0) != expect_stale) {
	FAIL("mode-specific source stale state should match expected state");
    }
    if (!expect_stale &&
	    summary.realizationStatus != SoBRLDatabaseSource::REALIZED)
	FAIL("mode-specific source should be realized when expected current");
    const char *instance_key = summary.instanceKey.getString();
    if (expect_vlist && !source->hasRealizedWireGeometry() &&
	    !scene_has_source_geometry(controller->getSceneRoot(), instance_key, 0))
	FAIL("mode-specific source should carry realized VLIST geometry");
    if (expect_mesh && !source->hasRealizedMeshGeometry() &&
	    !scene_has_source_geometry(controller->getSceneRoot(), instance_key, 1))
	FAIL("mode-specific source should carry realized mesh geometry");
    if (!summary.sourceBoundsValid || summary.sourceBounds.isEmpty())
	FAIL("mode-specific source should retain valid bounds");

    (void)label;
    return 0;
}

static int
apply_mode_value_transaction(struct ged *gedp,
	ged_scene_reducer_operation kind,
	const char *path,
	int mode,
	fastf_t value,
	const char *label)
{
    struct ged_scene_reducer_request txn =
	ged_scene_reducer_request_make_value(kind, path, value);
    txn.mode = mode;
    struct ged_scene_reducer_result result;
    ged_scene_reducer_result_init(&result);
    int ret = ged_scene_reduce(gedp, &txn, &result);
    ged_scene_reducer_result_free(&result);
    if (ret <= 0)
	FAIL("mode-specific value transaction should succeed");

    (void)label;
    return 0;
}

static int
apply_mode_path_transaction(struct ged *gedp,
	ged_scene_reducer_operation kind,
	const char *path,
	int mode,
	const char *label)
{
    struct ged_scene_reducer_request txn = ged_scene_reducer_request_make(kind, path);
    txn.mode = mode;
    if (kind == GED_SCENE_REDUCER_STALE_SOURCE)
	txn.stale_reason = GED_DRAW_STALE_SETTINGS_CHANGED;
    struct ged_scene_reducer_result result;
    ged_scene_reducer_result_init(&result);
    int ret = ged_scene_reduce(gedp, &txn, &result);
    ged_scene_reducer_result_free(&result);
    if (ret <= 0)
	FAIL("mode-specific path transaction should succeed");

    (void)label;
    return 0;
}

static int
apply_path_transaction(struct ged *gedp,
	ged_scene_reducer_operation kind,
	const char *path,
	struct ged_view_context *view_ctx,
	int mode,
	const char *label)
{
    struct ged_scene_reducer_request txn = ged_scene_reducer_request_make(kind, path);
    txn.view = view_ctx;
    txn.mode = mode;
    struct ged_scene_reducer_result result;
    ged_scene_reducer_result_init(&result);
    int ret = ged_scene_reduce(gedp, &txn, &result);
    const int status = result.status;
    ged_scene_reducer_result_free(&result);
    if (ret <= 0) {
	fprintf(stderr, "public path transaction failed: %s (kind=%d path=%s mode=%d ret=%d status=%d)\n",
	    label ? label : "", static_cast<int>(kind), path ? path : "",
	    mode, ret, status);
	FAIL("public path transaction should succeed");
    }

    (void)label;
    return 0;
}

static int
try_path_transaction(struct ged *gedp,
	ged_scene_reducer_operation kind,
	const char *path,
	struct ged_view_context *view_ctx,
	int mode)
{
    struct ged_scene_reducer_request txn = ged_scene_reducer_request_make(kind, path);
    txn.view = view_ctx;
    txn.mode = mode;
    struct ged_scene_reducer_result result;
    ged_scene_reducer_result_init(&result);
    int ret = ged_scene_reduce(gedp, &txn, &result);
    ged_scene_reducer_result_free(&result);
    if (ret < 0)
	FAIL("public path transaction should not fail");

    return ret;
}

static int
exercise_mode_specific_source_lifecycle(struct ged *gedp,
	BObolSceneController *controller,
	const char *path,
	int mode,
	int representation_mode,
	int expect_vlist,
	int expect_mesh,
	int expect_direct_provider,
	const char *label)
{
    if (!gedp || !controller || !path)
	FAIL("mode-specific lifecycle test needs GED and Obol scene state");
    (void)expect_direct_provider;

    (void)try_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, path,
	    ged_draw_active_view_ctx(gedp), mode);
    if (source_representation_count(controller, path, representation_mode))
	FAIL("mode-specific lifecycle setup should start without target representation");

    char mode_arg[16] = {0};
    snprintf(mode_arg, sizeof(mode_arg), "-m%d", mode);
    const char *draw_mode_cmd[4] = {"draw", mode_arg, path, NULL};
    if (ged_exec_draw(gedp, 3, draw_mode_cmd) != BRLCAD_OK)
	FAIL("mode-specific draw command should succeed");
    if (source_representation_count(controller, path, representation_mode) != 1)
	FAIL("mode-specific draw should create exactly one target representation source");
    if (verify_mode_source(controller, path, representation_mode,
	    expect_vlist, expect_mesh, 1, 0, label))
	return 1;

    if (apply_mode_value_transaction(gedp, GED_SCENE_REDUCER_VISIBILITY,
	    path, mode, 0.0, label))
	return 1;
    if (verify_mode_source(controller, path, representation_mode,
	    expect_vlist, expect_mesh, 0, 0, label))
	return 1;

    if (apply_mode_value_transaction(gedp, GED_SCENE_REDUCER_VISIBILITY,
	    path, mode, 1.0, label))
	return 1;
    if (apply_mode_path_transaction(gedp, GED_SCENE_REDUCER_STALE_SOURCE,
	    path, mode, label))
	return 1;
    if (verify_mode_source(controller, path, representation_mode,
	    expect_vlist, expect_mesh, 1, 1, label))
	return 1;

    if (apply_mode_path_transaction(gedp, GED_SCENE_REDUCER_REDRAW,
	    path, mode, label))
	return 1;
    if (verify_mode_source(controller, path, representation_mode,
	    expect_vlist, expect_mesh, 1, 0, label))
	return 1;

    const char *autoview_cmd[2] = {"autoview", NULL};
    if (ged_exec_autoview(gedp, 1, autoview_cmd) != BRLCAD_OK)
	FAIL("mode-specific autoview should succeed");
    if (verify_mode_source(controller, path, representation_mode,
	    expect_vlist, expect_mesh, 1, 0, label))
	return 1;

    if (!ged_draw_obol_database_source_ensure_for_path(gedp, path,
	    gedp->dbip, GED_DRAW_MODE_WIRE, 0))
	FAIL("mode-specific lifecycle should be able to ensure a shared wire source");
    if (source_path_count(controller, path) < 2)
	FAIL("mode-specific lifecycle should have shared and mode sources before scoped erase");

    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, path,
	    ged_draw_active_view_ctx(gedp), mode, "mode-specific erase"))
	return 1;
    if (source_representation_count(controller, path, representation_mode))
	FAIL("mode-specific erase should not leave target representation sources");
    if (!source_for_path(controller, path))
	FAIL("mode-specific erase should preserve the shared wire source");

    return 0;
}

static int
exercise_deferred_mode_replacement(struct ged *gedp,
	BObolSceneController *controller, const char *path)
{
    if (!gedp || !controller || !path)
	FAIL("deferred mode replacement test needs GED and Obol scene state");

    (void)try_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, path,
	ged_draw_active_view_ctx(gedp), -1);

    const char *draw_shaded[3] = {"draw", "-m2", path};
    if (ged_exec_draw(gedp, 3, draw_shaded) != BRLCAD_OK)
	FAIL("deferred shaded draw should succeed");
    if (source_representation_count(controller, path,
		SoBRLDatabaseSource::REPRESENTATION_SHADED) != 1 ||
	source_representation_count(controller, path,
		SoBRLDatabaseSource::REPRESENTATION_WIRE) != 0)
	FAIL("deferred shaded draw should install only its shaded representation");

    const char *draw_wire[3] = {"draw", "-m0", path};
    if (ged_exec_draw(gedp, 3, draw_wire) != BRLCAD_OK)
	FAIL("deferred wire draw should succeed");
    if (source_representation_count(controller, path,
		SoBRLDatabaseSource::REPRESENTATION_WIRE) != 1 ||
	source_representation_count(controller, path,
		SoBRLDatabaseSource::REPRESENTATION_SHADED) != 0)
	FAIL("normal deferred draw should replace an earlier representation");

    const char *add_shaded[4] = {"draw", "-m2", "--add-mode", path};
    if (ged_exec_draw(gedp, 4, add_shaded) != BRLCAD_OK)
	FAIL("deferred add-mode shaded draw should succeed");
    if (source_representation_count(controller, path,
		SoBRLDatabaseSource::REPRESENTATION_WIRE) != 1 ||
	source_representation_count(controller, path,
		SoBRLDatabaseSource::REPRESENTATION_SHADED) != 1)
	FAIL("deferred add-mode should retain the existing representation");

    (void)try_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, path,
	ged_draw_active_view_ctx(gedp), -1);
    return 0;
}

static int
box3f_near(const SbBox3f &box,
	float min_x,
	float min_y,
	float min_z,
	float max_x,
	float max_y,
	float max_z)
{
    const SbVec3f bmin = box.getMin();
    const SbVec3f bmax = box.getMax();
    return fabsf(bmin[0] - min_x) <= 0.001f &&
	fabsf(bmin[1] - min_y) <= 0.001f &&
	fabsf(bmin[2] - min_z) <= 0.001f &&
	fabsf(bmax[0] - max_x) <= 0.001f &&
	fabsf(bmax[1] - max_y) <= 0.001f &&
	fabsf(bmax[2] - max_z) <= 0.001f;
}

static int
exercise_multi_instance_transform_reuse(struct ged *gedp,
	BObolSceneController *controller)
{
    if (!gedp || !controller)
	FAIL("multi-instance transform test needs GED and Obol scene state");

    const int initial_source_count = controller->getDatabaseSourceCount();
    const char *path_a =
	"reuse_root.c/reuse_inst_a.c/reuse_shared.c/reuse_leaf.s";
    const char *path_b =
	"reuse_root.c/reuse_inst_b.c/reuse_shared.c@1/reuse_leaf.s@1";
    const char *draw_reuse_root[2] = {"draw", "reuse_root.c"};
    if (ged_exec_draw(gedp, 2, draw_reuse_root) != BRLCAD_OK)
	FAIL("GED multi-instance transform root draw should succeed");
    (void)controller->realizePending();

    auto compact_for_path = [](SoBRLDatabaseSource *source,
	    const char *path, BObolCompactInstanceHandle &handle,
	    BObolCompactInstanceSummary &summary) {
	if (!source || !path)
	    return false;
	for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	    BObolCompactInstanceHandle candidate;
	    BObolCompactInstanceSummary candidateSummary;
	    if (!source->getCompactInstanceHandle(i, candidate) ||
		!source->getCompactInstanceSummary(candidate, candidateSummary))
		continue;
	    if (path_equal(candidateSummary.path.getString(), path)) {
		handle = candidate;
		summary = candidateSummary;
		return true;
	    }
	}
	return false;
    };
    auto transformed_bounds = [](const BObolCompactInstanceSummary &summary) {
	SbBox3f bounds = summary.localBounds;
	bounds.transform(summary.localToSource);
	return bounds;
    };

    SoBRLDatabaseSource *wire_source = source_for_representation(controller,
	"reuse_root.c", SoBRLDatabaseSource::REPRESENTATION_WIRE);
    BObolCompactInstanceHandle handle_a;
    BObolCompactInstanceHandle handle_b;
    BObolCompactInstanceSummary compact_a;
    BObolCompactInstanceSummary compact_b;
    if (!wire_source || !wire_source->hasCompactInstanceIndex() ||
	    wire_source->getCompactInstanceCountForPath(path_a, FALSE) != 1 ||
	    wire_source->getCompactInstanceCountForPath(path_b, FALSE) != 1 ||
	    !compact_for_path(wire_source, path_a, handle_a, compact_a) ||
	    !compact_for_path(wire_source, path_b, handle_b, compact_b))
	FAIL("GED multi-instance root should retain both compact occurrences");
    if (handle_a.instanceWord0 == handle_b.instanceWord0 &&
	    handle_a.instanceWord1 == handle_b.instanceWord1)
	FAIL("GED multi-instance occurrences should have distinct stable handles");
    if (!compact_a.wireGeometry || !compact_b.wireGeometry ||
	    compact_a.geometryIdentity == 0 ||
	    compact_a.geometryIdentity != compact_b.geometryIdentity ||
	    !box3f_near(compact_a.localBounds, -1.0f, -1.0f, -1.0f,
		1.0f, 1.0f, 1.0f))
	FAIL("GED multi-instance occurrences should share source-local wire geometry");
    const SbBox3f bounds_a = transformed_bounds(compact_a);
    const SbBox3f bounds_b = transformed_bounds(compact_b);
    if (!box3f_near(bounds_a, -13.0f, 2.0f, -1.0f,
	    -11.0f, 4.0f, 1.0f) ||
	    !box3f_near(bounds_b, 17.0f, 2.0f, -1.0f,
		19.0f, 4.0f, 1.0f))
	FAIL("GED multi-instance occurrence transforms should place shared geometry");

    SbBox3f scene_bounds_a;
    SbBox3f scene_bounds_b;
    SbBox3f scene_bounds_root;
    const SbBool have_scene_bounds_a = controller->getSceneSubtreeBounds(
	path_a, TRUE, scene_bounds_a);
    if (!have_scene_bounds_a ||
	    !box3f_near(scene_bounds_a, -13.0f, 2.0f, -1.0f,
		-11.0f, 4.0f, 1.0f)) {
	fprintf(stderr, "first occurrence bounds valid=%d min=(%g,%g,%g) max=(%g,%g,%g)\n",
	    have_scene_bounds_a, scene_bounds_a.getMin()[0],
	    scene_bounds_a.getMin()[1], scene_bounds_a.getMin()[2],
	    scene_bounds_a.getMax()[0], scene_bounds_a.getMax()[1],
	    scene_bounds_a.getMax()[2]);
	FAIL("GED multi-instance first occurrence bounds should apply its transform");
    }
    if (!controller->getSceneSubtreeBounds(path_b, TRUE, scene_bounds_b) ||
	    !box3f_near(scene_bounds_b, 17.0f, 2.0f, -1.0f,
		19.0f, 4.0f, 1.0f))
	FAIL("GED multi-instance second occurrence bounds should apply its transform");
    if (!controller->getSceneSubtreeBounds("reuse_root.c", TRUE,
	    scene_bounds_root) ||
	    !box3f_near(scene_bounds_root, -13.0f, 2.0f, -1.0f,
		19.0f, 4.0f, 1.0f))
	FAIL("GED multi-instance root bounds should include both occurrences");

    if (apply_mode_value_transaction(gedp, GED_SCENE_REDUCER_VISIBILITY,
	    path_a, GED_DRAW_MODE_WIRE, 0.0, "multi-instance visibility"))
	return 1;
    if (!wire_source->getCompactInstanceSummary(handle_a, compact_a) ||
	    !wire_source->getCompactInstanceSummary(handle_b, compact_b) ||
	    compact_a.visible || !compact_b.visible)
	FAIL("GED multi-instance visibility should target one occurrence path");
    if (apply_mode_value_transaction(gedp, GED_SCENE_REDUCER_VISIBILITY,
	    path_a, GED_DRAW_MODE_WIRE, 1.0, "multi-instance visibility restore"))
	return 1;
    if (!wire_source->getCompactInstanceSummary(handle_a, compact_a) ||
	    !wire_source->getCompactInstanceSummary(handle_b, compact_b) ||
	    !compact_a.visible || !compact_b.visible)
	FAIL("GED multi-instance visibility restore should update repeated logical instances consistently");

    if (apply_mode_value_transaction(gedp, GED_SCENE_REDUCER_HIGHLIGHT,
	    "reuse_root.c", GED_DRAW_MODE_WIRE, 1.0,
	    "multi-instance highlight"))
	return 1;
    if (!wire_source->getCompactInstanceSummary(handle_a, compact_a) ||
	    !wire_source->getCompactInstanceSummary(handle_b, compact_b) ||
	    !compact_a.highlighted || !compact_b.highlighted)
	FAIL("GED multi-instance highlight should update repeated logical instances consistently");
    if (apply_mode_value_transaction(gedp, GED_SCENE_REDUCER_HIGHLIGHT,
	    "reuse_root.c", GED_DRAW_MODE_WIRE, 0.0,
	    "multi-instance highlight restore"))
	return 1;
    if (!wire_source->getCompactInstanceSummary(handle_a, compact_a) ||
	    !wire_source->getCompactInstanceSummary(handle_b, compact_b) ||
	    compact_a.highlighted || compact_b.highlighted)
	FAIL("GED multi-instance highlight restore should update repeated logical instances consistently");

    ged_draw_index_stats_reset(gedp);
    if (apply_mode_path_transaction(gedp, GED_SCENE_REDUCER_REDRAW,
	    path_a, GED_DRAW_MODE_WIRE, "multi-instance redraw"))
	return 1;
    struct ged_draw_index_stats redraw_stats;
    memset(&redraw_stats, 0, sizeof(redraw_stats));
    ged_draw_index_stats_get(gedp, &redraw_stats);
    if (redraw_stats.slow_path_shape_scans ||
	    redraw_stats.slow_path_group_scans)
	FAIL("GED multi-instance logical redraw should avoid registry/index slow-path scans");
    wire_source = source_for_representation(controller, "reuse_root.c",
	    SoBRLDatabaseSource::REPRESENTATION_WIRE);
    BObolCompactInstanceHandle redraw_handle_a;
    BObolCompactInstanceHandle redraw_handle_b;
    BObolCompactInstanceSummary redraw_a;
    BObolCompactInstanceSummary redraw_b;
    if (!compact_for_path(wire_source, path_a, redraw_handle_a, redraw_a) ||
	    !compact_for_path(wire_source, path_b, redraw_handle_b, redraw_b) ||
	    redraw_handle_a.instanceWord0 != handle_a.instanceWord0 ||
	    redraw_handle_a.instanceWord1 != handle_a.instanceWord1 ||
	    redraw_handle_b.instanceWord0 != handle_b.instanceWord0 ||
	    redraw_handle_b.instanceWord1 != handle_b.instanceWord1 ||
	    redraw_a.geometryIdentity != compact_a.geometryIdentity ||
	    redraw_b.geometryIdentity != compact_b.geometryIdentity)
	FAIL("GED multi-instance redraw should preserve handles and shared geometry");

    const char *draw_reuse_root_shaded[3] = {
	"draw", "-m2", "reuse_root.c"
    };
    if (ged_exec_draw(gedp, 3, draw_reuse_root_shaded) != BRLCAD_OK)
	FAIL("GED multi-instance shaded root draw should succeed");
    SoBRLDatabaseSource *mesh_source = source_for_representation(controller,
	"reuse_root.c", SoBRLDatabaseSource::REPRESENTATION_SHADED);
    BObolCompactInstanceHandle mesh_handle_a;
    BObolCompactInstanceHandle mesh_handle_b;
    BObolCompactInstanceSummary mesh_a;
    BObolCompactInstanceSummary mesh_b;
    if (!mesh_source || !mesh_source->hasCompactInstanceIndex() ||
	    !compact_for_path(mesh_source, path_a, mesh_handle_a, mesh_a) ||
	    !compact_for_path(mesh_source, path_b, mesh_handle_b, mesh_b) ||
	    !mesh_a.meshGeometry || !mesh_b.meshGeometry ||
	    mesh_a.geometryIdentity == 0 ||
	    mesh_a.geometryIdentity != mesh_b.geometryIdentity ||
	    !box3f_near(transformed_bounds(mesh_a),
		-13.0f, 2.0f, -1.0f, -11.0f, 4.0f, 1.0f) ||
	    !box3f_near(transformed_bounds(mesh_b),
		17.0f, 2.0f, -1.0f, 19.0f, 4.0f, 1.0f))
	FAIL("GED multi-instance shaded draw should retain transformed mesh occurrences");

    const char *autoview_cmd[2] = {"autoview", NULL};
    point_t autoview_min;
    point_t autoview_max;
    int autoview_empty = 1;
    if (!ged_draw_obol_scene_database_autoview_bounds(gedp,
	    &autoview_min, &autoview_max, &autoview_empty, 1) ||
	autoview_empty || fabs(autoview_min[X] + 13.0) > 0.01 ||
	fabs(autoview_min[Y] + 1.0) > 0.01 ||
	fabs(autoview_min[Z] + 1.0) > 0.01 ||
	fabs(autoview_max[X] - 19.0) > 0.01 ||
	fabs(autoview_max[Y] - 4.0) > 0.01 ||
	fabs(autoview_max[Z] - 1.0) > 0.01)
	FAIL("GED multi-instance autoview bounds should retain the raw scene AABB");
    if (ged_exec_autoview(gedp, 1, autoview_cmd) != BRLCAD_OK)
	FAIL("GED multi-instance autoview should succeed");
    mat_t view_center_mat;
    point_t view_center;
    bv_center_mat_get(view_center_mat, DRAW_TEST_BV_CONST(ged_draw_active_view_ctx(gedp)));
    MAT_DELTAS_GET_NEG(view_center, view_center_mat);

    /* bv_autoview_bounds applies one rotation-stable bounding-sphere fit to
     * the raw scene AABB.  Bounds delivery must not pre-expand another cube,
     * which would apply diagonal margin twice and over-scale this view by
     * sqrt(3). */
    const struct bv *active_view =
	DRAW_TEST_BV_CONST(ged_draw_active_view_ctx(gedp));
    point_t expected_center;
    vect_t half_diagonal;
    VADD2SCALE(expected_center, autoview_min, autoview_max, 0.5);
    VSUB2(half_diagonal, autoview_max, expected_center);
    fastf_t expected_view_size = 2.0 * MAGNITUDE(half_diagonal);
    const int view_width = bv_width_get(active_view);
    const int view_height = bv_height_get(active_view);
    if (view_width > view_height && view_height > 0)
	expected_view_size *= static_cast<fastf_t>(view_width) / view_height;

    const fastf_t reuse_view_size = bv_size_get(active_view);
    if (fabs(reuse_view_size - expected_view_size) > 0.1 ||
	fabs(view_center[X] - expected_center[X]) > 0.1 ||
	fabs(view_center[Y] - expected_center[Y]) > 0.1 ||
	fabs(view_center[Z] - expected_center[Z]) > 0.1) {
	fprintf(stderr,
	    "multi-instance autoview size=%g center=(%g,%g,%g)\n",
	    reuse_view_size, view_center[X], view_center[Y], view_center[Z]);
	FAIL("GED multi-instance autoview should use transformed scene bounds");
	}

    const char *erase_reuse_root[2] = {"erase", "reuse_root.c"};
    if (ged_exec_erase(gedp, 2, erase_reuse_root) != BRLCAD_OK)
	FAIL("GED multi-instance transform root erase should succeed");
    if (source_for_representation(controller, "reuse_root.c",
		SoBRLDatabaseSource::REPRESENTATION_WIRE) ||
	    source_for_representation(controller, "reuse_root.c",
		SoBRLDatabaseSource::REPRESENTATION_SHADED) ||
	    controller->getDatabaseSourceCount() != initial_source_count)
	FAIL("GED multi-instance transform cleanup should restore prior source state");

    return 0;
}

static int
exercise_duplicate_occurrence_pick_identity(struct ged *gedp,
	BObolSceneController *controller)
{
    if (!gedp || !controller)
	FAIL("duplicate occurrence pick test needs GED and Obol scene state");

    const int initial_source_count = controller->getDatabaseSourceCount();
    const char *path_a = "dup_twice.c/dup_leaf.s";
    const char *path_b = "dup_twice.c/dup_leaf.s@1";
    const char *draw_duplicate[2] = {"draw", "dup_twice.c"};
    if (ged_exec_draw(gedp, 2, draw_duplicate) != BRLCAD_OK)
	FAIL("GED duplicate occurrence draw should succeed");

    SoBRLDatabaseSource *aggregate = source_for_path(controller,
	"dup_twice.c");
    if (aggregate && aggregate->hasCompactInstanceIndex()) {
	if (aggregate->getCompactInstanceCountForPath(path_a, FALSE) != 1 ||
	    aggregate->getCompactInstanceCountForPath(path_b, FALSE) != 1)
	    FAIL("GED duplicate occurrences should have distinct registry paths");
	BObolCompactInstanceHandle handle_a;
	BObolCompactInstanceHandle handle_b;
	BObolCompactInstanceSummary summary_a;
	BObolCompactInstanceSummary summary_b;
	SbString key_a;
	SbString key_b;
	for (int i = 0; i < aggregate->getCompactInstanceCount(); i++) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    if (!aggregate->getCompactInstanceHandle(i, handle) ||
		!aggregate->getCompactInstanceSummary(handle, summary))
		FAIL("GED duplicate registry should expose valid handles");
	    if (path_equal(summary.path.getString(), path_a)) {
		handle_a = handle;
		summary_a = summary;
		key_a = summary.sourceInstanceKey;
	    }
	    if (path_equal(summary.path.getString(), path_b)) {
		handle_b = handle;
		summary_b = summary;
		key_b = summary.sourceInstanceKey;
	    }
	}
	if (key_a.getLength() == 0 || key_b.getLength() == 0) {
	    FAIL("GED duplicate registry occurrences should expose identities");
	}
	if (key_a == key_b)
	    FAIL("GED duplicate registry occurrence identities should be distinct");
	const int changed =
	    aggregate->setCompactInstanceDisplayStateForPath(path_a, FALSE,
		1, FALSE, 0, FALSE, 0, FALSE);
	const int got_a =
	    aggregate->getCompactInstanceSummary(handle_a, summary_a);
	const int got_b =
	    aggregate->getCompactInstanceSummary(handle_b, summary_b);
	if (changed <= 0 || !got_a || !got_b ||
	    summary_a.visible == summary_b.visible) {
	    fprintf(stderr,
		"duplicate visibility changed=%d got_a=%d got_b=%d "
		"visible_a=%d visible_b=%d path_a=%s path_b=%s\n",
		changed, got_a, got_b, summary_a.visible, summary_b.visible,
		summary_a.path.getString(), summary_b.path.getString());
	    FAIL("GED duplicate registry visibility should target one occurrence");
	}
	const char *visible_path = summary_a.visible ? path_a : path_b;
	const SbString &visible_key = summary_a.visible ? key_a : key_b;
	SoBRLExportAction export_action;
	export_action.apply(controller->getSceneRoot());
	SbVec3f midpoint;
	SbBool found_visible_segment = FALSE;
	for (int i = 0; i < export_action.getLineCount(); i++) {
	    const SoBRLExportAction::LineRecord &line = export_action.getLine(i);
	    if (!path_equal(line.path.getString(), visible_path))
		continue;
	    midpoint = (line.a + line.b) * 0.5f;
	    found_visible_segment = TRUE;
	    break;
	}
	if (!found_visible_segment)
	    FAIL("GED duplicate registry export should expose the visible occurrence");
	SbViewportRegion viewport(200, 200);
	SoRayPickAction pick_action(viewport);
	pick_action.setRay(SbVec3f(midpoint[0], midpoint[1],
	    midpoint[2] + 10.0f), SbVec3f(0.0f, 0.0f, -1.0f));
	pick_action.apply(controller->getSceneRoot());
	const SoPickedPoint *picked_point = pick_action.getPickedPoint();
	const SoDetail *raw_detail = picked_point ? picked_point->getDetail() : NULL;
	if (!raw_detail ||
	    !raw_detail->isOfType(SoBRLPickDetail::getClassTypeId()))
	    FAIL("GED duplicate registry pick should return BRL-CAD detail");
	const SoBRLPickDetail *pick_detail =
	    static_cast<const SoBRLPickDetail *>(raw_detail);
	if (!path_equal(pick_detail->getPath().getString(), visible_path) ||
	    !BU_STR_EQUAL(pick_detail->getSourceInstanceKey().getString(),
		visible_key.getString()))
	    FAIL("GED duplicate registry pick should identify the visible occurrence");

	const char *erase_duplicate[2] = {"erase", "dup_twice.c"};
	if (ged_exec_erase(gedp, 2, erase_duplicate) != BRLCAD_OK ||
	    controller->getDatabaseSourceCount() != initial_source_count)
	    FAIL("GED duplicate registry cleanup should restore prior state");
	return 0;
    }

    SoBRLDatabaseSource *source_a = source_for_representation(controller,
	    path_a, SoBRLDatabaseSource::REPRESENTATION_WIRE);
    SoBRLDatabaseSource *source_b = source_for_representation(controller,
	    path_b, SoBRLDatabaseSource::REPRESENTATION_WIRE);
    if (!source_a || !source_b)
	FAIL("GED duplicate occurrence draw should preserve both source instances");

    BObolDatabaseSourceSummary summary_a;
    BObolDatabaseSourceSummary summary_b;
    if (!source_a->getSummary(summary_a) || !summary_a.valid ||
	    !source_b->getSummary(summary_b) || !summary_b.valid ||
	    BU_STR_EQUAL(summary_a.instanceKey.getString(),
		summary_b.instanceKey.getString()))
	FAIL("GED duplicate occurrence summaries should report distinct instance keys");

    SoBRLVListShape *shape_a = source_a->getRealizedShape();
    SoBRLVListShape *shape_b = source_b->getRealizedShape();
    const SoBRLVListShape *geometry_a =
	shape_a ? shape_a->getGeometrySource() : NULL;
    const SoBRLVListShape *geometry_b =
	shape_b ? shape_b->getGeometrySource() : NULL;
    if (!shape_a || !shape_b || !geometry_a || !geometry_b ||
	    geometry_a == shape_a || geometry_b == shape_b ||
	    geometry_a != geometry_b)
	FAIL("GED duplicate occurrences should remain separate sources sharing one geometry node");
    if (bu_strcmp(shape_a->ownerSourceInstanceKey.getValue().getString(),
	    summary_a.instanceKey.getString()) != 0 ||
	    bu_strcmp(shape_b->ownerSourceInstanceKey.getValue().getString(),
		summary_b.instanceKey.getString()) != 0)
	FAIL("GED duplicate occurrence shapes should retain owner source instance keys");

    if (controller->setDatabaseSourceInstanceState(
	    summary_a.instanceKey.getString(),
	    TRUE, summary_a.sourceRevision, summary_a.inputsRevision,
	    FALSE, summary_a.selected, summary_a.highlighted,
	    summary_a.lineStyle,
	    summary_a.lineWidth, summary_a.transparency,
	    summary_a.colorOverride, summary_a.color,
	    summary_a.materialColorValid, summary_a.materialColor,
	    summary_a.materialRevision) < 0)
	FAIL("GED duplicate occurrence source visibility update should target one instance");
    if (!source_a->getSummary(summary_a) || !summary_a.valid ||
	    !source_b->getSummary(summary_b) || !summary_b.valid ||
	    summary_a.visible || !summary_b.visible ||
	    shape_a->visible.getValue() || !shape_b->visible.getValue())
	FAIL("GED duplicate occurrence visibility should hide only the targeted instance");

    SbVec3f segment_a;
    SbVec3f segment_b;
    if (!shape_b->getSegment(0, segment_a, segment_b))
	FAIL("GED duplicate occurrence pick fixture should expose a wire segment");
    SbVec3f midpoint(
	0.5f * (segment_a[0] + segment_b[0]),
	0.5f * (segment_a[1] + segment_b[1]),
	0.5f * (segment_a[2] + segment_b[2]));

    SbViewportRegion viewport(200, 200);
    SoRayPickAction pick_action(viewport);
    pick_action.setRay(
	SbVec3f(midpoint[0], midpoint[1], midpoint[2] + 10.0f),
	SbVec3f(0.0f, 0.0f, -1.0f));
    pick_action.apply(controller->getSceneRoot());
    const SoPickedPoint *picked_point = pick_action.getPickedPoint();
    if (!picked_point)
	FAIL("GED duplicate occurrence ray pick should hit the visible second instance");
    const SoDetail *raw_detail = picked_point->getDetail(shape_b);
    if (!raw_detail)
	raw_detail = picked_point->getDetail();
    if (!raw_detail ||
	    !raw_detail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("GED duplicate occurrence ray pick should return BRL-CAD pick detail");
    const SoBRLPickDetail *pick_detail =
	static_cast<const SoBRLPickDetail *>(raw_detail);
    if (bu_strcmp(pick_detail->getPath().getString(),
	    summary_b.path.getString()) != 0 ||
	    bu_strcmp(pick_detail->getSourceInstanceKey().getString(),
		summary_b.instanceKey.getString()) != 0 ||
	    BU_STR_EQUAL(pick_detail->getSourceInstanceKey().getString(),
		summary_a.instanceKey.getString()))
	FAIL("GED duplicate occurrence pick detail should identify the visible source instance");

    if (controller->setDatabaseSourceInstanceState(
	    summary_a.instanceKey.getString(),
	    TRUE, summary_a.sourceRevision, summary_a.inputsRevision,
	    TRUE, summary_a.selected, summary_a.highlighted,
	    summary_a.lineStyle,
	    summary_a.lineWidth, summary_a.transparency,
	    summary_a.colorOverride, summary_a.color,
	    summary_a.materialColorValid, summary_a.materialColor,
	    summary_a.materialRevision) < 0)
	FAIL("GED duplicate occurrence source visibility restore should succeed");

    const char *erase_duplicate[2] = {"erase", "dup_twice.c"};
    if (ged_exec_erase(gedp, 2, erase_duplicate) != BRLCAD_OK)
	FAIL("GED duplicate occurrence erase should succeed");
    if (source_for_representation(controller, path_a,
		SoBRLDatabaseSource::REPRESENTATION_WIRE) ||
	    source_for_representation(controller, path_b,
		SoBRLDatabaseSource::REPRESENTATION_WIRE) ||
	    controller->getDatabaseSourceCount() != initial_source_count)
	FAIL("GED duplicate occurrence cleanup should restore prior source state");

    return 0;
}

static int
exercise_progressive_occurrence_and_boolean_identity(struct ged *gedp,
	BObolSceneController *controller)
{
    if (!gedp || !controller)
	FAIL("progressive identity test needs GED and Obol scene state");
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx || !ged_view_context_display_endpoint_ensure(view_ctx))
	FAIL("progressive identity test needs an attached display endpoint");

    const int initial_source_count = controller->getDatabaseSourceCount();
    struct ged_draw_appearance_settings appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    appearance.defer_leaf_expansion = 1;
    struct ged_scene_reducer_request txn =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_DRAW, "progressive_root.c");
    txn.view = view_ctx;
    txn.mode = GED_DRAW_MODE_WIRE;
    txn.appearance = &appearance;
    struct ged_scene_reducer_result result;
    ged_scene_reducer_result_init(&result);
    int rejected_draw = 0;
    bool startup_escaped = false;
    {
	ScopedTransactionFault fault(
	    BObolTransactionFaultPoint::SOURCE_REALIZATION_WORKER_START);
	try {
	    rejected_draw = ged_scene_reduce(gedp, &txn, &result);
	} catch (const std::exception &) {
	    startup_escaped = true;
	}
    }
    ged_scene_reducer_result_free(&result);
    if (startup_escaped)
	FAIL("source worker startup failure escaped the GED draw transaction");
    SoBRLDatabaseSource *retained_proxy =
	source_for_path(controller, "progressive_root.c");
    BObolViewController *rejected_controller = ged_bobol_view_controller(view_ctx);
    BObolProgressiveOptions rejected_options;
    if (rejected_controller)
	(void)rejected_controller->advanceProgressiveWork(&rejected_options, NULL);
    BObolLodConvergenceStatus rejected_status;
    if (rejected_controller)
	rejected_controller->getLodConvergenceStatus(rejected_status);
    SbBox3f retained_bounds;
    if (rejected_draw <= 0 || !retained_proxy || !rejected_controller ||
	retained_proxy->getCompactInstanceCount() != 1 ||
	!retained_proxy->getSourceBounds(retained_bounds) || retained_bounds.isEmpty() ||
	rejected_status.sourcePreparationPending)
	FAIL("source startup denial should retain bounded coverage without a pending producer");
    if (try_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, "progressive_root.c",
	    view_ctx, -1) <= 0)
	FAIL("source startup denial should permit erasing its retained proxy");

    /* A completed worker may still have all its output queued. Denying that
     * owner's merge must retain the preview, without bypassing the denial by
     * installing the completed worker's full registry during adoption. */
    ged_scene_reducer_result_init(&result);
    const int denied_draw = ged_scene_reduce(gedp, &txn, &result);
    ged_scene_reducer_result_free(&result);
    SoBRLDatabaseSource *denied_source =
	source_for_path(controller, "progressive_root.c");
    BObolViewController *denied_controller = ged_bobol_view_controller(view_ctx);
    if (denied_draw <= 0 || !denied_source || !denied_controller ||
	denied_source->getCompactInstanceCount() != 1)
	FAIL("stream denial fixture should begin with one retained root preview");
    BObolCompactOccurrence denied_preview;
    if (!denied_source->getCompactOccurrence(0, denied_preview) ||
	!denied_preview.geometry || !denied_preview.summary.visible)
	FAIL("stream denial fixture should begin with visible preview geometry");
    const auto wire_bounds = [](const Obol::PartGeometry &geometry,
	const SbMatrix &transform) {
	SbBox3f bounds;
	if (geometry.wire) {
	    for (const SbVec3f &point : geometry.wire->segmentPoints) {
		SbVec3f transformed;
		transform.multVecMatrix(point, transformed);
		bounds.extendBy(transformed);
	    }
	}
	return bounds;
    };
    const SbBox3f intended_preview_bounds(denied_preview.summary.lodBoundsMin,
	denied_preview.summary.lodBoundsMax);
    const SbVec3f intended_min = intended_preview_bounds.getMin();
    const SbVec3f intended_max = intended_preview_bounds.getMax();
    const SbBox3f actual_preview_bounds = wire_bounds(*denied_preview.geometry,
	denied_preview.geometryTransform);
    if (actual_preview_bounds.isEmpty() || !box3f_near(actual_preview_bounds,
	intended_min[0], intended_min[1], intended_min[2],
	intended_max[0], intended_max[1], intended_max[2])) {
	fprintf(stderr, "preview extent: actual=(%g,%g,%g)-(%g,%g,%g) "
	    "intended=(%g,%g,%g)-(%g,%g,%g)\n",
	    actual_preview_bounds.getMin()[0], actual_preview_bounds.getMin()[1],
	    actual_preview_bounds.getMin()[2], actual_preview_bounds.getMax()[0],
	    actual_preview_bounds.getMax()[1], actual_preview_bounds.getMax()[2],
	    intended_min[0], intended_min[1], intended_min[2],
	    intended_max[0], intended_max[1], intended_max[2]);
	FAIL("retained overview wire vertices should preserve their intended extent");
    }
    const uint64_t preview_population = denied_source->getCompactPopulationEpoch();
    BObolViewLodState denied_view_state;
    const std::vector<const BObolViewLodState::CadPayload *> no_denied_payloads;
    (void)denied_source->compactViewLodAssembly(no_denied_payloads, &denied_view_state);
    SoCADAssembly *initial_presentation =
	denied_view_state.findCadPresentation(denied_source);
    const auto initial_ids = initial_presentation ? initial_presentation->instanceIds() :
	std::vector<Obol::InstanceId>();
    const auto initial_record = initial_ids.size() == 1 ?
	initial_presentation->getInstanceRecord(initial_ids[0]) :
	std::optional<Obol::InstanceRecord>();
    if (!initial_record || !box3f_near(wire_bounds(*denied_preview.geometry,
	initial_record->localToRoot), intended_min[0], intended_min[1],
	intended_min[2], intended_max[0], intended_max[1], intended_max[2]))
	FAIL("retained assembly should apply the initial overview transform exactly once");
    BObolProgressiveOptions denied_options;
    BObolLodConvergenceStatus denied_status;
    BObolCompactOccurrence certified_preview;
    struct StreamCoverageObserver {
	enum { callbackWidth = 31 };
	BObolSceneController &scene;
	SoBRLDatabaseSource &source;
	SbString key;
	uint64_t frame;
	bool called = false;
	bool partial = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<StreamCoverageObserver *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    SbBox3f bounds;
	    self.partial =
		self.scene.findDatabaseSourceInstance(self.key.getString()) !=
		    &self.source ||
		self.scene.getFrameRevision() <= self.frame ||
		!self.source.hasExactSourceBounds() ||
		!self.source.getSourceBounds(bounds) || bounds.isEmpty();
	    BObolDatabaseSourceDisplayPatch patch;
	    patch.lineWidthValid = TRUE;
	    patch.lineWidth = callbackWidth;
	    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(
		self.key.getString(), patch) <= 0;
	}
    } coverageObserver{*controller, *denied_source,
	denied_source->instanceKey.getValue(), controller->getFrameRevision()};
    SoFieldSensor coverageSensor(StreamCoverageObserver::changed,
	&coverageObserver);
    coverageSensor.setPriority(0);
    coverageSensor.attach(&denied_source->sourceBoundsExact);
    {
	ScopedTransactionFault fault(
	    BObolTransactionFaultPoint::SOURCE_STREAM_MERGE_AFTER_COMPLETION);
	/* Drain just the priority overview before attempting any leaf. This
	 * captures the last successful commit independently of the denied merge. */
	denied_options.maxProviderItems = 1;
	const auto coverage_deadline = std::chrono::steady_clock::now() +
	    std::chrono::seconds(2);
	do {
	    (void)denied_controller->advanceProgressiveWork(&denied_options, NULL);
	    denied_controller->getLodConvergenceStatus(denied_status);
	    if (denied_source->hasExactSourceBounds() || denied_status.failedSourceCount)
		break;
	    std::this_thread::sleep_for(std::chrono::milliseconds(1));
	} while (std::chrono::steady_clock::now() < coverage_deadline);
	if (denied_status.failedSourceCount || !denied_status.sourcePreparationPending ||
	    denied_source->getCompactInstanceCount() != 1 ||
	    !denied_source->hasExactSourceBounds() ||
	    !denied_source->getCompactOccurrence(0, certified_preview))
	    FAIL("priority overview should publish before leaf delivery is attempted");
	const auto deadline = std::chrono::steady_clock::now() +
	    std::chrono::seconds(2);
	do {
	    (void)denied_controller->advanceProgressiveWork(&denied_options, NULL);
	    denied_controller->getLodConvergenceStatus(denied_status);
	    if (!denied_status.sourcePreparationPending)
		break;
	    std::this_thread::sleep_for(std::chrono::milliseconds(1));
	} while (std::chrono::steady_clock::now() < deadline);
    }
    coverageSensor.detach();
    if (!coverageObserver.called || coverageObserver.partial ||
	coverageObserver.failed || denied_source->lineWidth.getValue() !=
	    StreamCoverageObserver::callbackWidth)
	FAIL("stream coverage callback should see complete state and retain its later edit");
    BObolCompactOccurrence after_denial;
    const bool have_denied_preview =
	denied_source->getCompactOccurrence(0, after_denial);
    if (denied_status.sourcePreparationPending ||
	denied_status.failedSourceCount != 1 || denied_status.viewReady ||
	denied_source->realizationStatus.getValue() != SoBRLDatabaseSource::FAILED ||
	denied_source->realizationDiagnostic.getValue().getLength() == 0 ||
	denied_source->getCompactInstanceCount() != 1 ||
	denied_source->getCompactExpectedInstanceCount() <= 1 ||
	denied_source->getCompactPopulationEpoch() != preview_population ||
	!have_denied_preview ||
	after_denial.geometry != certified_preview.geometry ||
	after_denial.geometryTransform != certified_preview.geometryTransform ||
	!after_denial.summary.visible ||
	!denied_source->getSourceBounds(retained_bounds) || retained_bounds.isEmpty()) {
	fprintf(stderr, "stream denial: pending=%d count=%d population=%llu/%llu "
	    "same-geometry=%d visible=%d\n",
	    denied_status.sourcePreparationPending,
	    denied_source->getCompactInstanceCount(),
	    static_cast<unsigned long long>(denied_source->getCompactPopulationEpoch()),
	    static_cast<unsigned long long>(preview_population),
	    after_denial.geometry == certified_preview.geometry ? 1 : 0,
	    after_denial.summary.visible);
	FAIL("completed worker adoption should preserve a memory-denied stream preview");
    }
    (void)denied_source->compactViewLodAssembly(no_denied_payloads, &denied_view_state);
    SoCADAssembly *denied_presentation =
	denied_view_state.findCadPresentation(denied_source);
    const std::vector<Obol::InstanceId> denied_ids = denied_presentation ?
	denied_presentation->instanceIds() : std::vector<Obol::InstanceId>();
    const std::optional<Obol::InstanceRecord> denied_record = denied_ids.size() == 1 ?
	denied_presentation->getInstanceRecord(denied_ids[0]) :
	std::optional<Obol::InstanceRecord>();
    if (!denied_record || denied_presentation->isInstanceHidden(denied_ids[0]) ||
	denied_presentation->partGeometry(denied_record->part) != certified_preview.geometry.get())
	FAIL("source delivery failure should keep its actual retained preview drawable");
    const SbVec3f source_min = retained_bounds.getMin();
    const SbVec3f source_max = retained_bounds.getMax();
    const SbBox3f actual_retained_bounds = wire_bounds(*after_denial.geometry,
	denied_record->localToRoot);
    if (!denied_source->hasExactSourceBounds() ||
	!box3f_near(actual_retained_bounds, source_min[0], source_min[1], source_min[2],
	    source_max[0], source_max[1], source_max[2])) {
	fprintf(stderr, "retained coverage: exact=%d actual=(%g,%g,%g)-(%g,%g,%g) "
	    "source=(%g,%g,%g)-(%g,%g,%g) identity-transform=%d\n",
	    denied_source->hasExactSourceBounds(),
	    actual_retained_bounds.getMin()[0], actual_retained_bounds.getMin()[1],
	    actual_retained_bounds.getMin()[2], actual_retained_bounds.getMax()[0],
	    actual_retained_bounds.getMax()[1], actual_retained_bounds.getMax()[2],
	    source_min[0], source_min[1], source_min[2],
	    source_max[0], source_max[1], source_max[2],
	    denied_record->localToRoot == SbMatrix::identity());
	FAIL("certified source overview should survive denied leaf delivery at its full extent");
    }

    ged_view_lod_policy original_policy = BV_LOD_POLICY_INIT;
    if (!ged_view_lod_policy_get(&original_policy, view_ctx))
	FAIL("stream denial fixture should expose its display policy");
    if (!denied_controller->syncCameraFromViewContext(view_ctx))
	FAIL("stream denial fixture should establish its display camera");
    denied_controller->setViewportSize(160, 120);
    denied_controller->setLodControlTransitionTracing(TRUE);
    for (const bool automatic : {false, true, false}) {
	ged_view_lod_policy policy = original_policy;
	policy.policy = automatic ? BV_LOD_AUTO : BV_LOD_OFF;
	policy.mesh_enabled = automatic ? 1 : 0;
	policy.csg_enabled = automatic ? 1 : 0;
	if (!ged_view_lod_policy_apply(view_ctx, &policy))
	    FAIL("stream denial fixture should permit display policy changes");
	denied_controller->setLodAutoSubmit(TRUE);
	const auto deadline = std::chrono::steady_clock::now() +
	    std::chrono::seconds(2);
	do {
	    (void)denied_controller->advanceProgressiveWork(&denied_options, NULL);
	    unsigned char *pixels = NULL;
	    const int rendered = denied_controller->renderToImage(&pixels, 0, 0,
		NULL, bobol_headless_context_manager());
	    const bool have_image = pixels != NULL;
	    if (pixels)
		bu_free(pixels, "stream denial preview frame");
	    if (rendered != BRLCAD_OK || !have_image)
		FAIL("source delivery failure should still permit rendering");
	    denied_controller->noteFramePresented();
	    denied_controller->getLodConvergenceStatus(denied_status);
	    if (denied_status.terminalError)
		break;
	    std::this_thread::sleep_for(std::chrono::milliseconds(1));
	} while (std::chrono::steady_clock::now() < deadline);
	if (!denied_status.terminalError || denied_status.viewReady ||
	    denied_status.failedSourceCount != 1 ||
	    denied_source->getCompactInstanceCount() != 1 ||
	    denied_source->realizationStatus.getValue() != SoBRLDatabaseSource::FAILED) {
	    fprintf(stderr, "stream denial display: automatic=%d terminal=%d error=%d "
		"pending=%d boxes=%zu failures=%zu work=%u facts=%u exact=%d executed=%d reason=%s\n", automatic,
		denied_status.terminal, denied_status.terminalError,
		denied_status.sourcePreparationPending,
		denied_status.presentedStructuralBoxCount,
		denied_status.terminalOccurrenceFailureCount,
		denied_status.controlObligationMask, denied_status.controlFactMask,
		denied_controller->getViewLodState()->lastCadPresentationFrameExact(),
		denied_controller->getViewLodState()->lastCadPresentationFrameExecuted(),
		denied_controller->getRenderReason().getString());
	    for (int i = 0; i < controller->getDatabaseSourceCount(); ++i) {
		SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
		fprintf(stderr, "source %s state=%d count=%d expected=%zu assembly=%p\n",
		    source->path.getValue().getString(), source->realizationStatus.getValue(),
		    source->getCompactInstanceCount(), source->getCompactExpectedInstanceCount(),
		    static_cast<void *>(denied_controller->getViewLodState()->findCadPresentation(source)));
	    }
	    FAIL("failed source delivery should settle honestly across display policies");
	}
	fprintf(stderr, "stream denial display settled: automatic=%d\n", automatic);
    }
    if (try_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, "progressive_root.c",
	    view_ctx, -1) <= 0)
	FAIL("stream denial should permit erasing its retained preview");
    std::vector<BObolLodControlTransitionRecord> policy_records;
    denied_controller->drainLodControlTransitions(policy_records);
    if (policy_records.empty() ||
	denied_controller->getDroppedLodControlTransitionCount())
	FAIL("policy/erase fixture should retain its transition evidence");
    for (const auto &record : policy_records) {
	if (record.event == BOBOL_LOD_CONTROL_TRANSITION_UNNAMED ||
	    record.before.convergence.controlViolationMask ||
	    record.after.convergence.controlViolationMask) {
	    fprintf(stderr, "policy/erase trace: serial=%llu event=%s violations=%u/%u\n",
		static_cast<unsigned long long>(record.serial),
		bobol_lod_control_transition_event_name(record.event),
		record.before.convergence.controlViolationMask,
		record.after.convergence.controlViolationMask);
	    FAIL("policy changes and policy-off erase need named, valid transitions");
	}
    }
    denied_controller->setLodControlTransitionTracing(FALSE);
    if (!ged_view_lod_policy_apply(view_ctx, &original_policy))
	FAIL("stream denial fixture should restore its display policy");

    struct TerminalDeliveryFailure {
	BObolTransactionFaultPoint point;
	int retainedCount;
	const char *diagnosticStage;
    };
    const TerminalDeliveryFailure deliveryFailures[] = {
	{BObolTransactionFaultPoint::SOURCE_STREAM_PARTIAL_MERGE_AFTER_COMPLETION, 2, "merge"},
	{BObolTransactionFaultPoint::SOURCE_TERMINAL_PREPARATION, 5, "adopt"}
    };
    for (const TerminalDeliveryFailure &failure : deliveryFailures) {
	ScopedTransactionFault fault(failure.point);
	ged_scene_reducer_result_init(&result);
	const int delivery_draw = ged_scene_reduce(gedp, &txn, &result);
	ged_scene_reducer_result_free(&result);
	SoBRLDatabaseSource *delivery_source =
	    source_for_path(controller, "progressive_root.c");
	if (delivery_draw <= 0 || !delivery_source)
	    FAIL("terminal delivery fixture should begin with a retained preview");
	BObolProgressiveOptions delivery_options;
	BObolLodConvergenceStatus delivery_status;
	const auto deadline = std::chrono::steady_clock::now() +
	    std::chrono::seconds(2);
	do {
	    (void)denied_controller->advanceProgressiveWork(&delivery_options, NULL);
	    denied_controller->getLodConvergenceStatus(delivery_status);
	    if (!delivery_status.sourcePreparationPending)
		break;
	    std::this_thread::sleep_for(std::chrono::milliseconds(1));
	} while (std::chrono::steady_clock::now() < deadline);
	if (delivery_status.sourcePreparationPending || delivery_status.failedSourceCount != 1 ||
	    delivery_status.viewReady ||
	    delivery_source->realizationStatus.getValue() != SoBRLDatabaseSource::FAILED ||
	    delivery_source->getCompactInstanceCount() != failure.retainedCount ||
	    !strstr(delivery_source->realizationDiagnostic.getValue().getString(), failure.diagnosticStage) ||
	    !delivery_source->hasExactSourceBounds())
	    FAIL("failed delivery should preserve its committed leaves and whole-root coverage");
	for (int i = 0; i < delivery_source->getCompactInstanceCount(); ++i) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    if (!delivery_source->getCompactInstanceHandle(i, handle) ||
		!delivery_source->getCompactInstanceSummary(handle, summary))
		FAIL("failed delivery must retain usable occurrence identities");
	}
	BObolViewLodState delivery_view;
	(void)delivery_source->compactViewLodAssembly(no_denied_payloads, &delivery_view);
	SoCADAssembly *delivery_presentation = delivery_view.findCadPresentation(delivery_source);
	if (!delivery_presentation || delivery_presentation->instanceCount() != size_t(failure.retainedCount))
	    FAIL("failed delivery should publish the same population to the retained renderer");
	if (try_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, "progressive_root.c",
		view_ctx, -1) <= 0)
	    FAIL("terminal delivery failure should permit erase");
    }

    ged_scene_reducer_result_init(&result);
    const int draw_ret = ged_scene_reduce(gedp, &txn, &result);
    ged_scene_reducer_result_free(&result);
    if (draw_ret <= 0)
	FAIL("deferred progressive root draw should succeed");

    SoBRLDatabaseSource *root_source =
	source_for_path(controller, "progressive_root.c");
    BObolViewController *view_controller = ged_bobol_view_controller(view_ctx);
    if (!root_source || !view_controller)
	FAIL("deferred progressive draw should retain its source and view controller");
    struct TerminalAdoptionObserver {
	enum { callbackWidth = 43 };
	BObolSceneController &scene;
	SoBRLDatabaseSource &source;
	SbString key;
	uint64_t frame;
	bool called = false;
	bool partial = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<TerminalAdoptionObserver *>(data);
	    if (self.called || self.source.realizationStatus.getValue() !=
		    SoBRLDatabaseSource::REALIZED ||
		self.source.getCompactInstanceCount() != 5)
		return;
	    self.called = true;
	    self.partial =
		self.scene.findDatabaseSourceInstance(self.key.getString()) !=
		    &self.source ||
		self.scene.getFrameRevision() <= self.frame ||
		self.source.needsRealization();
	    BObolDatabaseSourceDisplayPatch patch;
	    patch.lineWidthValid = TRUE;
	    patch.lineWidth = callbackWidth;
	    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(
		self.key.getString(), patch) <= 0;
	}
    } terminalObserver{*controller, *root_source,
	root_source->instanceKey.getValue(),
	controller->getFrameRevision()};
    SoNodeSensor terminalSensor(TerminalAdoptionObserver::changed,
	&terminalObserver);
    terminalSensor.setPriority(0);
    if (root_source)
	terminalSensor.attach(root_source);
    BObolProgressiveOptions progressive_options;
    BObolProgressiveStatus progressive_status;
    for (int attempt = 0;
	 view_controller && root_source &&
	 root_source->getCompactInstanceCount() != 4 && attempt < 2000;
	 attempt++) {
	(void)view_controller->advanceProgressiveWork(&progressive_options,
	    &progressive_status);
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    terminalSensor.detach();
    BObolDatabaseSourceSummary root_summary;
    const auto authoritative_occurrence_count =
	[](SoBRLDatabaseSource *source) {
	    size_t count = 0;
	    if (!source)
		return count;
	    for (int i = 0; i < source->getCompactInstanceCount(); ++i) {
		BObolCompactInstanceHandle handle;
		BObolCompactInstanceSummary summary;
		if (!source->getCompactInstanceHandle(i, handle) ||
		    !source->getCompactInstanceSummary(handle, summary) ||
		    !summary.valid || BU_STR_EQUAL(summary.geometryKind.getString(),
			"overview-aabb"))
		    continue;
		count++;
	    }
	    return count;
	};
    if (!root_source || !root_source->getSummary(root_summary) ||
	!root_summary.valid ||
	root_summary.realizationStatus != SoBRLDatabaseSource::REALIZED ||
	!root_source->isCompactOccurrenceRegistry() ||
	authoritative_occurrence_count(root_source) != 4 ||
	controller->getDatabaseSourceCount() != initial_source_count + 1)
	FAIL("deferred progressive draw should initially publish a compact occurrence registry");
    if (!terminalObserver.called || terminalObserver.partial ||
	terminalObserver.failed || root_source->lineWidth.getValue() !=
	    TerminalAdoptionObserver::callbackWidth) {
	fprintf(stderr, "terminal acceptance called=%d partial=%d failed=%d "
	    "width=%d count=%d status=%d frame=%llu/%llu current=%p source=%p\n",
	    terminalObserver.called, terminalObserver.partial,
	    terminalObserver.failed, root_source->lineWidth.getValue(),
	    root_source->getCompactInstanceCount(),
	    root_source->realizationStatus.getValue(),
	    static_cast<unsigned long long>(controller->getFrameRevision()),
	    static_cast<unsigned long long>(terminalObserver.frame),
	    static_cast<void *>(controller->findDatabaseSourceInstance(
		terminalObserver.key.getString())), static_cast<void *>(root_source));
	FAIL("terminal adoption callback should see complete state and retain its later edit");
    }
    BObolLodConvergenceStatus recovered_status;
    view_controller->getLodConvergenceStatus(recovered_status);
    if (recovered_status.failedSourceCount != 0)
	FAIL("fresh source realization should retire the preceding failure evidence");

    int duplicate_count = 0;
    int saw_first = 0;
    int saw_duplicate = 0;
    int saw_subtract = 0;
    int saw_intersect = 0;
    int saw_inherited_material = 0;
    SbString first_key;
	for (int i = 0; i < root_source->getCompactInstanceCount(); i++) {
	BObolCompactInstanceHandle handle;
	BObolCompactInstanceSummary summary;
	if (!root_source->getCompactInstanceHandle(i, handle) ||
	    !root_source->getCompactInstanceSummary(handle, summary) ||
	    !summary.valid)
	    continue;
	const char *summary_path = summary.path.getString();
	while (summary_path && *summary_path == '/')
	    summary_path++;
	if (BU_STR_EQUAL(summary_path,
		"progressive_root.c/dup_leaf.s") ||
	    BU_STR_EQUAL(summary_path,
		"progressive_root.c/dup_leaf.s@1")) {
	    duplicate_count++;
	    if (summary.occurrenceIndex == 0 &&
		fabs(summary.localToSource[3][0] - 11.0f) < 0.001f &&
		fabs(summary.localToSource[3][1] - 12.0f) < 0.001f &&
		fabs(summary.localToSource[3][2] - 13.0f) < 0.001f) {
		first_key = summary.sourceInstanceKey;
		saw_first = 1;
		if (summary.materialColorValid &&
		    fabs(summary.materialColor[0] - 42.0f / 255.0f) < 0.001f &&
		    fabs(summary.materialColor[1] - 84.0f / 255.0f) < 0.001f &&
		    fabs(summary.materialColor[2] - 126.0f / 255.0f) < 0.001f)
		    saw_inherited_material = 1;
	    }
	    if (summary.occurrenceIndex == 1 &&
		fabs(summary.localToSource[3][0] - 21.0f) < 0.001f &&
		fabs(summary.localToSource[3][1] - 22.0f) < 0.001f &&
		fabs(summary.localToSource[3][2] - 23.0f) < 0.001f &&
		summary.sourceInstanceKey != first_key)
		saw_duplicate = 1;
	} else if (BU_STR_EQUAL(summary_path,
		"progressive_root.c/box.s") &&
	    summary.booleanOperation == SoBRLDatabaseSource::BOOLEAN_SUBTRACT &&
	    summary.lineStyle == 1) {
	    saw_subtract = 1;
	} else if (BU_STR_EQUAL(summary_path,
		"progressive_root.c/ball.s") &&
	    summary.booleanOperation == SoBRLDatabaseSource::BOOLEAN_INTERSECT) {
	    saw_intersect = 1;
	}
    }
    if (duplicate_count != 2 || !saw_first || !saw_duplicate ||
	!saw_subtract || !saw_intersect || !saw_inherited_material) {
	fprintf(stderr,
	    "progressive identity duplicate=%d first=%d duplicate_key=%d subtract=%d intersect=%d material=%d\n",
	    duplicate_count, saw_first, saw_duplicate, saw_subtract,
	    saw_intersect, saw_inherited_material);
	for (int i = 0; i < root_source->getCompactInstanceCount(); i++) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    if (root_source->getCompactInstanceHandle(i, handle) &&
		root_source->getCompactInstanceSummary(handle, summary) &&
		summary.valid)
		fprintf(stderr,
		    "  [%d] path=%s occurrence=%u op=%d style=%d key=%s translation=(%.3f %.3f %.3f) material=%d color=(%.3f %.3f %.3f)\n",
		    i, summary.path.getString(), summary.occurrenceIndex,
		    summary.booleanOperation, summary.lineStyle,
		    summary.sourceInstanceKey.getString(),
		    summary.localToSource[3][0], summary.localToSource[3][1],
		    summary.localToSource[3][2], summary.materialColorValid,
		    summary.materialColor[0], summary.materialColor[1],
		    summary.materialColor[2]);
	}
	FAIL("compact proxy occurrences should preserve transforms, boolean identity, and inherited material");
    }

    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE,
	    "progressive_root.c", NULL, -1, "progressive identity cleanup"))
	return 1;
    if (controller->getDatabaseSourceCount() != initial_source_count)
	FAIL("progressive identity cleanup should restore source count");
    std::puts("PASS GED deferred-result acceptance: coverage and terminal callbacks remain current");
    return 0;
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
    unsigned long long pathHash;
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
    state->pathHash = record->path_hash;
    return 0;
}

static int
exercise_deferred_source_replacement(struct ged *gedp,
	BObolSceneController *scene)
{
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx || !ged_view_context_display_endpoint_ensure(view_ctx))
	FAIL("source replacement fixture needs a production display endpoint");
    scene = ged_draw_obol_scene_controller(gedp);
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!scene || !controller)
	FAIL("source replacement fixture needs an attached scene");
    const int initial_count = scene->getDatabaseSourceCount();
    enum class Replacement {
	Owner, Path, Population, Database, DatabaseRoundTrip, Representation,
	TessellationAbs, TessellationRel, TessellationNorm
    };
    std::unique_ptr<struct db_i, decltype(&db_close)> replacement_database(
	db_create_inmem(), db_close);
    if (!replacement_database)
	FAIL("source replacement fixture needs an independent database identity");
    for (const Replacement kind : {Replacement::Owner, Replacement::Path,
	    Replacement::Population, Replacement::Database,
	    Replacement::DatabaseRoundTrip, Replacement::Representation,
	    Replacement::TessellationAbs, Replacement::TessellationRel,
	    Replacement::TessellationNorm}) {
	const bool replace_owner = kind == Replacement::Owner;
	struct ged_draw_appearance_settings appearance =
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
	appearance.defer_leaf_expansion = 1;
	struct ged_scene_reducer_request txn =
	    ged_scene_reducer_request_make(GED_SCENE_REDUCER_DRAW, "progressive_root.c");
	txn.view = view_ctx;
	txn.mode = GED_DRAW_MODE_WIRE;
	txn.appearance = &appearance;
	if (ged_scene_reduce(gedp, &txn, NULL) <= 0)
	    FAIL("replacement fixture should launch its original source worker");
	SoBRLDatabaseSource *original = source_for_path(scene, "progressive_root.c");
	BObolDatabaseSourceSummary original_summary;
	if (!original || !original->getSummary(original_summary))
	    FAIL("replacement fixture should capture its original source stamp");
	const uint64_t original_route = original->getCompactSourceRoutingId();
	const std::string key = original_summary.instanceKey.getString();
	const char *replacement_path = kind == Replacement::Path ?
	    "ball.s" : "progressive_root.c";
	BObolDatabaseSourcePublishState replacement;
	replacement.sourceInstanceKey = key.c_str();
	replacement.sourcePath = replacement_path;
	replacement.sourceRepresentationKey = original_summary.representationKey.getString();
	replacement.database = gedp->dbip;
	replacement.drawMode = original_summary.drawMode;
	replacement.representationMode = original_summary.representationMode;
	replacement.sourceRevisionValid = TRUE;
	replacement.sourceRevision = original_summary.sourceRevision;
	replacement.inputsRevision = original_summary.inputsRevision;
	replacement.viewPolicyValid = TRUE;
	replacement.viewDependent = original->realizationViewDependent.getValue();
	replacement.csgLodEnabled = original->realizationCsgLodEnabled.getValue();
	replacement.meshLodEnabled = original->realizationMeshLodEnabled.getValue();
	replacement.viewScale = original->realizationViewScale.getValue();
	replacement.lodScale = original->realizationLodScale.getValue();
	replacement.viewWidth = original->realizationViewWidth.getValue();
	replacement.viewHeight = original->realizationViewHeight.getValue();
	replacement.botThreshold = original->realizationBotThreshold.getValue();
	replacement.curveScale = original->realizationCurveScale.getValue();
	replacement.pointScale = original->realizationPointScale.getValue();
	if (replace_owner && scene->removeDatabaseSourceInstance(key.c_str()) <= 0)
	    FAIL("replacement fixture should retire the original source owner");
	if (scene->publishDatabaseSourceInstance(replacement) <= 0)
	    FAIL("replacement fixture should publish its new source contract");
	SoBRLDatabaseSource *live = scene->findDatabaseSourceInstance(key.c_str());
	if (!live)
	    FAIL("replacement fixture should find the retained owner");
	const uint64_t original_population = live->getCompactPopulationEpoch();
	switch (kind) {
	    case Replacement::Owner:
	    case Replacement::Path:
	    case Replacement::Population:
		if (!live->realizeDatabaseWireframe())
		    FAIL("replacement fixture should install its own drawable population");
		break;
	    case Replacement::Database:
		live->setDatabase(replacement_database.get());
		break;
	    case Replacement::DatabaseRoundTrip:
		live->setDatabase(replacement_database.get());
		live->setDatabase(gedp->dbip);
		break;
	    case Replacement::Representation:
		(void)live->setRepresentationState("replacement-representation",
		    live->representationMode.getValue());
		break;
	    case Replacement::TessellationAbs:
		live->tessellationAbsTol = live->tessellationAbsTol.getValue() + 1.0f;
		live->markStale(SoBRLDatabaseSource::STALE_TESSELLATION);
		break;
	    case Replacement::TessellationRel:
		live->tessellationRelTol = live->tessellationRelTol.getValue() + 0.01f;
		live->markStale(SoBRLDatabaseSource::STALE_TESSELLATION);
		break;
	    case Replacement::TessellationNorm:
		live->tessellationNormTol = live->tessellationNormTol.getValue() + 0.1f;
		live->markStale(SoBRLDatabaseSource::STALE_TESSELLATION);
		break;
	}
	if (live->sourceRevision.getValue() != original_summary.sourceRevision ||
	    live->inputsRevision.getValue() != original_summary.inputsRevision ||
	    live->drawMode.getValue() != original_summary.drawMode ||
	    live->representationMode.getValue() != original_summary.representationMode ||
	    (live->getCompactSourceRoutingId() != original_route) != replace_owner)
	    FAIL("replacement fixture must preserve the old numeric admission stamp");
	const bool replaced_population = kind == Replacement::Owner ||
	    kind == Replacement::Path || kind == Replacement::Population;
	if ((live->getCompactPopulationEpoch() != original_population) != replaced_population)
	    FAIL("configuration-only replacement must retain its previous population");
	controller->notifyProgressiveSourceInputsChanged();
	const int replacement_count = live->getCompactInstanceCount();
	const uint64_t replacement_population = live->getCompactPopulationEpoch();
	const int replacement_status = live->realizationStatus.getValue();
	SbBox3f replacement_bounds;
	if (!replacement_count || !live->getSourceBounds(replacement_bounds))
	    FAIL("replacement fixture should realize its own drawable geometry");
	std::vector<BObolCompactOccurrence> retained(replacement_count);
	for (int i = 0; i < replacement_count; ++i)
	    if (!live->getCompactOccurrence(i, retained[i]))
		FAIL("replacement fixture should capture its retained occurrences");
	BObolLodConvergenceStatus convergence;
	BObolProgressiveOptions options;
	options.maxProviderItems = 1;
	{
	    /* The original small worker can complete with its output queued.
	     * A stale result must be rejected before even reaching leaf denial. */
	    ScopedTransactionFault fault(
		BObolTransactionFaultPoint::SOURCE_STREAM_MERGE_AFTER_COMPLETION);
	    const auto deadline = std::chrono::steady_clock::now() +
		std::chrono::seconds(2);
	    do {
		(void)controller->advanceProgressiveWork(&options, NULL);
		controller->getLodConvergenceStatus(convergence);
		SbBox3f current_bounds;
		if (live->getCompactInstanceCount() != replacement_count ||
		    live->getCompactPopulationEpoch() != replacement_population ||
		    live->realizationStatus.getValue() != replacement_status ||
		    !live->getSourceBounds(current_bounds) || current_bounds != replacement_bounds) {
		    fprintf(stderr, "source replacement: case=%d new-owner=%d path=%s count=%d/%d "
			"state=%d pending=%d\n", static_cast<int>(kind), replace_owner, replacement_path,
			live->getCompactInstanceCount(), replacement_count,
			live->realizationStatus.getValue(), convergence.sourcePreparationPending);
		    FAIL("stale streamed realization must not mutate the replacement source");
		}
		for (int i = 0; i < replacement_count; ++i) {
		    BObolCompactOccurrence current;
		    if (!live->getCompactOccurrence(i, current) ||
			current.geometry != retained[i].geometry ||
			current.geometryTransform != retained[i].geometryTransform ||
			current.localTransform != retained[i].localTransform ||
			current.summary.path != retained[i].summary.path)
			FAIL("stale delivery must preserve replacement geometry and placement");
		}
		if (!convergence.sourcePreparationPending)
		    break;
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
	    } while (std::chrono::steady_clock::now() < deadline);
	}
	if (convergence.sourcePreparationPending || convergence.failedSourceCount)
	    FAIL("superseded source work must retire without failing its replacement");
	if (scene->removeDatabaseSourceInstance(key.c_str()) <= 0)
	    FAIL("replacement fixture should remove its temporary source");
	controller->notifyProgressiveSourceInputsChanged();
	if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE,
		"progressive_root.c", view_ctx, -1, "source replacement cleanup") ||
	    scene->getDatabaseSourceCount() != initial_count)
	    FAIL("source replacement cleanup should restore the original scene");
    }
    return 0;
}

static int
exercise_progressive_autoview_lifecycle(struct ged *gedp,
	BObolViewController *controller, struct ged_view_context *view_ctx,
	int draw_mode, const char *draw_path, const char *renamed_path,
	size_t expected_authoritative_count, bool measure_worker_retirement)
{
    if (!gedp || !controller || !view_ctx || !draw_path || !draw_path[0] ||
	!renamed_path || !renamed_path[0] || !expected_authoritative_count)
	FAIL("progressive autoview test needs an attached view");

    constexpr int lifecycle_width = 160;
    constexpr int lifecycle_height = 120;
    constexpr size_t lifecycle_pixel_bytes =
	static_cast<size_t>(lifecycle_width) * lifecycle_height * 3;
    controller->setViewportSize(lifecycle_width, lifecycle_height);
    if (!controller->syncCameraFromViewContext(view_ctx))
	FAIL("progressive lifecycle should establish its initial camera");

    struct bv *view = DRAW_TEST_BV(view_ctx);
    const uint64_t initial_revision = bv_frame_revision_get(view);
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    const int initial_scene_source_count = scene ?
	scene->getDatabaseSourceCount() : 0;

    struct ged_draw_appearance_settings appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    appearance.draw_mode = draw_mode;
    appearance.defer_leaf_expansion = 1;
    struct ged_scene_reducer_request txn =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_DRAW, draw_path);
    txn.view = view_ctx;
    txn.mode = draw_mode;
    txn.appearance = &appearance;
    txn.autoview = 1;
    struct ged_scene_reducer_result result;
    ged_scene_reducer_result_init(&result);
    const int draw_ret = ged_scene_reduce(gedp, &txn, &result);
    ged_scene_reducer_result_free(&result);
    if (draw_ret <= 0)
	FAIL("progressive autoview deferred draw should succeed");

    /* An explicit autoview while the root is still realizing must replace
     * the transaction's initial fit and continue following source-bound
     * growth until the final extent is certified. */
    const char *autoview_cmd[1] = {"autoview"};
    if (ged_exec_autoview(gedp, 1, autoview_cmd) != BRLCAD_OK)
	FAIL("explicit autoview should arm deferred Obol bound tracking");

    uint64_t observed_autoview_revision = initial_revision;
    size_t autoview_application_ticks = 0;
    const auto note_autoview_application = [&]() {
	const uint64_t revision = bv_frame_revision_get(view);
	if (revision == observed_autoview_revision)
	    return;
	observed_autoview_revision = revision;
	autoview_application_ticks++;
    };
    note_autoview_application();

    BObolProgressiveOptions options;
    BObolProgressiveStatus status;
    auto compact_counts = [](SoBRLDatabaseSource *source,
	    size_t &authoritative, size_t &overviews,
	    size_t &visible_overviews) {
	authoritative = 0;
	overviews = 0;
	visible_overviews = 0;
	if (!source)
	    return;
	for (int i = 0; i < source->getCompactInstanceCount(); ++i) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    if (!source->getCompactInstanceHandle(i, handle) ||
		!source->getCompactInstanceSummary(handle, summary) ||
		!summary.valid)
		continue;
	    if (BU_STR_EQUAL(summary.geometryKind.getString(),
		    "overview-aabb")) {
		overviews++;
		if (summary.visible)
		    visible_overviews++;
	    } else {
		authoritative++;
	    }
	}
    };
    SoBRLDatabaseSource *initial_source = scene ? source_for_path(scene,
	draw_path) : NULL;
    SbBox3f proxy_bounds;
    proxy_bounds.makeEmpty();
    size_t authoritative_count = 0;
    size_t overview_count = 0;
    size_t visible_overview_count = 0;
    for (int attempt = 0; attempt < 2000; attempt++) {
	compact_counts(initial_source, authoritative_count, overview_count,
	    visible_overview_count);
	const bool proxy_ready =
	    initial_source &&
	    initial_source->isCompactOccurrenceRegistry() &&
	    authoritative_count == expected_authoritative_count &&
	    initial_source->getEffectiveSourceBounds(proxy_bounds) &&
	    !proxy_bounds.isEmpty();
	if (proxy_ready)
	    break;
	(void)controller->advanceProgressiveWork(&options, &status);
	note_autoview_application();
	initial_source = scene ? source_for_path(scene,
	    draw_path) : NULL;
	proxy_bounds.makeEmpty();
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    compact_counts(initial_source, authoritative_count, overview_count,
	visible_overview_count);
    if (!initial_source || !initial_source->isCompactOccurrenceRegistry() ||
	authoritative_count != expected_authoritative_count ||
	!initial_source->getEffectiveSourceBounds(proxy_bounds) ||
	proxy_bounds.isEmpty()) {
	fprintf(stderr,
	    "progressive proxy source=%p compact=%d count=%d authoritative=%zu "
	    "overviews=%zu visible_overviews=%zu bounds=%d "
	    "empty=%d more=%d changed=%d providers=%zu advanced=%zu "
	    "remaining=%zu pending=%zu\n",
	    (void *)initial_source,
	    initial_source ?
		initial_source->isCompactOccurrenceRegistry() : -1,
	    initial_source ? initial_source->getCompactInstanceCount() : -1,
	    authoritative_count, overview_count, visible_overview_count,
	    initial_source ?
		initial_source->getEffectiveSourceBounds(proxy_bounds) : 0,
	    proxy_bounds.isEmpty() ? 1 : 0,
	    status.hasMore, status.changed, status.providerCount,
	    status.providerAdvanced, status.remaining, status.pendingTasks);
	FAIL("deferred root should publish compact conservative proxy bounds");
    }

    /* Initial AABBs and final detail share the same root source. */
    if (scene->getDatabaseSourceCount() != initial_scene_source_count + 1)
	FAIL("deferred root should not create proxy child sources");

    /* A redraw can advance the aggregate draw revision while this source is
     * realizing.  It must not invalidate an otherwise matching snapshot. */
    const uint64_t revision_before_redraw = ged_draw_scene_revision(gedp);
    struct ged_scene_reducer_request redraw_txn =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_REDRAW, draw_path);
    redraw_txn.view = view_ctx;
    if (ged_scene_reduce(gedp, &redraw_txn, NULL) < 0 ||
	ged_draw_scene_revision(gedp) <= revision_before_redraw)
	FAIL("progressive redraw should advance scene bookkeeping without cancelling refinement");

    int initial_progress = 0;
    int settled = 0;
    SoBRLDatabaseSource *settled_source = NULL;
    for (int attempt = 0; attempt < 2000; attempt++) {
	initial_progress = controller->advanceProgressiveWork(&options, &status);
	note_autoview_application();
	settled_source = scene ? source_for_path(scene,
	    draw_path) : NULL;
	compact_counts(settled_source, authoritative_count, overview_count,
	    visible_overview_count);
	if (settled_source && settled_source->isCompactOccurrenceRegistry() &&
	    authoritative_count == expected_authoritative_count &&
	    visible_overview_count == 0 &&
	    !status.hasMore) {
	    settled = 1;
	    break;
	}
	if (!status.hasMore)
	    break;
	/* Production hosts present every published PoP cut before admitting the
	 * next one.  This headless contract test has no paint loop, so explicitly
	 * acknowledge one requested frame rather than spinning forever behind
	 * the controller's completed-frame gate. */
	SbString render_reason;
	if (controller->consumeRenderRequest(&render_reason)) {
	    const uint64_t frame_started = controller->beginRenderTiming();
	    std::this_thread::sleep_for(std::chrono::microseconds(1));
	    controller->completeRenderTiming(
		frame_started, test_capacity_cad_timing());
	}
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    SbBox3f settled_bounds;
    settled_bounds.makeEmpty();
    if (settled_source)
	(void)settled_source->getEffectiveSourceBounds(settled_bounds);
    compact_counts(settled_source, authoritative_count, overview_count,
	visible_overview_count);
    /* A fast worker or warm cache may finish before the first explicit pump.
     * The contract is at least one useful frame and the stable final autoview;
     * a cold producer is permitted to apply several throttled partial fits. */
    if (!settled || !settled_source ||
	!settled_source->isCompactOccurrenceRegistry() ||
	authoritative_count != expected_authoritative_count ||
	visible_overview_count != 0 ||
	settled_bounds.isEmpty() ||
	scene->getDatabaseSourceCount() != initial_scene_source_count + 1 ||
	bv_frame_revision_get(view) <= initial_revision ||
	autoview_application_ticks < 1) {
	fprintf(stderr, "progressive settle ret=%d changed=%d settled=%d "
	    "bounds_empty=%d frame=%llu initial=%llu providers=%zu "
	    "advanced=%zu remaining=%zu pending=%zu scene=%d/%d source=%p "
	    "compact=%d count=%d authoritative=%zu overviews=%zu "
	    "visible_overviews=%zu autoview_ticks=%zu\n",
	    initial_progress, status.changed,
	    settled, settled_bounds.isEmpty(),
	    static_cast<unsigned long long>(bv_frame_revision_get(view)),
	    static_cast<unsigned long long>(initial_revision),
	    status.providerCount,
	    status.providerAdvanced, status.remaining, status.pendingTasks,
	    scene ? scene->getDatabaseSourceCount() : -1,
	    initial_scene_source_count, (void *)settled_source,
	    settled_source ? settled_source->isCompactOccurrenceRegistry() : -1,
	    settled_source ? settled_source->getCompactInstanceCount() : -1,
	    authoritative_count, overview_count, visible_overview_count,
	    autoview_application_ticks);
	if (scene) {
	    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
		BObolDatabaseSourceSummary summary;
		SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
		if (scene->getDatabaseSourceSummary(i, summary) && summary.valid)
		    fprintf(stderr, "  source[%d] key=%s path=%s rep=%d compact=%d count=%d\n",
			i, summary.instanceKey.getString(), summary.path.getString(),
			summary.representationMode,
			source ? source->isCompactOccurrenceRegistry() : -1,
			source ? source->getCompactInstanceCount() : -1);
	    }
	}
	FAIL("background progressive refinement should settle compact detail without retaining structural proxy geometry");
    }

    uint64_t stable_autoview_revision = bv_frame_revision_get(view);
    int stable_autoview_ticks = 0;
    for (int attempt = 0; attempt < 64 && stable_autoview_ticks < 2;
	attempt++) {
	(void)controller->advanceProgressiveWork(&options, &status);
	SbString render_reason;
	if (status.hasMore &&
	    controller->consumeRenderRequest(&render_reason)) {
	    const uint64_t frame_started = controller->beginRenderTiming();
	    std::this_thread::sleep_for(std::chrono::microseconds(1));
	    controller->completeRenderTiming(
		frame_started, test_capacity_cad_timing());
	}
	const uint64_t current_revision = bv_frame_revision_get(view);
	if (current_revision == stable_autoview_revision) {
	    stable_autoview_ticks++;
	} else {
	    stable_autoview_revision = current_revision;
	    stable_autoview_ticks = 0;
	}
	if (stable_autoview_ticks < 2)
	    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (stable_autoview_ticks < 2)
	FAIL("settled progressive autoview should stop changing the view");

    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE,
	    draw_path, view_ctx, -1,
	    "progressive autoview settle cleanup"))
	return 1;

    ged_scene_reducer_result_init(&result);
    const int redraw_ret = ged_scene_reduce(gedp, &txn, &result);
    ged_scene_reducer_result_free(&result);
    if (redraw_ret <= 0)
	FAIL("progressive autoview cancellation draw should succeed");
    bv_size_set(view, 1234.0);
    const fastf_t user_size = bv_size_get(view);
    const char *rename_progressive[3] = {"move", draw_path, renamed_path};
    if (ged_exec(gedp, 3, rename_progressive) != BRLCAD_OK)
	FAIL("active background refinement database rename should succeed");
    for (int attempt = 0; attempt < 20; attempt++) {
	(void)controller->advanceProgressiveWork(&options, &status);
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (scene && (source_for_path(scene, draw_path) ||
	!source_for_path(scene, renamed_path)))
	FAIL("database rename should cancel stale refinement and retarget only the live proxy");
    const char *restore_progressive[3] = {"move", renamed_path, draw_path};
    if (ged_exec(gedp, 3, restore_progressive) != BRLCAD_OK)
	FAIL("active background refinement database rename restore should succeed");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE,
	    draw_path, view_ctx, -1,
	    "progressive autoview cancellation cleanup"))
	return 1;
    for (int attempt = 0; attempt < 20; attempt++) {
	(void)controller->advanceProgressiveWork(&options, &status);
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!NEAR_EQUAL(bv_size_get(view), user_size, SMALL_FASTF))
	FAIL("user view change should cancel pending progressive autoview");
    if (scene && source_for_path(scene, draw_path))
	FAIL("cancelled background refinement must not republish an erased root");

    /* Closing a production endpoint owns cancellation and worker retirement.
     * Recreate the endpoint against the same GED scene, redraw from the cold
     * fixture, and require the replacement controller to reach an observable
     * terminal frame with no queued service work. */
    const size_t closing_worker_count = controller->getManagedLodWorkerCount();
#if defined(__linux__)
    const size_t closing_thread_count = measure_worker_retirement ?
	bu_file_list("/proc/self/task", "[0-9]*", NULL) : 0;
#endif
    if (!ged_view_context_obol_endpoint_set(view_ctx, NULL, 0))
	FAIL("worker-active production endpoint should close cleanly");
#if defined(__linux__)
    bool workers_retired = !measure_worker_retirement ||
	closing_worker_count == 0;
    const auto worker_deadline = std::chrono::steady_clock::now() +
	std::chrono::seconds(2);
    while (!workers_retired &&
	std::chrono::steady_clock::now() < worker_deadline) {
	const size_t current_threads =
	    bu_file_list("/proc/self/task", "[0-9]*", NULL);
	workers_retired = current_threads + closing_worker_count <=
	    closing_thread_count;
	if (!workers_retired)
	    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!workers_retired)
	FAIL("closing a production endpoint should retire its managed LoD workers");
#endif

    if (!ged_view_context_display_endpoint_ensure(view_ctx))
	FAIL("closed production view should recreate its display endpoint");
    controller = ged_bobol_view_controller(view_ctx);
    scene = ged_draw_obol_scene_controller(gedp);
    if (!controller || !scene ||
	scene->getDatabaseSourceCount() != initial_scene_source_count)
	FAIL("reopened production endpoint should preserve the current GED scene");
    controller->setViewportSize(lifecycle_width, lifecycle_height);
    if (!controller->syncCameraFromViewContext(view_ctx))
	FAIL("reopened production endpoint should accept the retained camera");

    unsigned char *baseline_pixels = NULL;
    if (controller->renderToImage(&baseline_pixels, 0, 0, NULL,
	    bobol_headless_context_manager()) != BRLCAD_OK || !baseline_pixels) {
	if (baseline_pixels)
	    bu_free(baseline_pixels, "production lifecycle baseline frame");
	FAIL("reopened production endpoint should render its retained baseline");
    }
    std::vector<unsigned char> baseline(
	baseline_pixels, baseline_pixels + lifecycle_pixel_bytes);
    bu_free(baseline_pixels, "production lifecycle baseline frame");

    ged_scene_reducer_result_init(&result);
    const int reopened_draw = ged_scene_reduce(gedp, &txn, &result);
    ged_scene_reducer_result_free(&result);
    if (reopened_draw <= 0)
	FAIL("reopened production endpoint should accept a deferred redraw");

    BObolProgressiveOptions terminal_options;
    terminal_options.forceTerminalLodRefinement = TRUE;
    BObolProgressiveStatus terminal_progress;
    BObolLodConvergenceStatus terminal_status;
    std::vector<unsigned char> terminal_frame;
    const auto terminal_deadline = std::chrono::steady_clock::now() +
	std::chrono::seconds(8);
    do {
	(void)controller->advanceProgressiveWork(&terminal_options,
	    &terminal_progress);
	unsigned char *pixels = NULL;
	const int rendered = controller->renderToImage(&pixels, 0, 0, NULL,
	    bobol_headless_context_manager(), &terminal_progress);
	if (rendered != BRLCAD_OK || !pixels) {
	    if (pixels)
		bu_free(pixels, "production lifecycle terminal frame");
	    FAIL("reopened production endpoint should render progressive frames");
	}
	terminal_frame.assign(pixels, pixels + lifecycle_pixel_bytes);
	bu_free(pixels, "production lifecycle terminal frame");
	controller->noteFramePresented();
	controller->getLodConvergenceStatus(terminal_status);
	if (terminal_status.terminal)
	    break;
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    } while (std::chrono::steady_clock::now() < terminal_deadline);

    const bool image_changed = terminal_frame.size() == baseline.size() &&
	!std::equal(terminal_frame.begin(), terminal_frame.end(),
	    baseline.begin());
    SoBRLDatabaseSource *reopened_source =
	source_for_path(scene, draw_path);
    BObolLodService *terminal_service = controller->getLodService();
    const bool service_idle = terminal_service &&
	terminal_service->workStatus().isIdle();
    if (!terminal_status.terminal || !terminal_status.viewReady ||
	terminal_status.terminalError || !image_changed || !reopened_source ||
	!reopened_source->isCompactOccurrenceRegistry() || !service_idle) {
	fprintf(stderr,
	    "production lifecycle terminal=%d ready=%d error=%d image=%d "
	    "source=%p compact=%d workers=%zu service=%p "
	    "pending=%zu executing=%zu in_flight=%zu reservations=%zu active=%zu "
	    "results=%zu cache_writes=%zu delayed=%zu\n",
	    terminal_status.terminal, terminal_status.viewReady,
	    terminal_status.terminalError, image_changed,
	    static_cast<void *>(reopened_source),
	    reopened_source ? reopened_source->isCompactOccurrenceRegistry() : -1,
	    closing_worker_count, static_cast<void *>(terminal_service),
	    terminal_service ? terminal_service->pendingTaskCountForDiagnostics() : 0,
	    terminal_service ? terminal_service->executingTaskCountForDiagnostics() : 0,
	    terminal_service ? terminal_service->inFlightCount() : 0,
	    terminal_service ? terminal_service->resultReservationCountForDiagnostics() : 0,
	    terminal_service ? terminal_service->activeRequestCountForDiagnostics() : 0,
	    terminal_service ? terminal_service->queuedResultCountForDiagnostics() : 0,
	    terminal_service ? terminal_service->queuedCacheWriteCountForDiagnostics() : 0,
	    terminal_service ? terminal_service->delayedTaskCountForDiagnostics() : 0);
	FAIL("reopened production endpoint should publish one terminal image and release transient work");
    }

    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE,
	    draw_path, view_ctx, -1,
	    "production lifecycle terminal cleanup"))
	return 1;
    if (scene->getDatabaseSourceCount() != initial_scene_source_count)
	FAIL("production lifecycle cleanup should restore its starting scene");
    std::puts("PASS GED production flow lifecycle: cold draw, camera cancellation, close, reopen, terminal image, resource release");
    return 0;
}

static int
exercise_delayed_mesh_camera_close(struct ged *gedp,
	struct ged_view_context *view_ctx)
{
    if (!gedp || !view_ctx ||
	!ged_view_context_display_endpoint_ensure(view_ctx))
	FAIL("delayed mesh lifecycle needs a production display endpoint");
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!controller || !scene)
	FAIL("delayed mesh lifecycle needs attached production owners");

    constexpr const char *mesh_path = "mesh_owner.bot";
    constexpr int viewport_width = 160;
    constexpr int viewport_height = 120;
    constexpr fastf_t changed_camera_size = 641.0;
    const int initial_source_count = scene->getDatabaseSourceCount();
    controller->setViewportSize(viewport_width, viewport_height);
    if (!controller->syncCameraFromViewContext(view_ctx))
	FAIL("delayed mesh lifecycle should establish its camera");
    controller->setLodAutoSubmit(TRUE);

    const char *draw_mesh[2] = {"draw", mesh_path};
    if (ged_exec_draw(gedp, 2, draw_mesh) != BRLCAD_OK)
	FAIL("delayed mesh lifecycle cold draw should succeed");
    record_source_state draw_record = {};
    draw_record.matchPath = mesh_path;
    ged_draw_foreach_shape_record(gedp, record_source_state_cb, &draw_record);
    struct ged_view_context *lod_views[1] = {view_ctx};
    if (!draw_record.found ||
	!ged_draw_shape_ref_lod_ensure(gedp, draw_record.ref,
	    view_ctx, lod_views, 1))
	FAIL("delayed mesh lifecycle should start its source mesh provider");

    unsigned char *first_pixels = NULL;
    if (controller->renderToImage(&first_pixels, 0, 0, NULL,
	    bobol_headless_context_manager()) != BRLCAD_OK || !first_pixels) {
	if (first_pixels)
	    bu_free(first_pixels, "delayed mesh first frame");
	FAIL("delayed mesh lifecycle should present its first frame");
    }
    bu_free(first_pixels, "delayed mesh first frame");
    controller->noteFramePresented();

    BObolProgressiveOptions options;
    BObolProgressiveStatus progress;
    BObolLodService *service = controller->getLodService();
    bool observed_delayed_task = false;
    const auto delay_deadline = std::chrono::steady_clock::now() +
	std::chrono::seconds(4);
    do {
	(void)controller->advanceProgressiveWork(&options, &progress);
	service = controller->getLodService();
	if (service && service->delayedTaskCountForDiagnostics() > 0) {
	    observed_delayed_task = true;
	    break;
	}
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    } while (std::chrono::steady_clock::now() < delay_deadline);
    SoBRLDatabaseSource *mesh_source = source_for_path(scene, mesh_path);
    if (!observed_delayed_task || !service || !mesh_source ||
	!mesh_source->getMeshLod())
	FAIL("delayed mesh lifecycle should expose active delayed production work");

    struct bv *view = DRAW_TEST_BV(view_ctx);
    const uint64_t camera_revision = bv_frame_revision_get(view);
    bv_size_set(view, changed_camera_size);
    if (bv_frame_revision_get(view) <= camera_revision ||
	!controller->syncCameraFromViewContext(view_ctx) ||
	!NEAR_EQUAL(bv_size_get(view), changed_camera_size, SMALL_FASTF))
	FAIL("delayed mesh lifecycle should accept camera input during loading");

    const size_t closing_worker_count = controller->getManagedLodWorkerCount();
#if defined(__linux__)
    const size_t closing_thread_count =
	bu_file_list("/proc/self/task", "[0-9]*", NULL);
#endif
    if (!ged_view_context_obol_endpoint_set(view_ctx, NULL, 0))
	FAIL("delayed mesh lifecycle should close its worker-active endpoint");
#if defined(__linux__)
    bool workers_retired = closing_worker_count == 0;
    const auto worker_deadline = std::chrono::steady_clock::now() +
	std::chrono::seconds(2);
    while (!workers_retired &&
	std::chrono::steady_clock::now() < worker_deadline) {
	const size_t current_threads =
	    bu_file_list("/proc/self/task", "[0-9]*", NULL);
	workers_retired = current_threads + closing_worker_count <=
	    closing_thread_count;
	if (!workers_retired)
	    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!workers_retired)
	FAIL("delayed mesh endpoint close should retire managed workers");
#endif

    if (!ged_view_context_display_endpoint_ensure(view_ctx))
	FAIL("delayed mesh lifecycle should reopen its display endpoint");
    controller = ged_bobol_view_controller(view_ctx);
    scene = ged_draw_obol_scene_controller(gedp);
    if (!controller || !scene || !source_for_path(scene, mesh_path))
	FAIL("delayed mesh reopen should retain the current GED scene");
    const char *erase_mesh[2] = {"erase", mesh_path};
    if (ged_exec_erase(gedp, 2, erase_mesh) != BRLCAD_OK ||
	scene->getDatabaseSourceCount() != initial_source_count ||
	source_for_path(scene, mesh_path))
	FAIL("delayed mesh cancellation should not republish its erased source");

    std::puts("PASS GED delayed mesh lifecycle: camera input and worker-active endpoint close");
    return 0;
}

static SoBRLVListShape *
auxiliary_for_path_variant(SoBRLDatabaseSource *source, const char *path)
{
    if (!source || !path || !path[0])
	return NULL;

    SoBRLVListShape *shape = source->findAuxiliaryVListShape(path);
    if (shape)
	return shape;

    if (path[0] != '/') {
	std::string slash_path = "/";
	slash_path += path;
	shape = source->findAuxiliaryVListShape(slash_path.c_str());
	if (shape)
	    return shape;
    } else {
	shape = source->findAuxiliaryVListShape(skip_leading_slash(path));
	if (shape)
	    return shape;
    }

    for (int i = 0; i < source->getNumChildren(); i++) {
	SoNode *child = source->getChild(i);
	if (!child || !child->isOfType(SoBRLVListShape::getClassTypeId()))
	    continue;
	SoBRLVListShape *candidate = static_cast<SoBRLVListShape *>(child);
	if (bu_strcmp(candidate->recordRole.getValue().getString(),
		"auxiliary") == 0)
	    return candidate;
    }

    return NULL;
}

static int
exercise_typed_pick_result(void)
{
    struct ged_pick_result *result = ged_pick_result_create();
    if (!result)
	FAIL("typed pick result should allocate");

    struct ged_pick_detail input = ged_pick_detail_default();
    input.source_id = 17;
    input.primitive_kind = 3;
    input.primitive_index = 9;
    input.material_id = 23;
    input.face_vertex_index[0] = 4;
    input.face_vertex_index[1] = 5;
    input.face_vertex_index[2] = 6;
    input.nearest_face_vertex_index = 5;
    VSET(input.model_point, 1.0, 2.0, 3.0);
    input.model_point_valid = 1;
    if (!ged_pick_result_append_detail(result, "bot.s", 4.5,
	    &input)) {
	ged_pick_result_free(result);
	FAIL("typed pick result should append a detailed hit");
    }

    struct ged_pick_result *first =
	ged_pick_result_filter_first(result);
    struct ged_pick_detail output = ged_pick_detail_default();
    if (!first || ged_pick_result_count(first) != 1 ||
	!ged_pick_result_detail(first, 0, &output) ||
	output.source_id != input.source_id ||
	output.primitive_kind != input.primitive_kind ||
	output.primitive_index != input.primitive_index ||
	output.material_id != input.material_id ||
	output.face_vertex_index[0] != input.face_vertex_index[0] ||
	output.face_vertex_index[1] != input.face_vertex_index[1] ||
	output.face_vertex_index[2] != input.face_vertex_index[2] ||
	output.nearest_face_vertex_index !=
	    input.nearest_face_vertex_index ||
	!output.model_point_valid ||
	!VNEAR_EQUAL(output.model_point, input.model_point, SMALL_FASTF)) {
	ged_pick_result_free(first);
	ged_pick_result_free(result);
	FAIL("typed pick filtering should preserve primitive edit detail");
    }

    ged_pick_result_free(first);
    ged_pick_result_free(result);
    return 0;
}

struct record_source_mode_state {
    record_source_state recordState;
    int matchDrawMode;
};

static int
record_source_mode_state_cb(const struct ged_draw_shape_record *record,
	void *userdata)
{
    record_source_mode_state *state =
	static_cast<record_source_mode_state *>(userdata);
    if (!state || !record)
	return 1;
    if (record->draw_mode != state->matchDrawMode)
	return 1;
    return record_source_state_cb(record, &state->recordState);
}

static int
exercise_evaluated_wire_shape_ref_realize_context(struct ged *gedp,
	BObolSceneController *controller,
	const char *path)
{
    if (!gedp || !controller || !path)
	FAIL("evaluated-wire realize-context test needs GED and Obol scene state");

    const int mode = GED_DRAW_MODE_EVAL_WIRE;
    const int representation = SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE;
    (void)try_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, path,
	    ged_draw_active_view_ctx(gedp), mode);

    char mode_arg[16] = {0};
    snprintf(mode_arg, sizeof(mode_arg), "-m%d", mode);
    const char *draw_mode_cmd[4] = {"draw", mode_arg, path, NULL};
    if (ged_exec_draw(gedp, 3, draw_mode_cmd) != BRLCAD_OK)
	FAIL("evaluated-wire realize-context setup draw should succeed");

    SoBRLDatabaseSource *source =
	source_for_representation(controller, path, representation);
    BObolDatabaseSourceSummary summary;
    if (!source || !source->getSummary(summary) || !summary.valid ||
	    summary.realizationStatus != SoBRLDatabaseSource::REALIZED ||
	    summary.realizedShapeCount <= 0 ||
	    summary.instanceKey.getLength() == 0 ||
	    !summary.sourceBoundsValid || !summary.sourceBoundsExact)
	FAIL("evaluated-wire realize-context setup should create a realized mode source");

    const std::string instance_key(summary.instanceKey.getString());
    if (controller->markDatabaseSourceInstanceStale(instance_key.c_str(),
	    SoBRLDatabaseSource::STALE_SOURCE) < 0)
	FAIL("evaluated-wire realize-context setup should mark the mode source stale");
    if (!source->getSummary(summary) || !summary.valid ||
	    !summary.stale ||
	    summary.realizationStatus != SoBRLDatabaseSource::UNREALIZED)
	FAIL("evaluated-wire source should be stale before realize-context");

    record_source_mode_state eval_state = {};
    eval_state.recordState.ref = GED_DRAW_SHAPE_REF_NULL;
    eval_state.recordState.group = GED_DRAW_GROUP_REF_NULL;
    eval_state.recordState.matchPath = path;
    eval_state.matchDrawMode = mode;
    ged_draw_foreach_shape_record(gedp, record_source_mode_state_cb,
	    &eval_state);
    if (!eval_state.recordState.found ||
	    ged_draw_shape_ref_is_null(eval_state.recordState.ref))
	FAIL("evaluated-wire shape record should produce a mode-specific ref");

    struct ged_draw_shape_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_draw_shape_record_get(gedp, eval_state.recordState.ref,
	    &record) || record.draw_mode != mode)
	FAIL("evaluated-wire shape ref should retain its mode identity");

    if (!ged_draw_obol_database_source_realize_for_path(gedp, path))
	FAIL("evaluated-wire Obol source realization should succeed");
    if (!source->getSummary(summary) || !summary.valid ||
	    summary.stale ||
	    summary.realizationStatus != SoBRLDatabaseSource::REALIZED ||
	    summary.realizedShapeCount <= 0 ||
	    !summary.sourceBoundsValid || !summary.sourceBoundsExact)
	FAIL("evaluated-wire shape-ref realize-context should refresh the mode source");

    (void)try_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, path,
	    ged_draw_active_view_ctx(gedp), mode);
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

struct shape_index_state {
    int count;
    ged_draw_shape_ref ref;
};

static int
shape_index_state_cb(ged_draw_shape_ref ref, void *userdata)
{
    shape_index_state *state = static_cast<shape_index_state *>(userdata);
    if (!state || ged_draw_shape_ref_is_null(ref))
	return 1;

    state->count++;
    if (ged_draw_shape_ref_is_null(state->ref))
	state->ref = ref;
    return 1;
}

struct group_index_state {
    int count;
    ged_draw_group_ref ref;
};

static int
group_index_state_cb(ged_draw_group_ref ref, void *userdata)
{
    group_index_state *state = static_cast<group_index_state *>(userdata);
    if (!state || ged_draw_group_ref_is_null(ref))
	return 1;

    state->count++;
    if (ged_draw_group_ref_is_null(state->ref))
	state->ref = ref;
    return 1;
}


static std::set<std::string>
frontier_paths(struct ged *gedp)
{
    std::set<std::string> paths;
    struct bu_vls listing = BU_VLS_INIT_ZERO;
    (void)ged_scene_paths_append(gedp, ged_view_active_ctx(gedp),
	GED_SCENE_DRAW_DEFAULT, GED_SCENE_PATHS_DRAW_INTENTS, &listing);
    const char *start = bu_vls_cstr(&listing);
    for (const char *p = start; ; p++) {
	if (*p != '\n' && *p != '\r' && *p != '\0')
	    continue;
	if (p > start)
	    paths.insert(std::string(start, static_cast<size_t>(p - start)));
	if (*p == '\0')
	    break;
	start = p + 1;
    }
    bu_vls_free(&listing);
    return paths;
}


static int
frontier_draw_one(struct ged *gedp, const char *path)
{
    const char *args[2] = {"draw", path};
    return ged_exec_draw(gedp, 2, args);
}


static int
frontier_draw_one_mode(struct ged *gedp, int mode, const char *path)
{
    struct bu_vls mode_arg = BU_VLS_INIT_ZERO;
    bu_vls_printf(&mode_arg, "-m%d", mode);
    const char *args[3] = {"draw", bu_vls_cstr(&mode_arg), path};
    const int ret = ged_exec_draw(gedp, 3, args);
    bu_vls_free(&mode_arg);
    return ret;
}


static void
frontier_clear(struct ged *gedp)
{
    struct ged_scene_reducer_request clear =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_CLEAR, NULL);
    (void)ged_scene_reduce(gedp, &clear, NULL);
}


static int
frontier_visible_occurrences(SoBRLDatabaseSource *source)
{
    if (!source)
	return -1;
    int visible = 0;
    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactInstanceHandle handle;
	BObolCompactInstanceSummary summary;
	if (source->getCompactInstanceHandle(i, handle) &&
	    source->getCompactInstanceSummary(handle, summary) &&
	    summary.valid && summary.visible)
	    visible++;
    }
    return visible;
}


static int
frontier_selected_occurrences_for_path(SoBRLDatabaseSource *source,
	const char *path)
{
    if (!source || !path)
	return -1;
    int selected = 0;
    const auto semantic_path = [](const char *input) {
	std::string value;
	if (!input)
	    return value;
	while (*input == '/')
	    input++;
	for (const char *cursor = input; *cursor; cursor++) {
	    if (*cursor == '@' && cursor[1] >= '0' && cursor[1] <= '9') {
		while (cursor[1] && cursor[1] != '/')
		    cursor++;
		continue;
	    }
	    value.push_back(*cursor);
	}
	return value;
    };
    const std::string target_path = semantic_path(path);
    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactInstanceHandle handle;
	BObolCompactInstanceSummary summary;
	if (!source->getCompactInstanceHandle(i, handle) ||
	    !source->getCompactInstanceSummary(handle, summary) ||
	    !summary.valid || !summary.selected)
	    continue;
	const std::string candidate = semantic_path(summary.path.getString());
	if (candidate == target_path ||
	    (candidate.size() > target_path.size() &&
	     candidate.compare(0, target_path.size(), target_path) == 0 &&
	     candidate[target_path.size()] == '/'))
	    selected++;
    }
    return selected;
}


static int
exercise_draw_frontier(struct ged *gedp, BObolSceneController *scene)
{
    const char *nested_root = "nested_parent.c";
    const char *nested_target =
	"nested_parent.c/nested_child.c/nested_leaf.s";
    const char *nested_sibling =
	"nested_parent.c/nested_sibling.s";

    frontier_clear(gedp);
    if (frontier_draw_one(gedp, nested_root) != BRLCAD_OK)
	FAIL("frontier test should draw its nested root");
    if (!ged_selection_select_path(gedp, NULL, nested_target, 1))
	FAIL("frontier test should select the nested target");

    struct ged_scene_edit_request request;
    ged_scene_edit_request_init(&request);
    request.path = nested_target;
    request.view = ged_view_active_ctx(gedp);
    request.occurrences = GED_SCENE_EDIT_EXACT_OCCURRENCE;
    request.purpose = "frontier-regression";
    ged_scene_edit_scope_ref edit_scope = GED_SCENE_EDIT_SCOPE_REF_NULL;
    struct ged_scene_result *edit_result = ged_scene_result_create();
    if (ged_scene_edit_acquire(gedp, &request, &edit_scope,
	    edit_result) != GED_SCENE_OK ||
	    ged_scene_result_path_count(edit_result) != 1 ||
	    ged_scene_result_group_count(edit_result) != 1)
	FAIL("exact promotion should split one drawn root");
    ged_scene_result_destroy(edit_result);
    if (!scene || scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(scene, nested_root))
	FAIL("promotion must retain one LoD/data-owning database source");

    const std::set<std::string> exact_expected = {nested_root};
    const std::set<std::string> exact_actual = frontier_paths(gedp);
    if (exact_actual != exact_expected) {
	for (const std::string &path : exact_actual)
	    fprintf(stderr, "exact promotion path: %s\n", path.c_str());
	FAIL("exact promotion should retain the compact owning draw root");
    }
    if (ged_selection_count(gedp, NULL) != 1)
	FAIL("promotion must not change semantic selection membership");

    edit_result = ged_scene_result_create();
    if (ged_scene_edit_release(gedp, edit_scope,
	    GED_SCENE_EDIT_CANCEL, edit_result) != GED_SCENE_OK ||
	    ged_scene_result_conflict_count(edit_result) != 0 ||
	    frontier_paths(gedp) != std::set<std::string>{nested_root})
	FAIL("unchanged promotion should collapse to its original draw intent");
    if (ged_selection_count(gedp, NULL) != 1)
	FAIL("promotion collapse must preserve semantic selection membership");
    ged_scene_result_destroy(edit_result);

    /* The common librt edit buffer owns the same presentation lifecycle:
     * staging a primitive edit promotes every already-drawn occurrence, while
     * abandoning it collapses the retained owner without creating sources. */
    frontier_clear(gedp);
    if (frontier_draw_one(gedp, nested_root) != BRLCAD_OK)
	FAIL("edit-buffer promotion should draw its nested root");
    struct db_full_path edit_path;
    db_full_path_init(&edit_path);
    if (db_string_to_path(&edit_path, gedp->dbip, "nested_leaf.s") < 0)
	FAIL("edit-buffer promotion should resolve its primitive");
    struct bn_tol edit_tol = BN_TOL_INIT_TOL;
    struct rt_edit *edit_state =
	rt_edit_create(&edit_path, gedp->dbip, &edit_tol, NULL);
    if (!edit_state)
	FAIL("edit-buffer promotion should create an librt edit state");
    ged_edit_buf_set(gedp, &edit_path, edit_state);
    const std::set<std::string> edit_expected = {nested_root};
    if (ged_edit_buf_get(gedp, &edit_path) != edit_state ||
	frontier_paths(gedp) != edit_expected ||
	!scene || scene->getDatabaseSourceCount() != 1)
	FAIL("staged librt edit should promote all drawn occurrences on one retained source");

    ged_edit_session_ref edit_session = GED_EDIT_SESSION_REF_NULL;
    if (ged_edit_session_find(gedp, "nested_leaf.s", &edit_session) !=
	    GED_EDIT_OK)
	FAIL("staged librt edit should expose its session identity");
    char edit_feature_name[96] = {0};
    snprintf(edit_feature_name, sizeof(edit_feature_name),
	"_ged_edit_preview_%016llx_%016llx",
	(unsigned long long)edit_session.owner,
	(unsigned long long)edit_session.id);
    struct ged_view_context *edit_view = ged_view_active_ctx(gedp);
    if (!edit_view || !ged_view_feature_exists(edit_view,
	    edit_feature_name) || !ged_view_feature_visible(edit_view,
	    edit_feature_name))
	FAIL("a staged edit of a drawn occurrence should show its preview");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE,
	    nested_target, edit_view, -1, "erase active edit occurrence"))
	return 1;
    if (!ged_view_feature_visible(edit_view, edit_feature_name))
	FAIL("an edit preview should remain visible while another drawn occurrence remains");
    const char *renamed_occurrence =
	"nested_parent.c/nested_child_renamed.c/nested_leaf.s";
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE,
	    renamed_occurrence, edit_view, -1,
	    "erase remaining active edit occurrence"))
	return 1;
    if (ged_edit_buf_get(gedp, &edit_path) != edit_state ||
	!ged_view_feature_exists(edit_view, edit_feature_name) ||
	ged_view_feature_visible(edit_view, edit_feature_name))
	FAIL("erasing an edited occurrence should retain but hide its preview");
    if (frontier_draw_one(gedp, nested_target) != BRLCAD_OK ||
	!ged_view_feature_visible(edit_view, edit_feature_name))
	FAIL("redrawing an edited occurrence should restore its retained preview");
    ged_edit_buf_abandon(gedp, &edit_path);
    if (ged_edit_buf_get(gedp, &edit_path) ||
	frontier_paths(gedp) != std::set<std::string>{nested_root} ||
	scene->getDatabaseSourceCount() != 1)
	FAIL("abandoned librt edit should collapse its retained source frontier");
    db_free_full_path(&edit_path);

    /* A broad draw must adopt resident payloads and retire their narrow owner
     * as one publication.  Give the donor one probe occurrence that cannot be
     * produced by the database walk, so this callback contract does not depend
     * on whether asynchronous realization has already delivered every real
     * leaf to the broad source. */
    constexpr const char *resident_probe_path =
	"nested_parent.c/resident_adoption_probe.s";
    constexpr uint32_t resident_probe_occurrence =
	std::numeric_limits<uint32_t>::max();
    frontier_clear(gedp);
    if (frontier_draw_one(gedp, nested_root) != BRLCAD_OK)
	FAIL("resident-adoption test should draw its broad target");
    SoBRLDatabaseSource *broad_source = source_for_path(scene, nested_root);
    if (!broad_source || scene->getDatabaseSourceCount() != 1 ||
	broad_source->getCompactInstanceCount() <= 0)
	FAIL("resident-adoption test needs one broad compact source");
    BObolCompactOccurrence resident_probe;

    if (!broad_source->getCompactOccurrence(0, resident_probe) ||
	!resident_probe.geometry)
	FAIL("resident-adoption test needs reusable resident geometry");
    resident_probe.summary.path = resident_probe_path;
    resident_probe.summary.sourceName = "resident_adoption_probe.s";
    resident_probe.sourceMeshRequestValid = FALSE;
    resident_probe.viewDependentCsgGeometry = FALSE;
    resident_probe.lodBacked = FALSE;
    resident_probe.occurrenceIndex = resident_probe_occurrence;

    struct ged_bobol_publication_context resident_publication;
    if (!ged_bobol_publication_begin(&resident_publication, gedp,
	    ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_WIRE))
	FAIL("resident-adoption test should begin narrow-source publication");
    if (!ged_bobol_database_source_ensure_for_path(&resident_publication,
	    nested_target, gedp->dbip, GED_DRAW_MODE_WIRE,
	    broad_source->sourceRevision.getValue())) {
	ged_bobol_publication_end(&resident_publication);
	FAIL("resident-adoption test should stage its narrow donor");
    }
    ged_bobol_publication_end(&resident_publication);
    SoBRLDatabaseSource *narrow_source = source_for_path(scene, nested_target);
    if (!narrow_source || broad_source == narrow_source ||
	scene->getDatabaseSourceCount() != 2 ||
	frontier_paths(gedp) != std::set<std::string>{nested_root})
	FAIL("resident-adoption test should retain distinct broad and narrow owners");
    if (narrow_source->mergeCompactOccurrences({resident_probe}, TRUE) != 1 ||
	narrow_source->getCompactInstanceCountForPath(
	    resident_probe_path, TRUE) != 1)
	FAIL("resident-adoption test should add its unique donor occurrence");
    SbModernUtils::SoNodeRef broad_owner(broad_source);
    broad_source->clearSourceBounds();

    struct AdoptionObserver {
	BObolSceneController &scene;
	SoBRLDatabaseSource &narrow;
	const char *adoptedPath;
	SbString narrowKey;
	SoBRLDatabaseSource *broad;
	bool called = false;
	bool partial = false;
	bool failed = false;
	static void sourceChanged(void *data, SoSensor *)
	{
	    auto &self = *static_cast<AdoptionObserver *>(data);
	    if (self.called || !self.broad ||
		self.broad->getCompactInstanceCountForPath(
		    self.adoptedPath, TRUE) <= 0)
		return;
	    self.called = true;
	    self.partial = self.scene.getDatabaseSourceCount() != 1 ||
		self.scene.findDatabaseSourceInstance(
		    self.narrowKey.getString()) != nullptr;
	    if (self.partial)
		return;
	    BObolDatabaseSourceDisplayPatch patch;
	    patch.lineWidthValid = TRUE;
	    patch.lineWidth = 41;
	    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(
		self.broad->instanceKey.getValue().getString(), patch) <= 0;
	}
    } adoptionObserver{*scene, *narrow_source, resident_probe_path,
	narrow_source->instanceKey.getValue(), broad_source};
    SoNodeSensor adoptionSourceSensor(AdoptionObserver::sourceChanged,
	&adoptionObserver);
    adoptionSourceSensor.setPriority(0);
    adoptionSourceSensor.attach(broad_source);
    const char *resident_donors[] = {
	narrow_source->instanceKey.getValue().getString()
    };
    if (scene->subsumeDatabaseSourceInstances(
	    broad_source->instanceKey.getValue().getString(), resident_donors,
	    sizeof(resident_donors) / sizeof(resident_donors[0])) != 1)
	FAIL("resident-adoption test should subsume its narrow donor");
    adoptionSourceSensor.detach();
    broad_source = source_for_path(scene, nested_root);
    if (!adoptionObserver.called || adoptionObserver.partial ||
	adoptionObserver.failed ||
	broad_source != adoptionObserver.broad ||
	!broad_source || broad_source->lineWidth.getValue() != 41) {
	fprintf(stderr, "resident adoption called=%d partial=%d "
	    "failed=%d broad=%p observed=%p width=%d sources=%d adopted=%d "
	    "narrow-current=%p broad-key=%s observed-key=%s\n",
	    adoptionObserver.called,
	    adoptionObserver.partial,
	    adoptionObserver.failed, (void *)broad_source,
	    (void *)adoptionObserver.broad,
	    broad_source ? broad_source->lineWidth.getValue() : -1,
	    scene->getDatabaseSourceCount(), broad_source ?
	    broad_source->getCompactInstanceCountForPath(
		resident_probe_path, TRUE) : -1,
	    (void *)scene->findDatabaseSourceInstance(
		adoptionObserver.narrowKey.getString()),
	    broad_source ? broad_source->instanceKey.getValue().getString() : "<none>",
	    adoptionObserver.broad->instanceKey.getValue().getString());
	FAIL("resident adoption must retire narrow owners before broad-source observers");
    }
    std::puts("PASS GED resident-source adoption: complete handoff and current callback target");

    /* Repeat the real narrow-to-broad transition without the synthetic probe.
     * A later broad draw owns the intent, while selection and re-drawing a
     * covered narrow path remain semantic operations on that one source. */
    frontier_clear(gedp);
    if (frontier_draw_one(gedp, nested_target) != BRLCAD_OK)
	FAIL("subsumption test should draw its narrow path");
    narrow_source = source_for_path(scene, nested_target);
    if (!narrow_source || scene->getDatabaseSourceCount() != 1 ||
	narrow_source->getCompactInstanceCount() <= 0)
	FAIL("a direct narrow draw should create only a narrow lazy source");
    if (!ged_selection_select_path(gedp, NULL, nested_target, 1))
	FAIL("subsumption test should establish compact selection");
    (void)ged_selection_present_private(gedp);
    if (frontier_selected_occurrences_for_path(narrow_source,
	    nested_target) <= 0)
	FAIL("narrow source should expose the established compact selection");
    if (frontier_draw_one(gedp, nested_root) != BRLCAD_OK)
	FAIL("subsumption test should draw the broader root");
    broad_source = source_for_path(scene, nested_root);
    if (!broad_source || scene->getDatabaseSourceCount() != 1 ||
	broad_source->getCompactInstanceCountForPath(nested_target, TRUE) <= 0)
	FAIL("broad draw should adopt resident narrow data and retire its owner");
    if (frontier_selected_occurrences_for_path(broad_source,
	    nested_target) <= 0) {
	for (int i = 0; i < broad_source->getCompactInstanceCount(); i++) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    BObolCompactOccurrence direct;
	    const bool have_direct =
		broad_source->getCompactOccurrence(i, direct) ? true : false;
	    if (broad_source->getCompactInstanceHandle(i, handle) &&
		broad_source->getCompactInstanceSummary(handle, summary))
		fprintf(stderr, "broad occurrence: handle=%016llx:%016llx path=%s direct=%s selected=%d visible=%d\n",
		    (unsigned long long)handle.instanceWord0,
		    (unsigned long long)handle.instanceWord1,
		    summary.path.getString(), have_direct ?
		    direct.summary.path.getString() : "<none>",
		    summary.selected ? 1 : 0,
		    summary.visible ? 1 : 0);
	}
	FAIL("broad source adoption should preserve semantic selection visuals");
    }

    /* Compact roots do not manufacture one libged record per leaf.  A client
     * that genuinely needs an addressable occurrence (MGED sed/oed, for
     * example) must nevertheless resolve an exact or first-descendant path
     * without expanding siblings or duplicating the retained source. */
    struct ged_scene_path_request occurrence_request;
    ged_scene_path_request_init(&occurrence_request);
    occurrence_request.view = ged_view_active_ctx(gedp);
    occurrence_request.path = nested_target;
    occurrence_request.match = GED_SCENE_PATH_MATCH_EXACT;
    ged_scene_occurrence_ref exact_occurrence =
	ged_scene_occurrence_resolve(gedp, &occurrence_request);
    struct ged_scene_occurrence_info occurrence_info;
    if (ged_scene_occurrence_ref_is_null(exact_occurrence) ||
	!ged_scene_occurrence_get(gedp, exact_occurrence, &occurrence_info) ||
	!occurrence_info.fullpath || scene->getDatabaseSourceCount() != 1)
	FAIL("exact compact occurrence resolution should retain one owning source");
    char *resolved_occurrence_path = db_path_to_string(
	(struct db_full_path *)occurrence_info.fullpath);
    if (!resolved_occurrence_path ||
	!path_equal(resolved_occurrence_path, nested_target)) {
	if (resolved_occurrence_path)
	    bu_free(resolved_occurrence_path, "resolved compact occurrence path");
	FAIL("exact compact occurrence resolution should preserve semantic path identity");
    }
    bu_free(resolved_occurrence_path, "resolved compact occurrence path");

    occurrence_request.path = "nested_parent.c/nested_child.c";
    occurrence_request.match = GED_SCENE_PATH_MATCH_SUBTREE;
    ged_scene_occurrence_ref descendant_occurrence =
	ged_scene_occurrence_resolve(gedp, &occurrence_request);
    if (ged_scene_occurrence_ref_is_null(descendant_occurrence) ||
	!ged_scene_occurrence_get(gedp, descendant_occurrence,
	    &occurrence_info) || !occurrence_info.fullpath)
	FAIL("compact subtree occurrence resolution should find one visible descendant");
    resolved_occurrence_path = db_path_to_string(
	(struct db_full_path *)occurrence_info.fullpath);
    if (!resolved_occurrence_path ||
	!path_equal(resolved_occurrence_path, nested_target)) {
	if (resolved_occurrence_path)
	    bu_free(resolved_occurrence_path, "resolved compact subtree path");
	FAIL("compact subtree occurrence resolution should use deterministic scene order");
    }
    bu_free(resolved_occurrence_path, "resolved compact subtree path");

    occurrence_request.path = "nested_parent.c/missing.s";
    occurrence_request.match = GED_SCENE_PATH_MATCH_EXACT;
    if (!ged_scene_occurrence_ref_is_null(
	    ged_scene_occurrence_resolve(gedp, &occurrence_request)))
	FAIL("missing compact semantic paths must not create occurrence records");

    if (frontier_draw_one(gedp, nested_target) != BRLCAD_OK ||
	scene->getDatabaseSourceCount() != 1 ||
	!source_for_path(scene, nested_root))
	FAIL("drawing a path already covered by a broad intent must be idempotent");
    if (frontier_paths(gedp) != std::set<std::string>{nested_root})
	FAIL("a covered narrow draw must remain subsumed by the broad intent");

    struct ged_scene_erase_request erase;
    ged_scene_erase_request_init(&erase);
    erase.path = nested_target;
    erase.view = ged_view_active_ctx(gedp);
    struct ged_scene_result *erase_result = ged_scene_result_create();
    if (ged_scene_erase(gedp, &erase, erase_result) != GED_SCENE_OK ||
	    frontier_paths(gedp) != std::set<std::string>{nested_root} ||
	    ged_scene_path_state_get(gedp, ged_draw_active_view_ctx(gedp),
		nested_root, GED_SCENE_DRAW_DEFAULT) !=
		GED_SCENE_PATH_PARTIALLY_DRAWN ||
	    ged_scene_path_state_get(gedp, ged_draw_active_view_ctx(gedp),
		nested_target, GED_SCENE_DRAW_DEFAULT) !=
		GED_SCENE_PATH_NOT_DRAWN ||
	    ged_scene_path_state_get(gedp, ged_draw_active_view_ctx(gedp),
		nested_sibling, GED_SCENE_DRAW_DEFAULT) != GED_SCENE_PATH_DRAWN)
	FAIL("nested erase should retain one compact owner with partial/hidden/full path states");
    ged_scene_result_destroy(erase_result);
    SoBRLDatabaseSource *nested_source =
	source_for_path(scene, nested_root);
    if (!nested_source || nested_source->getCompactInstanceCount() != 3 ||
	    frontier_visible_occurrences(nested_source) != 2 ||
	    scene->getDatabaseSourceCount() != 1)
	FAIL("nested erase should mask one occurrence without duplicating or rebuilding its source");
    const int redraw_target_ret = frontier_draw_one(gedp, nested_target);
    const std::set<std::string> collapsed_paths = frontier_paths(gedp);
    if (redraw_target_ret != BRLCAD_OK ||
	collapsed_paths != std::set<std::string>{nested_root} ||
	scene->getDatabaseSourceCount() != 1 ||
	frontier_visible_occurrences(nested_source) != 3) {
	fprintf(stderr, "frontier collapse: ret=%d sources=%d visible=%d paths=",
	    redraw_target_ret, scene->getDatabaseSourceCount(),
	    frontier_visible_occurrences(nested_source));
	for (const std::string &collapsed_path : collapsed_paths)
	    fprintf(stderr, " %s", collapsed_path.c_str());
	fprintf(stderr, " frontier_active=%d\n",
	    nested_source->hasCompactInstanceVisibilityFrontier() ? 1 : 0);
	for (int i = 0; i < nested_source->getCompactInstanceCount(); i++) {
	    BObolCompactOccurrence occurrence;
	    if (nested_source->getCompactOccurrence(i, occurrence))
		fprintf(stderr, "  direct %s visible=%d\n",
		    occurrence.summary.path.getString(),
		    occurrence.summary.visible ? 1 : 0);
	}
	FAIL("re-drawing an erased subtree should collapse the retained source frontier");
    }
    if (ged_selection_count(gedp, NULL) != 1 ||
	frontier_selected_occurrences_for_path(nested_source,
	    nested_target) <= 0)
	FAIL("frontier collapse should preserve semantic and visual selection");

    /* Shaded-BOT roots use the same compact occurrence owner but retain a
     * distinct representation-mode identity.  Exercise the public -m1 path
     * explicitly: a mode mismatch in the frontier bridge otherwise reports a
     * successful semantic erase while leaving the shaded occurrence visible. */
    frontier_clear(gedp);
    if (frontier_draw_one_mode(gedp, 1, nested_root) != BRLCAD_OK)
	FAIL("shaded frontier test should draw its nested root");
    nested_source = source_for_path(scene, nested_root);
    BObolDatabaseSourceSummary shaded_summary;
    if (!nested_source || !nested_source->getSummary(shaded_summary) ||
	!shaded_summary.valid ||
	shaded_summary.representationMode !=
	    SoBRLDatabaseSource::REPRESENTATION_SHADED_BOTS ||
	frontier_visible_occurrences(nested_source) != 3)
	FAIL("shaded frontier root should retain three compact occurrences");
    ged_scene_erase_request_init(&erase);
    erase.path = nested_target;
    erase.view = ged_view_active_ctx(gedp);
    erase_result = ged_scene_result_create();
    if (ged_scene_erase(gedp, &erase, erase_result) != GED_SCENE_OK ||
	frontier_visible_occurrences(nested_source) != 2)
	FAIL("shaded nested erase should mask one retained occurrence");
    ged_scene_result_destroy(erase_result);
    if (frontier_draw_one_mode(gedp, 1, nested_target) != BRLCAD_OK ||
	frontier_visible_occurrences(nested_source) != 3)
	FAIL("shaded nested redraw should restore the retained occurrence");

    frontier_clear(gedp);
    if (frontier_draw_one(gedp, "reuse_root.c") != BRLCAD_OK)
	FAIL("all-occurrence promotion should draw the reuse root");
    request.path =
	"reuse_root.c/reuse_inst_a.c/reuse_shared.c/reuse_leaf.s";
    request.occurrences = GED_SCENE_EDIT_ALL_DRAWN_OCCURRENCES;
    edit_scope = GED_SCENE_EDIT_SCOPE_REF_NULL;
    edit_result = ged_scene_result_create();
    if (ged_scene_edit_acquire(gedp, &request, &edit_scope,
	    edit_result) != GED_SCENE_OK ||
	    ged_scene_result_path_count(edit_result) != 2 ||
	    ged_scene_result_group_count(edit_result) != 1 ||
	    frontier_paths(gedp) != std::set<std::string>{"reuse_root.c"})
	FAIL("all-occurrence promotion should retain one compact owner per root");
    ged_scene_result_destroy(edit_result);
    edit_result = ged_scene_result_create();
    if (ged_scene_edit_release(gedp, edit_scope,
	    GED_SCENE_EDIT_COMMIT, edit_result) != GED_SCENE_OK ||
	    frontier_paths(gedp) != std::set<std::string>{"reuse_root.c"})
	FAIL("all-occurrence promotion should collapse without duplicate roots");
    ged_scene_result_destroy(edit_result);

    frontier_clear(gedp);
    if (frontier_draw_one(gedp, "dup_twice.c") != BRLCAD_OK)
	FAIL("duplicate-occurrence promotion should draw its root");
    request.path = "dup_twice.c/dup_leaf.s@1";
    request.occurrences = GED_SCENE_EDIT_EXACT_OCCURRENCE;
    edit_scope = GED_SCENE_EDIT_SCOPE_REF_NULL;
    edit_result = ged_scene_result_create();
    if (ged_scene_edit_acquire(gedp, &request, &edit_scope,
	    edit_result) != GED_SCENE_OK ||
	    frontier_paths(gedp) != std::set<std::string>{"dup_twice.c"})
	FAIL("exact promotion should preserve duplicate child occurrences");
    ged_scene_result_destroy(edit_result);

    ged_scene_erase_request_init(&erase);
    erase.path = request.path;
    erase.view = ged_view_active_ctx(gedp);
    erase_result = ged_scene_result_create();
    if (ged_scene_erase(gedp, &erase, erase_result) != GED_SCENE_OK ||
	!ged_scene_result_changed(erase_result))
	FAIL("draw change inside a promotion should succeed");
    ged_scene_result_destroy(erase_result);
    edit_result = ged_scene_result_create();
    const enum ged_scene_status duplicate_release = ged_scene_edit_release(
	gedp, edit_scope, GED_SCENE_EDIT_CANCEL, edit_result);
    const std::set<std::string> duplicate_frontier = frontier_paths(gedp);
    const enum ged_scene_path_state duplicate_state =
	ged_scene_path_state_get(gedp, ged_draw_active_view_ctx(gedp),
	    request.path, GED_SCENE_DRAW_DEFAULT);
    if (duplicate_release != GED_SCENE_OK ||
	    ged_scene_result_conflict_count(edit_result) != 1 ||
	    duplicate_frontier != std::set<std::string>{"dup_twice.c"} ||
	    duplicate_state != GED_SCENE_PATH_NOT_DRAWN) {
	fprintf(stderr,
	    "duplicate conflict: release=%d conflicts=%zu state=%d paths=",
	    static_cast<int>(duplicate_release),
	    ged_scene_result_conflict_count(edit_result), duplicate_state);
	for (const std::string &path : duplicate_frontier)
	    fprintf(stderr, " %s", path.c_str());
	fprintf(stderr, "\n");
	FAIL("same-root conflict should preserve the user's current frontier");
    }
    ged_scene_result_destroy(edit_result);

    (void)ged_selection_clear(gedp, NULL);
    frontier_clear(gedp);
    return 0;
}


static int
exercise_metadata_effects(struct ged *gedp, BObolSceneController *scene)
{
    struct GedMetadataObserver {
	struct Failure {
	};
	BObolSceneController &scene;
	uint64_t frame;
	bool throwing, called = false, partial = false;
	static void
	changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<GedMetadataObserver *>(data);
	    self.called = true;
	    self.partial = self.scene.getFrameRevision() != self.frame + 1;
	    if (self.throwing)
		throw Failure();
	}
    };

    bool passed = true;
    for (bool compact : {false, true})
	for (bool throwing : {false, true}) {
	    if (compact && throwing)
		continue;
	    scene->replaceDatabaseSourceInstance(
	        "metadata-owner", "box.s", gedp->dbip, SoBRLDatabaseSource::WIREFRAME, 1);
	    auto *source = scene->findDatabaseSourceInstance("metadata-owner");
	    if (!source)
		FAIL("GED metadata fixture preparation failed");
	    if (compact) {
		if (!source->realizeDatabaseWireframe() || !source->hasCompactInstanceIndex())
		    FAIL("GED metadata fixture preparation failed");
	    } else {
		auto *wire = new SoBRLVListShape;
		wire->sourcePath = "box.s";
		source->addChild(wire);
	    }
	    GedMetadataObserver observed{*scene, scene->getFrameRevision(), throwing};
	    SoFieldSensor sensor(GedMetadataObserver::changed, &observed);
	    sensor.setPriority(0);
	    sensor.attach(&source->databaseRegionId);
	    BObolDrawMetadataRecord record;
	    bobol_draw_metadata_record_init(&record);
	    record.directoryFound = 1;
	    record.hasRegionId = 1;
	    record.regionId = 7;
	    int changed = 0;
	    bool threw = false;
	    try {
		changed = ged_draw_obol_database_source_apply_draw_metadata_for_path(
		    gedp, "box.s", &record);
	    } catch (const GedMetadataObserver::Failure &) {
		threw = true;
	    }
	    sensor.detach();
	    const bool frame = scene->getFrameRevision() == observed.frame + 1;
	    std::printf(
	        "ged-metadata compact=%d applied=%d observer_partial=%d frame_advanced=%d\n",
	        compact, changed > 0, observed.partial, frame);
	    passed = passed && (throwing || changed > 0) && threw == throwing &&
	             observed.called == !compact && !observed.partial && frame &&
	             ged_draw_obol_database_source_apply_draw_metadata_for_path(
	                 gedp, "box.s", &record) == 1 &&
	             scene->getFrameRevision() == observed.frame + 1;
	    scene->removeDatabaseSourceInstance("metadata-owner");
	}
    for (bool replace : {false, true}) {
	constexpr std::array<const char *, 2> keys{{"refresh-first", "refresh-second"}};
	std::vector<SbModernUtils::SoNodeRef> owners;
	for (size_t i = 0; i < keys.size(); ++i) {
	    scene->replaceDatabaseSourceInstance(
	        keys[i], "box.s", gedp->dbip, SoBRLDatabaseSource::WIREFRAME, 1);
	    owners.emplace_back(scene->findDatabaseSourceInstance(keys[i]));
	}
	struct State {
	    BObolSceneController &scene;
	    db_i *database;
	    const std::array<const char *, 2> &keys;
	    bool replace, called = false, valid = true;
	    size_t later = 0;
	} state{*scene, gedp->dbip, keys, replace};
	struct Observer {
	    State &state;
	    size_t index;
	    static void
	    changed(void *data, SoSensor *)
	    {
		auto &self = *static_cast<Observer *>(data);
		auto &s = self.state;
		if (s.called)
		    return;
		s.called = true;
		s.later = 1 - self.index;
		if (s.replace) {
		    s.valid = s.scene.removeDatabaseSourceInstance(s.keys[s.later]) > 0 &&
		              s.scene.replaceDatabaseSourceInstance(s.keys[s.later], "box.s",
		                  s.database, SoBRLDatabaseSource::WIREFRAME, 1) > 0;
		} else
		    s.valid = s.scene.refreshDatabaseSourceInstanceMaterialColorFromDatabase(
		                  s.keys[s.later], 2, s.database) > 0;
	    }
	} observers[2] = {{state, 0}, {state, 1}};
	SoFieldSensor a(Observer::changed, &observers[0]), b(Observer::changed, &observers[1]);
	a.setPriority(0);
	b.setPriority(0);
	a.attach(&static_cast<SoBRLDatabaseSource *>(owners[0].get())->materialRevision);
	b.attach(&static_cast<SoBRLDatabaseSource *>(owners[1].get())->materialRevision);
	const int changed = ged_draw_obol_database_source_refresh_material_color_for_path(
	    gedp, "box.s", gedp->dbip, 1);
	a.detach();
	b.detach();
	auto *later = scene->findDatabaseSourceInstance(keys[state.later]);
	const bool retained = state.called && state.valid && later &&
	                      later->materialRevision.getValue() == (replace ? 0u : 2u);
	std::printf("ged-refresh replace=%d changed=%d later_retained=%d\n", replace, changed > 0,
	    retained);
	passed = passed && changed > 0 && retained;
	for (auto *key : keys)
	    scene->removeDatabaseSourceInstance(key);
    }
    if (!passed)
	FAIL("GED metadata effects must commit frames and retain later refresh targets");
    std::puts(
        "PASS GED metadata effects: source/compact frames and retained later edits/replacements");
    return 0;
}

static int
exercise_evaluated_region_publication(struct ged *gedp, BObolSceneController *scene)
{
    constexpr const char *key = "evaluated-region-probe-owner";
    constexpr const char *path = "evaluated-region-probe.c";
    for (bool throwing : {false, true}) {
	if (scene->replaceDatabaseSourceInstance(key, path, gedp->dbip, SoBRLDatabaseSource::WIREFRAME, 1) < 0)
	    FAIL("evaluated region probe source creation failed");
	auto *source = scene->findDatabaseSourceInstance(key);
	if (!source) FAIL("evaluated region probe source missing");
	auto *wire0 = new SoBRLVListShape;
	auto *wire1 = new SoBRLVListShape;
	auto *mesh0 = new SoBRLMeshShape;
	auto *mesh1 = new SoBRLMeshShape;
	source->addChild(wire0); source->addChild(wire1); source->addChild(mesh0); source->addChild(mesh1);
	wire0->sourcePath = wire1->sourcePath = mesh0->sourcePath = mesh1->sourcePath = path;
	struct ObserverFailure {};
	struct Observer {
	    static void changed(void *data, SoSensor *) {
		auto &self = *static_cast<Observer *>(data);
		++self.calls;
		for (auto *field : self.fields) self.partial = self.partial || field->getValue() != 1;
		self.partial = self.partial || self.scene->getFrameRevision() != self.frame + 1;
		if (self.throwing) throw ObserverFailure();
	    }
	    std::array<SoSFInt32 *, 4> fields;
	    BObolSceneController *scene;
	    uint64_t frame;
	    bool throwing, partial = false;
	    size_t calls = 0;
	} observer{{&wire0->regionId, &wire1->regionId, &mesh0->regionId, &mesh1->regionId},
	    scene, scene->getFrameRevision(), throwing};
	SoFieldSensor sensor(Observer::changed, &observer); sensor.setPriority(0); sensor.attach(&wire0->regionId);
	bool threw = false;
	int changed = 0;
	try { changed = ged_draw_obol_database_source_set_evaluated_region_for_path(gedp, path, 1); }
	catch (const ObserverFailure &) { threw = true; }
	sensor.detach();
	bool complete = observer.calls && !observer.partial && threw == throwing && (throwing || changed == 1) &&
	    scene->getFrameRevision() == observer.frame + 1;
	for (auto *field : observer.fields) complete = complete && field->getValue() == 1;
	const auto frame = scene->getFrameRevision();
	complete = complete && ged_draw_obol_database_source_set_evaluated_region_for_path(gedp, path, 1) == 0 &&
	    scene->getFrameRevision() == frame;
	scene->removeDatabaseSourceInstance(key);
	if (!complete) FAIL("GED evaluated region must publish every wire/mesh and the scene frame before observers");
    }
    std::puts("PASS GED evaluated-region publication: complete multi-shape state, throwing observers and no-op frame retention");
    return 0;
}

static int
exercise_redraw_appearance(struct ged *gedp, BObolSceneController *scene)
{
    enum class Edit { Observe, Change, Replace, Remove, GroupChange };
    constexpr int initialWidth = 3, editedWidth = 7;
    const auto matches = [](const SoBRLDatabaseSource &source, const BObolDatabaseSourceSummary &expected) {
	return source.visible.getValue() == expected.visible && source.selected.getValue() == expected.selected &&
	    source.highlighted.getValue() == expected.highlighted && source.lineStyle.getValue() == expected.lineStyle &&
	    source.lineWidth.getValue() == expected.lineWidth &&
	    std::abs(source.transparency.getValue() - expected.transparency) < std::numeric_limits<float>::epsilon() &&
	    source.colorOverride.getValue() == expected.colorOverride && source.color.getValue() == expected.color &&
	    source.materialColorValid.getValue() == expected.materialColorValid &&
	    source.materialColor.getValue() == expected.materialColor && source.materialRevision.getValue() == expected.materialRevision;
    };
    for (int mode : {GED_DRAW_MODE_WIRE, GED_DRAW_MODE_SHADED_BOTS, GED_DRAW_MODE_SHADED,
	GED_DRAW_MODE_EVAL_WIRE, GED_DRAW_MODE_HIDDEN_LINE, GED_DRAW_MODE_EVAL_POINTS})
    for (auto edit : {Edit::Observe, Edit::Change, Edit::Replace, Edit::Remove, Edit::GroupChange}) {
	char option[16]; std::snprintf(option, sizeof(option), "-m%d", mode);
	const char *draw[] = {"draw", option, "box.s"};
	if (ged_exec_draw(gedp, 3, draw) != BRLCAD_OK) FAIL("redraw appearance fixture draw");
	auto *source = source_for_representation(scene, "box.s", mode);
	if (!source) FAIL("redraw appearance fixture source");
	SbModernUtils::SoNodeRef owner(source);
	const SbString key = source->instanceKey.getValue();
	if (scene->setDatabaseSourceInstanceState(key.getString(), TRUE, source->sourceRevision.getValue(),
	    source->inputsRevision.getValue(), TRUE, TRUE, TRUE, 1, initialWidth, 0.375f,
	    TRUE, SbColor(0.25f, 0.5f, 0.75f), TRUE, SbColor(0.125f, 0.25f, 0.5f), 17) < 0 ||
	    source->clearRealizedGeometry(TRUE) < 0) FAIL("redraw appearance fixture state");
	source->markStale(SoBRLDatabaseSource::STALE_SOURCE);
	// A real configuration change guarantees a callback even when appearance
	// itself remains unchanged throughout the redraw.
	source->setMaterialPolicyState(SoBRLDatabaseSource::MATERIAL_INHERIT);
	BObolDatabaseSourceSummary before;
	if (!source->getSummary(before) || !before.valid) FAIL("redraw appearance fixture summary");
	struct Observer {
	    BObolSceneController &scene;
	    SoBRLDatabaseSource &source;
	    const BObolDatabaseSourceSummary &before;
	    const decltype(matches) &compare;
	    struct db_i *database;
	    Edit edit;
	    bool called = false, partial = false, failed = false;
	    static void changed(void *data, SoSensor *)
	    {
		auto &self = *static_cast<Observer *>(data);
		if (self.called) return;
		self.called = true;
		self.partial = !self.compare(self.source, self.before);
		if (self.edit == Edit::Observe) return;
		const char *instance = self.before.instanceKey.getString();
		if (self.edit == Edit::GroupChange) {
		    auto *group = static_cast<SoBRLSceneGroup *>(self.scene.findGroup("box.s"));
		    self.failed = !group;
		    if (!group) return;
		    self.partial = self.partial || group->lineWidth.getValue() != self.before.lineWidth;
		    self.failed = self.scene.setGroupDisplayState("box.s", group->visible.getValue(), group->selected.getValue(),
			group->highlighted.getValue(), group->lineStyle.getValue(), editedWidth, group->transparency.getValue(),
			group->colorOverride.getValue(), group->color.getValue(), group->materialColorValid.getValue(),
			group->materialColor.getValue(), group->materialRevision.getValue()) < 0;
		    return;
		}
		if (self.edit == Edit::Change) {
		    BObolDatabaseSourceDisplayPatch patch;
		    patch.lineWidthValid = TRUE; patch.lineWidth = editedWidth;
		    patch.transparencyValid = TRUE; patch.transparency = 0.625f;
		    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(instance, patch) < 0;
		    return;
		}
		self.failed = self.scene.removeDatabaseSourceInstance(instance) <= 0;
		if (self.edit == Edit::Remove) return;
		BObolDatabaseSourcePublishState replacement;
		replacement.sourceInstanceKey = instance; replacement.sourcePath = self.before.path.getString();
		replacement.sourceRepresentationKey = self.before.representationKey.getString();
		replacement.database = self.database; replacement.drawMode = self.before.drawMode;
		replacement.representationMode = self.before.representationMode;
		replacement.sourceRevisionValid = TRUE; replacement.sourceRevision = self.before.sourceRevision + 1;
		replacement.lineWidth = editedWidth;
		replacement.roleFlagsValid = TRUE; replacement.roleFlags = self.before.realizationRoleFlags;
		self.failed = self.scene.publishDatabaseSourceInstance(replacement) < 0 || self.failed;
	    }
	} observer{*scene, *source, before, matches, gedp->dbip, edit};
	SoNodeSensor sensor(Observer::changed, &observer); sensor.setPriority(0); sensor.attach(source);
	auto request = ged_scene_reducer_request_make(GED_SCENE_REDUCER_SOURCE_UPDATED, "box.s");
	request.view = ged_draw_active_view_ctx(gedp); request.mode = mode; request.redraw = 1;
	request.stale_reason = GED_DRAW_STALE_SETTINGS_CHANGED;
	const int changed = ged_draw_obol_scene_sync_transaction(gedp, &request, nullptr, scene);
	sensor.detach();
	auto *current = scene->findDatabaseSourceInstance(key.getString());
	if (changed <= 0 || !observer.called || observer.partial || observer.failed)
	    FAIL("redraw exposed temporary source appearance or failed a callback edit");
	if (edit == Edit::Remove) {
	    if (current) FAIL("redraw resurrected a removed callback target");
	} else {
	    if (!current || current->representationMode.getValue() != mode ||
		(!current->hasRealizedWireGeometry() && !current->hasRealizedMeshGeometry()))
		FAIL("redraw lost current representation or geometry");
	    auto expected = before;
	    if (edit == Edit::Change) { expected.lineWidth = editedWidth; expected.transparency = 0.625f; }
	    if (edit == Edit::Replace) {
		if (current == source || current->lineWidth.getValue() != editedWidth ||
		    current->sourceRevision.getValue() != before.sourceRevision + 1)
		    FAIL("redraw modified a replacement callback target");
	    } else if (!matches(*current, expected)) FAIL("redraw overwrote current source appearance");
	}
	if (edit == Edit::GroupChange) {
	    auto *group = static_cast<SoBRLSceneGroup *>(scene->findGroup("box.s"));
	    if (!group || group->lineWidth.getValue() != editedWidth) FAIL("redraw overwrote a newer callback group edit");
	}
	const char *erase[] = {"erase", "box.s"};
	if (ged_exec_erase(gedp, 2, erase) != BRLCAD_OK) FAIL("redraw appearance fixture cleanup");
	std::printf("PASS GED redraw appearance mode%d edit%d: complete appearance and current callback target\n", mode, int(edit));
    }
    return 0;
}

static int
exercise_redraw_erase_publication(struct ged *gedp,
	BObolSceneController *scene)
{
    if (!gedp || !scene)
	FAIL("redraw/erase publication test needs GED and Obol scene state");

    constexpr const char *root_path = "nested_parent.c";
    frontier_clear(gedp);
    const char *draw[] = {"draw", root_path};
    if (ged_exec_draw(gedp, 2, draw) != BRLCAD_OK)
	FAIL("redraw/erase publication fixture draw");
    (void)scene->realizePending();

    SoBRLDatabaseSource *source = source_for_representation(scene,
	root_path, SoBRLDatabaseSource::REPRESENTATION_WIRE);
    BObolDatabaseSourceSummary before;
    if (!source || !source->getSummary(before) || !before.valid ||
	source->getCompactInstanceCount() <= 0)
	FAIL("redraw/erase publication fixture needs a compact root source");
    SbModernUtils::SoNodeRef source_owner(source);

    struct Observer {
	BObolSceneController &scene;
	SbString sourceKey;
	const char *groupPath;
	bool called = false;
	bool partial = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<Observer *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    self.partial =
		self.scene.findDatabaseSourceInstance(
		    self.sourceKey.getString()) != NULL ||
		self.scene.findGroup(self.groupPath) != NULL;
	}
    } observer{*scene, before.instanceKey, root_path};
    SoNodeSensor sensor(Observer::changed, &observer);
    sensor.setPriority(0);
    sensor.attach(scene->getSceneRoot());

    struct ged_scene_reducer_request request =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_ERASE, root_path);
    request.view = ged_draw_active_view_ctx(gedp);
    request.mode = GED_DRAW_MODE_WIRE;
    const int changed = ged_draw_obol_scene_sync_transaction(gedp,
	&request, NULL, scene);
    sensor.detach();

    if (changed <= 0 || !observer.called)
	FAIL("compact root erase should publish an observable scene change");
    if (observer.partial)
	FAIL("compact root erase exposed authored visibility before source and group retirement");
    if (scene->findDatabaseSourceInstance(before.instanceKey.getString()) ||
	scene->findGroup(root_path))
	FAIL("compact root erase should retire its source and owning group");

    frontier_clear(gedp);
    const char *draw_first[] = {"draw", "nested_parent.c"};
    const char *draw_second[] = {"draw", "reuse_root.c"};
    if (ged_exec_draw(gedp, 2, draw_first) != BRLCAD_OK ||
	ged_exec_draw(gedp, 2, draw_second) != BRLCAD_OK)
	FAIL("redraw callback fixture draws");
    (void)scene->realizePending();
    if (scene->getDatabaseSourceCount() != 2)
	FAIL("redraw callback fixture needs two compact sources");
    SoBRLDatabaseSource *first = scene->getDatabaseSource(0);
    SoBRLDatabaseSource *second = scene->getDatabaseSource(1);
    BObolDatabaseSourceSummary first_summary;
    BObolDatabaseSourceSummary second_summary;
    if (!first || !second || !first->getSummary(first_summary) ||
	!second->getSummary(second_summary) ||
	first->getCompactInstanceCount() <= 0 ||
	second->getCompactInstanceCount() <= 0)
	FAIL("redraw callback fixture compact sources");
    if (first->clearCompactInstanceVisibilityFrontier() < 0 ||
	second->clearCompactInstanceVisibilityFrontier() < 0 ||
	first->setCompactInstanceDisplayStateForPath("", TRUE,
	    1, FALSE, 0, FALSE, 0, FALSE) <= 0 ||
	second->setCompactInstanceDisplayStateForPath("", TRUE,
	    1, FALSE, 0, FALSE, 0, FALSE) <= 0)
	FAIL("redraw callback fixture hidden state");
    SbModernUtils::SoNodeRef first_owner(first), second_owner(second);
    struct RedrawObserver {
	BObolSceneController &scene;
	SbString sourceKey;
	bool called = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<RedrawObserver *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    self.failed = self.scene.removeDatabaseSourceInstance(
		self.sourceKey.getString()) <= 0;
	}
    } redraw_observer{*scene, first_summary.instanceKey};
    SoNodeSensor redraw_sensor(RedrawObserver::changed, &redraw_observer);
    redraw_sensor.setPriority(0);
    redraw_sensor.attach(first);
    struct ged_scene_reducer_request redraw =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_REDRAW, NULL);
    redraw.mode = -1;
    const int redrawn = ged_draw_obol_scene_sync_transaction(gedp,
	&redraw, NULL, scene);
    redraw_sensor.detach();
    SoBRLDatabaseSource *current_second =
	scene->findDatabaseSourceInstance(second_summary.instanceKey.getString());
    if (redrawn <= 0 || !redraw_observer.called || redraw_observer.failed ||
	scene->findDatabaseSourceInstance(first_summary.instanceKey.getString()) ||
	current_second != second ||
	frontier_visible_occurrences(current_second) !=
	    current_second->getCompactInstanceCount())
	FAIL("redraw skipped a compact source after callback removal changed scene order");

    if (current_second->clearCompactInstanceVisibilityFrontier() < 0 ||
	current_second->setCompactInstanceVisibilityOverrideForPathMatch("",
	    BOBOL_COMPACT_PATH_SUBTREE, FALSE) <= 0 ||
	scene->setDatabaseSourceInstanceRealizationRoleFlags(
	    second_summary.instanceKey.getString(),
	    SoBRLDatabaseSource::REALIZATION_ROLE_NONE) < 0)
	FAIL("redraw realization fixture hidden state");
    current_second->markStale(SoBRLDatabaseSource::STALE_SOURCE);
    struct ged_scene_reducer_request realization_redraw =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_REDRAW,
	    second_summary.path.getString());
    realization_redraw.view = ged_draw_active_view_ctx(gedp);
    realization_redraw.mode = -1;
    const int accepted_redraw = ged_draw_obol_scene_sync_transaction(gedp,
	&realization_redraw, NULL, scene);
    if (accepted_redraw <= 0 || current_second->needsRealization() ||
	frontier_visible_occurrences(current_second) !=
	    current_second->getCompactInstanceCount()) {
	std::fprintf(stderr, "accepted redraw ret=%d needs=%d visible=%d total=%d visited=%u realized=%u failed=%u diagnostic=%s key=%s path=%s roles=%d\n",
	    accepted_redraw, current_second->needsRealization(),
	    frontier_visible_occurrences(current_second),
	    current_second->getCompactInstanceCount(),
	    scene->getLastVisitedSourceCount(), scene->getLastRealizedSourceCount(),
	    scene->getLastFailedSourceCount(), scene->getLastDiagnostics().getString(),
	    second_summary.instanceKey.getString(), second_summary.path.getString(),
	    current_second->realizationRoleFlags.getValue());
	FAIL("view-scoped redraw did not realize its accepted compact source");
    }

    if (current_second->clearCompactInstanceVisibilityFrontier() < 0 ||
	current_second->setCompactInstanceVisibilityOverrideForPathMatch("",
	    BOBOL_COMPACT_PATH_SUBTREE, FALSE) <= 0)
	FAIL("redraw realization acceptance hidden state");
    current_second->markStale(SoBRLDatabaseSource::STALE_SOURCE);
    BObolDatabaseSourceSummary realization_before;
    if (!current_second->getSummary(realization_before) ||
	!realization_before.valid || !current_second->needsRealization())
	FAIL("redraw realization acceptance needs a pending source");

    constexpr int replacement_line_width = 37;
    struct ReplacementObserver {
	int calls = 0;
	static void changed(void *data, SoSensor *)
	{
	    ++static_cast<ReplacementObserver *>(data)->calls;
	}
    } replacement_observer;
    SoNodeSensor replacement_sensor(ReplacementObserver::changed,
	&replacement_observer);
    replacement_sensor.setPriority(0);
    struct RealizationObserver {
	BObolSceneController &scene;
	SoBRLDatabaseSource &source;
	BObolDatabaseSourceSummary before;
	struct db_i *database;
	SoNodeSensor &replacementSensor;
	bool called = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<RealizationObserver *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    const char *key = self.before.instanceKey.getString();
	    self.failed = self.scene.removeDatabaseSourceInstance(key) <= 0;
	    if (self.failed)
		return;
	    BObolDatabaseSourcePublishState replacement;
	    replacement.sourceInstanceKey = key;
	    replacement.sourcePath = self.before.path.getString();
	    replacement.sourceRepresentationKey =
		self.before.representationKey.getString();
	    replacement.database = self.database;
	    replacement.drawMode = self.before.drawMode;
	    replacement.representationMode = self.before.representationMode;
	    replacement.sourceRevisionValid = TRUE;
	    replacement.sourceRevision = self.before.sourceRevision + 1;
	    replacement.lineWidth = replacement_line_width;
	    replacement.roleFlagsValid = TRUE;
	    replacement.roleFlags = SoBRLDatabaseSource::REALIZATION_ROLE_NONE;
	    self.failed = self.scene.publishDatabaseSourceInstance(replacement) <= 0;
	    SoBRLDatabaseSource *published =
		self.scene.findDatabaseSourceInstance(key);
	    if (!self.failed && published)
		self.replacementSensor.attach(published);
	}
    } realization_observer{*scene, *current_second, realization_before,
	gedp->dbip, replacement_sensor};
    SoNodeSensor realization_sensor(RealizationObserver::changed,
	&realization_observer);
    realization_sensor.setPriority(0);
    realization_sensor.attach(current_second);
    const int realization_redrawn = ged_draw_obol_scene_sync_transaction(gedp,
	&realization_redraw, NULL, scene);
    realization_sensor.detach();
    replacement_sensor.detach();
    SoBRLDatabaseSource *replacement = scene->findDatabaseSourceInstance(
	realization_before.instanceKey.getString());
    if (realization_redrawn <= 0 || !realization_observer.called ||
	realization_observer.failed || replacement_observer.calls != 0 || !replacement ||
	replacement == current_second ||
	replacement->lineWidth.getValue() != replacement_line_width ||
	!replacement->needsRealization() || replacement->hasRealizedWireGeometry() ||
	replacement->hasRealizedMeshGeometry())
	FAIL("redraw realized a same-key source published by its visibility callback");

    (void)scene->clearDatabaseSources();
    frontier_clear(gedp);
    const char *draw_wire[] = {"draw", "-m0", root_path};
    const char *add_shaded[] = {"draw", "-m2", "--add-mode", root_path};
    if (ged_exec_draw(gedp, 3, draw_wire) != BRLCAD_OK ||
	ged_exec_draw(gedp, 4, add_shaded) != BRLCAD_OK)
	FAIL("mode-scoped compact erase fixture draws");
    (void)scene->realizePending();
    SoBRLDatabaseSource *wire = source_for_representation(scene, root_path,
	SoBRLDatabaseSource::REPRESENTATION_WIRE);
    SoBRLDatabaseSource *shaded = source_for_representation(scene, root_path,
	SoBRLDatabaseSource::REPRESENTATION_SHADED);
    if (!wire || !shaded || wire == shaded ||
	shaded->clearCompactInstanceVisibilityFrontier() < 0 ||
	frontier_visible_occurrences(shaded) != shaded->getCompactInstanceCount())
	FAIL("mode-scoped compact erase fixture sources");
    if (shaded->setCompactInstanceVisibilityOverrideForPathMatch(root_path,
	    BOBOL_COMPACT_PATH_SUBTREE, FALSE) <= 0)
	FAIL("mode-scoped redraw fixture hidden sibling representation");
    struct ged_scene_reducer_request redraw_wire =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_REDRAW, root_path);
    redraw_wire.view = ged_draw_active_view_ctx(gedp);
    redraw_wire.mode = GED_DRAW_MODE_WIRE;
    if (ged_draw_obol_scene_sync_transaction(gedp, &redraw_wire, NULL,
	    scene) <= 0 || frontier_visible_occurrences(shaded) != 0)
	FAIL("mode-scoped redraw changed a sibling representation");
    if (shaded->setCompactInstanceVisibilityOverrideForPathMatch(root_path,
	    BOBOL_COMPACT_PATH_SUBTREE, TRUE) <= 0)
	FAIL("mode-scoped redraw fixture sibling restore");
    const SbString shaded_key = shaded->instanceKey.getValue();
    struct ged_scene_reducer_request erase_wire =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_ERASE, root_path);
    erase_wire.view = ged_draw_active_view_ctx(gedp);
    erase_wire.mode = GED_DRAW_MODE_WIRE;
    if (ged_draw_obol_scene_sync_transaction(gedp, &erase_wire, NULL,
	    scene) <= 0 || source_for_representation(scene, root_path,
		SoBRLDatabaseSource::REPRESENTATION_WIRE) ||
	scene->findDatabaseSourceInstance(shaded_key.getString()) != shaded ||
	frontier_visible_occurrences(shaded) != shaded->getCompactInstanceCount())
	FAIL("mode-scoped erase hid the surviving compact representation");

    (void)scene->clearDatabaseSources();
    frontier_clear(gedp);
    constexpr const char *nested_target =
	"nested_parent.c/nested_child.c/nested_leaf.s";
    if (ged_exec_draw(gedp, 2, draw) != BRLCAD_OK)
	FAIL("erase acceptance fixture broad draw");
    (void)scene->realizePending();
    SoBRLDatabaseSource *broad = source_for_representation(scene, root_path,
	SoBRLDatabaseSource::REPRESENTATION_WIRE);
    if (!broad)
	FAIL("erase acceptance fixture broad source");
    struct ged_bobol_publication_context narrow_publication;
    if (!ged_bobol_publication_begin(&narrow_publication, gedp,
	    ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_WIRE))
	FAIL("erase acceptance fixture narrow publication");
    if (!ged_bobol_database_source_ensure_for_path(&narrow_publication,
	    nested_target, gedp->dbip, GED_DRAW_MODE_WIRE,
	    broad->sourceRevision.getValue())) {
	ged_bobol_publication_end(&narrow_publication);
	FAIL("erase acceptance fixture narrow draw");
    }
    ged_bobol_publication_end(&narrow_publication);
    SoBRLDatabaseSource *narrow = source_for_path(scene, nested_target);
    if (!narrow || narrow == broad ||
	broad->clearCompactInstanceVisibilityFrontier() < 0 ||
	broad->setCompactInstanceVisibilityOverrideForPathMatch(nested_target,
	    BOBOL_COMPACT_PATH_SUBTREE, FALSE) <= 0)
	FAIL("erase acceptance fixture retained compact owner");
    const SbString narrow_key = narrow->instanceKey.getValue();
    SbModernUtils::SoNodeRef broad_owner(broad), narrow_owner(narrow);
    struct EraseObserver {
	SoBRLDatabaseSource &broad;
	const char *target;
	bool called = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<EraseObserver *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    self.failed =
		self.broad.setCompactInstanceVisibilityOverrideForPathMatch(
		    self.target, BOBOL_COMPACT_PATH_SUBTREE, TRUE) <= 0;
	}
    } erase_observer{*broad, nested_target};
    SoNodeSensor erase_sensor(EraseObserver::changed, &erase_observer);
    erase_sensor.setPriority(0);
    erase_sensor.attach(scene->getSceneRoot());
    struct ged_scene_reducer_request erase_nested =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_ERASE,
	    nested_target);
    erase_nested.view = ged_draw_active_view_ctx(gedp);
    erase_nested.mode = GED_DRAW_MODE_WIRE;
    const int erased_nested = ged_draw_obol_scene_sync_transaction(gedp,
	&erase_nested, NULL, scene);
    erase_sensor.detach();
    BObolCompactInstanceHandle target_handle;
    BObolCompactInstanceSummary target_summary;
    if (erased_nested <= 0 || !erase_observer.called ||
	erase_observer.failed ||
	scene->findDatabaseSourceInstance(narrow_key.getString()) ||
	source_for_representation(scene, root_path,
	    SoBRLDatabaseSource::REPRESENTATION_WIRE) != broad ||
	!broad->getCompactInstanceForPath(nested_target, FALSE, FALSE,
	    target_handle, target_summary) || !target_summary.valid ||
	!target_summary.visible)
	FAIL("erase overwrote a newer callback edit on a retained compact owner");

    (void)scene->clearDatabaseSources();
    frontier_clear(gedp);
    if (ged_exec_draw(gedp, 2, draw) != BRLCAD_OK ||
	ged_exec_draw(gedp, 2, draw_second) != BRLCAD_OK)
	FAIL("erase-prefix realization acceptance fixture draws");
    (void)scene->realizePending();
    SoBRLDatabaseSource *retained = source_for_representation(scene,
	"reuse_root.c", SoBRLDatabaseSource::REPRESENTATION_WIRE);
    if (!retained || retained->needsRealization())
	FAIL("erase-prefix realization acceptance fixture current source");
    const SbString retained_key = retained->instanceKey.getValue();
    SbModernUtils::SoNodeRef retained_owner(retained);
    struct ErasePrefixObserver {
	BObolSceneController &scene;
	SbString retainedKey;
	bool called = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<ErasePrefixObserver *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    self.failed = self.scene.markDatabaseSourceInstanceStale(
		self.retainedKey.getString(),
		SoBRLDatabaseSource::STALE_SOURCE) <= 0;
	}
    } erase_prefix_observer{*scene, retained_key};
    SoNodeSensor erase_prefix_sensor(ErasePrefixObserver::changed,
	&erase_prefix_observer);
    erase_prefix_sensor.setPriority(0);
    erase_prefix_sensor.attach(scene->getSceneRoot());
    struct ged_scene_reducer_request erase_prefix =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_ERASE_PREFIX,
	    root_path);
    erase_prefix.view = ged_draw_active_view_ctx(gedp);
    const int erased_prefix = ged_draw_obol_scene_sync_transaction(gedp,
	&erase_prefix, NULL, scene);
    erase_prefix_sensor.detach();
    if (erased_prefix <= 0 || !erase_prefix_observer.called ||
	erase_prefix_observer.failed || source_for_path(scene, root_path) ||
	scene->findDatabaseSourceInstance(retained_key.getString()) != retained ||
	!retained->needsRealization())
	FAIL("erase-prefix started unrelated realization created by its callback");

    (void)scene->clearDatabaseSources();
    frontier_clear(gedp);
    std::puts("PASS GED redraw/erase publication: compact state, hierarchy and callback traversal remain current");
    return 0;
}

static int
exercise_base_source_promotion(struct ged *gedp, BObolSceneController *scene)
{
    if (!gedp || !gedp->dbip || !scene)
	FAIL("base-source promotion test needs GED and Obol scene state");

    const char *path = "box.s";
    frontier_clear(gedp);
    if (!ged_draw_obol_database_source_ensure_for_path(gedp, path,
	    gedp->dbip, GED_DRAW_MODE_SHADED, 41))
	FAIL("base-source promotion setup should publish a default source");

    SoBRLDatabaseSource *source = source_for_path(scene, path);
    BObolDatabaseSourceSummary before;
    if (!source || !source->getSummary(before) || !before.valid ||
	before.representationMode !=
	    SoBRLDatabaseSource::REPRESENTATION_DEFAULT)
	FAIL("base-source promotion setup should retain a default representation");
    SbModernUtils::SoNodeRef source_owner(source);

    constexpr int callback_line_width = 19;
    struct Observer {
	BObolSceneController &scene;
	SoBRLDatabaseSource &source;
	SbString originalKey;
	bool called = false;
	bool partial = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<Observer *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    BObolDatabaseSourceSummary observed;
	    const char *current_key =
		self.source.instanceKey.getValue().getString();
	    self.partial = !self.source.getSummary(observed) ||
		!observed.valid ||
		observed.representationMode !=
		    SoBRLDatabaseSource::REPRESENTATION_SHADED ||
		!current_key || !current_key[0] ||
		BU_STR_EQUAL(current_key, self.originalKey.getString()) ||
		self.scene.findDatabaseSourceInstance(current_key) !=
		    &self.source;
	    if (self.partial)
		return;
	    BObolDatabaseSourceDisplayPatch patch;
	    patch.lineWidthValid = TRUE;
	    patch.lineWidth = callback_line_width;
	    self.failed =
		self.scene.setDatabaseSourceInstanceDisplayPatch(
		    current_key, patch) < 0;
	}
    } observer{*scene, *source, before.instanceKey};

    SoNodeSensor sensor(Observer::changed, &observer);
    sensor.setPriority(0);
    sensor.attach(source);
    struct ged_bobol_publication_context publication;
    const bool begun = ged_bobol_publication_begin(&publication, gedp,
	ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_SHADED);
    const int drawn = begun ? ged_bobol_database_source_ensure_for_path(
	&publication, path, gedp->dbip, GED_DRAW_MODE_SHADED, 42) : 0;
    if (begun)
	ged_bobol_publication_end(&publication);
    sensor.detach();

    BObolDatabaseSourceSummary after;
    SoBRLDatabaseSource *current = source_for_representation(scene, path,
	SoBRLDatabaseSource::REPRESENTATION_SHADED);
    if (drawn <= 0 || !observer.called || observer.partial ||
	observer.failed || current != source ||
	!source->getSummary(after) || !after.valid ||
	after.lineWidth != callback_line_width ||
	scene->findDatabaseSourceInstance(before.instanceKey.getString())) {
	fprintf(stderr, "base promotion drawn=%d called=%d partial=%d "
	    "failed=%d original=%s current=%s same=%d rep=%d width=%d "
	    "old-present=%d\n", drawn, observer.called, observer.partial,
	    observer.failed, before.instanceKey.getString(),
	    after.instanceKey.getString(), current == source,
	    after.representationMode, after.lineWidth,
	    scene->findDatabaseSourceInstance(
		before.instanceKey.getString()) != NULL);
	FAIL("base-source promotion must publish complete state and preserve a newer callback edit");
    }

    frontier_clear(gedp);
    std::puts("PASS GED base-source promotion: one complete publication and current callback target");
    return 0;
}

static int
exercise_deferred_proxy_transition(struct ged *gedp,
	BObolSceneController *scene)
{
    if (!gedp || !gedp->dbip || !scene)
	FAIL("deferred-proxy transition test needs GED and Obol scene state");

    constexpr const char *path = "box.s";
    constexpr int callback_line_width = 23;
    frontier_clear(gedp);
    (void)scene->clearDatabaseSources();

    struct ged_bobol_publication_context publication;
    if (!ged_bobol_publication_begin(&publication, gedp,
	ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_SHADED))
	FAIL("deferred-proxy transition should start a GED publication");
    const int created = ged_bobol_database_source_ensure_for_path(
	&publication, path, gedp->dbip, GED_DRAW_MODE_SHADED, 51);
    ged_bobol_publication_end(&publication);
    SoBRLDatabaseSource *source = source_for_representation(scene, path,
	SoBRLDatabaseSource::REPRESENTATION_SHADED);
    if (created <= 0 || !source)
	FAIL("deferred-proxy transition should create its source");

    SbVec3f points[2] = {
	SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 1.0f, 1.0f)
    };
    int32_t commands[2] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW
    };
    BObolExternalLineSet proxy;
    proxy.points = points;
    proxy.commands = commands;
    proxy.count = 2;
    proxy.sourceType = "proxy";
    proxy.geometryKind = "aabb";
    const SbString instance_key = source->instanceKey.getValue();
    if (scene->publishDatabaseSourceInstanceExternalLineSet(
	    instance_key.getString(), proxy) <= 0)
	FAIL("deferred-proxy transition should publish its proxy");

    BObolDatabaseSourceSummary before;
    if (!source->getSummary(before) || !before.valid ||
	before.realizationStatus != SoBRLDatabaseSource::REALIZED ||
	!(before.realizationRoleFlags &
	    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL) ||
	before.realizedShapeCount != 1)
	FAIL("deferred-proxy transition needs one current external proxy");
    SbModernUtils::SoNodeRef source_owner(source);

    struct Observer {
	BObolSceneController &scene;
	SoBRLDatabaseSource &source;
	SbString instanceKey;
	uint64_t frame;
	bool called = false;
	bool partial = false;
	bool failed = false;
	bool retainedProxy = false;
	int editedWidth = -1;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<Observer *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    BObolDatabaseSourceSummary observed;
	    self.partial = !self.source.getSummary(observed) ||
		!observed.valid || !observed.stale ||
		observed.realizationStatus !=
		    SoBRLDatabaseSource::UNREALIZED ||
		!(observed.staleReason &
		    SoBRLDatabaseSource::STALE_SOURCE) ||
		(observed.realizationRoleFlags &
		    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL) ||
		self.scene.findDatabaseSourceInstance(
		    self.instanceKey.getString()) != &self.source;
	    self.retainedProxy = self.source.getRealizedShape() != NULL;
	    self.partial = self.partial || !self.retainedProxy;
	    BObolDatabaseSourceDisplayPatch patch;
	    patch.lineWidthValid = TRUE;
	    patch.lineWidth = callback_line_width;
	    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(
		self.instanceKey.getString(), patch) <= 0;
	    self.editedWidth = self.source.lineWidth.getValue();
	}
    } observer{*scene, *source, instance_key, scene->getFrameRevision()};

    SoNodeSensor sensor(Observer::changed, &observer);
    sensor.setPriority(0);
    sensor.attach(source);
    if (!ged_bobol_publication_begin(&publication, gedp,
	ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_SHADED))
	FAIL("deferred-proxy transition should restart its GED publication");
    const int changed = ged_bobol_database_source_ensure_for_path(
	&publication, path, gedp->dbip, GED_DRAW_MODE_SHADED, 51);
    ged_bobol_publication_end(&publication);
    sensor.detach();

    BObolDatabaseSourceSummary after;
    SoBRLDatabaseSource *current =
	scene->findDatabaseSourceInstance(instance_key.getString());
    const bool after_valid = source->getSummary(after) && after.valid;
    if (changed <= 0 || !observer.called || observer.partial ||
	observer.failed || current != source ||
	!after_valid || !after.stale ||
	after.realizationStatus != SoBRLDatabaseSource::UNREALIZED ||
	(after.realizationRoleFlags &
	    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL) ||
	!source->getRealizedShape() ||
	after.lineWidth != callback_line_width ||
	scene->getFrameRevision() != observer.frame + 1) {
	fprintf(stderr, "deferred proxy changed=%d called=%d partial=%d "
	    "failed=%d retained=%d edited-width=%d same=%d stale=%d "
	    "status=%d roles=%d shapes=%d children=%d width=%d "
	    "frame=%" PRIu64 " expected-frame=%" PRIu64 "\n",
	    changed, observer.called, observer.partial, observer.failed,
	    observer.retainedProxy, observer.editedWidth, current == source,
	    after.stale, after.realizationStatus, after.realizationRoleFlags,
	    after.realizedShapeCount, source->getNumChildren(),
	    after.lineWidth, scene->getFrameRevision(), observer.frame + 1);
	FAIL("deferred proxy must remain visible through one complete invalidation publication");
    }

    /* Exercise the real deferred draw adapter with this standing external
     * proxy. It may notify the preceding AABB publication first; the compact
     * transition callback must see the complete replacement and may publish a
     * newer display edit without that edit being overwritten on unwind. */
    constexpr int compact_callback_line_width = 27;
    struct CompactObserver {
	BObolSceneController &scene;
	const char *path;
	SoBRLDatabaseSource *source = nullptr;
	bool called = false;
	bool partial = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<CompactObserver *>(data);
	    SoBRLDatabaseSource *current = source_for_path(
		&self.scene, self.path);
	    if (self.called || !current ||
		!current->isCompactOccurrenceRegistry() ||
		current->getCompactInstanceCount() <= 0)
		return;
	    self.called = true;
	    self.source = current;
	    BObolDatabaseSourceSummary observed;
	    SbBox3f bounds;
	    self.partial = !current->getSummary(observed) ||
		!observed.valid ||
		observed.realizationStatus != SoBRLDatabaseSource::REALIZED ||
		observed.stale ||
		!(observed.realizationRoleFlags &
		    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL) ||
		observed.realizedShapeCount != 0 ||
		current->getRealizedShape() != NULL ||
		!current->getSourceBounds(bounds) || bounds.isEmpty() ||
		self.scene.findDatabaseSourceInstance(
		    observed.instanceKey.getString()) != current;
	    if (self.partial)
		return;
	    BObolDatabaseSourceDisplayPatch patch;
	    patch.lineWidthValid = TRUE;
	    patch.lineWidth = compact_callback_line_width;
	    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(
		observed.instanceKey.getString(), patch) <= 0;
	    throw std::bad_alloc();
	}
    } compactObserver{*scene, path};
    SoNodeSensor compactSensor(CompactObserver::changed, &compactObserver);
    compactSensor.setPriority(0);
    compactSensor.attach(scene->getSceneRoot());
    struct ged_draw_appearance_settings appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    appearance.defer_leaf_expansion = 1;
    struct ged_scene_reducer_request draw =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_DRAW, path);
    draw.view = ged_draw_active_view_ctx(gedp);
    draw.mode = GED_DRAW_MODE_SHADED;
    draw.appearance = &appearance;
    struct ged_scene_reducer_result result;
    ged_scene_reducer_result_init(&result);
    int drawn = 0;
    bool escaped = false;
    const uint64_t compactFrame = scene->getFrameRevision();
    {
	ScopedTransactionFault fault(
	    BObolTransactionFaultPoint::SOURCE_REALIZATION_WORKER_START);
	try {
	    drawn = ged_scene_reduce(gedp, &draw, &result);
	} catch (const std::exception &) {
	    escaped = true;
	}
    }
    ged_scene_reducer_result_free(&result);
    compactSensor.detach();
    current = source_for_path(scene, path);
    if (escaped || drawn <= 0 || !compactObserver.called ||
	compactObserver.partial || compactObserver.failed ||
	current != compactObserver.source || !current ||
	current->lineWidth.getValue() != compact_callback_line_width ||
	scene->getFrameRevision() <= compactFrame) {
	fprintf(stderr, "deferred compact drawn=%d escaped=%d called=%d "
	    "partial=%d failed=%d current=%p source=%p compact=%d count=%d "
	    "width=%d frame=%" PRIu64 "/%" PRIu64 "\n",
	    drawn, escaped, compactObserver.called, compactObserver.partial,
	    compactObserver.failed, (void *)current, (void *)source,
	    current ? current->isCompactOccurrenceRegistry() : 0,
	    current ? current->getCompactInstanceCount() : 0,
	    current ? current->lineWidth.getValue() : -1,
	    scene->getFrameRevision(), compactFrame);
	for (int i = 0; i < scene->getDatabaseSourceCount(); ++i) {
	    BObolDatabaseSourceSummary candidate;
	    SoBRLDatabaseSource *candidateSource = scene->getDatabaseSource(i);
	    if (candidateSource && candidateSource->getSummary(candidate))
		fprintf(stderr, "  source[%d] path=%s key=%s rep=%d status=%d "
		    "stale=%d roles=%d shapes=%d direct=%p exact=%d compact=%d "
		    "count=%d width=%d\n", i,
		    candidate.path.getString(), candidate.instanceKey.getString(),
		    candidate.representationMode, candidate.realizationStatus,
		    candidate.stale, candidate.realizationRoleFlags,
		    candidate.realizedShapeCount,
		    (void *)candidateSource->getRealizedShape(),
		    candidateSource->hasExactSourceBounds(),
		    candidateSource->isCompactOccurrenceRegistry(),
		    candidateSource->getCompactInstanceCount(), candidate.lineWidth);
	}
	FAIL("deferred compact proxy must publish complete state and retain its callback edit");
    }

    (void)scene->removeDatabaseSourceInstance(
	current->instanceKey.getValue().getString());
    frontier_clear(gedp);
    std::puts("PASS GED deferred-proxy transition: complete invalidation and compact replacement callbacks");
    return 0;
}

static int
exercise_retained_current_publication(struct ged *gedp,
	BObolSceneController *scene)
{
    if (!gedp || !gedp->dbip || !scene)
	FAIL("retained-current test needs GED and Obol scene state");

    enum class CallbackEdit { Change, Replace, Remove };
    constexpr const char *path = "box.s";
    constexpr int changed_line_width = 29;
    constexpr int replacement_line_width = 31;
    constexpr float policy_epsilon = 0.000001f;

    struct ged_view_context *view_ctx = ged_draw_active_view_ctx(gedp);
    ged_view_lod_policy original_policy = BV_LOD_POLICY_INIT;
    if (!view_ctx || !ged_view_lod_policy_get(&original_policy, view_ctx))
	FAIL("retained-current test needs its view policy");
    ged_view_lod_policy publication_policy = original_policy;
    publication_policy.csg_enabled = 1;
    publication_policy.mesh_enabled = 0;
    if (!ged_view_lod_policy_apply(view_ctx, &publication_policy))
	FAIL("retained-current test should enable a view-dependent policy");

    for (CallbackEdit edit : {CallbackEdit::Change,
	CallbackEdit::Replace, CallbackEdit::Remove}) {
	frontier_clear(gedp);
	(void)scene->clearDatabaseSources();

	struct ged_bobol_publication_context publication;
	if (!ged_bobol_publication_begin(&publication, gedp,
		view_ctx, GED_DRAW_MODE_WIRE))
	    FAIL("retained-current test should start a GED publication");
	const int created = ged_bobol_database_source_ensure_for_path(
	    &publication, path, gedp->dbip, GED_DRAW_MODE_WIRE, 61);
	ged_bobol_publication_end(&publication);
	SoBRLDatabaseSource *source = source_for_representation(scene, path,
	    SoBRLDatabaseSource::REPRESENTATION_WIRE);
	if (created <= 0 || !source)
	    FAIL("retained-current test should create its source");

	SbVec3f points[2] = {
	    SbVec3f(0.0f, 0.0f, 0.0f),
	    SbVec3f(1.0f, 1.0f, 1.0f)
	};
	int32_t commands[2] = {
	    SoBRLVListShape::MOVE,
	    SoBRLVListShape::DRAW
	};
	BObolExternalLineSet proxy;
	proxy.points = points;
	proxy.commands = commands;
	proxy.count = 2;
	proxy.sourceType = "retained";
	proxy.geometryKind = "wire";
	const SbString instance_key = source->instanceKey.getValue();
	if (scene->publishDatabaseSourceInstanceExternalLineSet(
		instance_key.getString(), proxy) <= 0)
	    FAIL("retained-current test should publish its retained geometry");

	BObolDatabaseSourceSummary expected_policy;
	if (!source->getSummary(expected_policy) || !expected_policy.valid ||
	    !expected_policy.realizationViewDependent ||
	    !expected_policy.realizationCsgLodEnabled ||
	    expected_policy.realizationMeshLodEnabled)
	    FAIL("retained-current test should publish its target view policy");
	if (scene->setDatabaseSourceInstanceRealizationViewPolicy(
		instance_key.getString(), FALSE, FALSE, FALSE,
		0.0f, 1.0f, 0, 0, 0, 0.0f, 0.0f) <= 0)
	    FAIL("retained-current test should prepare a mismatched view policy");
	BObolDatabaseSourceSummary stale_summary;
	if (!source->getSummary(stale_summary) || !stale_summary.valid ||
	    scene->setDatabaseSourceInstanceRealizationState(
		instance_key.getString(), SoBRLDatabaseSource::REALIZED,
		stale_summary.sourceRevision, stale_summary.inputsRevision,
		SoBRLDatabaseSource::STALE_NONE, NULL,
		stale_summary.realizationRoleFlags) <= 0)
	    FAIL("retained-current test should make its mismatched source current");
	BObolDatabaseSourceSummary before;
	SoNode *retained_shape = source->getRealizedShape();
	if (!source->getSummary(before) || !before.valid || before.stale ||
	    before.realizationStatus != SoBRLDatabaseSource::REALIZED ||
	    !(before.realizationRoleFlags &
		SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL) ||
	    before.realizationViewDependent ||
	    !retained_shape ||
	    !retained_shape->isOfType(SoBRLVListShape::getClassTypeId()))
	    FAIL("retained-current test needs current retained line geometry with an old policy");
	SbModernUtils::SoNodeRef source_owner(source);

	struct Observer {
	    BObolSceneController &scene;
	    SoBRLDatabaseSource &source;
	    BObolDatabaseSourceSummary before;
	    BObolDatabaseSourceSummary expectedPolicy;
	    struct db_i *database;
	    CallbackEdit edit;
	    float policyEpsilon;
	    bool called = false;
	    bool partial = false;
	    bool failed = false;
	    static void changed(void *data, SoSensor *)
	    {
		auto &self = *static_cast<Observer *>(data);
		if (self.called)
		    return;
		self.called = true;
		BObolDatabaseSourceSummary observed;
		SoNode *shape_node = self.source.getRealizedShape();
		auto *shape = shape_node && shape_node->isOfType(
		    SoBRLVListShape::getClassTypeId()) ?
		    static_cast<SoBRLVListShape *>(shape_node) : NULL;
		self.partial = !self.source.getSummary(observed) ||
		    !observed.valid || observed.stale ||
		    observed.staleReason != SoBRLDatabaseSource::STALE_NONE ||
		    observed.realizationStatus != SoBRLDatabaseSource::REALIZED ||
		    !(observed.realizationRoleFlags &
			SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL) ||
		    observed.realizedSourceRevision != observed.sourceRevision ||
		    observed.realizedInputsRevision != observed.inputsRevision ||
		    observed.realizationViewDependent !=
			self.expectedPolicy.realizationViewDependent ||
		    observed.realizationCsgLodEnabled !=
			self.expectedPolicy.realizationCsgLodEnabled ||
		    observed.realizationMeshLodEnabled !=
			self.expectedPolicy.realizationMeshLodEnabled ||
		    std::abs(observed.realizationViewScale -
			self.expectedPolicy.realizationViewScale) >
			self.policyEpsilon ||
		    std::abs(observed.realizationLodScale -
			self.expectedPolicy.realizationLodScale) >
			self.policyEpsilon ||
		    observed.realizationViewWidth !=
			self.expectedPolicy.realizationViewWidth ||
		    observed.realizationViewHeight !=
			self.expectedPolicy.realizationViewHeight ||
		    observed.realizationBotThreshold !=
			self.expectedPolicy.realizationBotThreshold ||
		    std::abs(observed.realizationCurveScale -
			self.expectedPolicy.realizationCurveScale) >
			self.policyEpsilon ||
		    std::abs(observed.realizationPointScale -
			self.expectedPolicy.realizationPointScale) >
			self.policyEpsilon ||
		    self.scene.findDatabaseSourceInstance(
			self.before.instanceKey.getString()) != &self.source ||
		    !shape ||
		    shape->ownerRealizationStatus.getValue() !=
			SoBRLDatabaseSource::REALIZED ||
		    shape->ownerSourceStale.getValue();
		if (self.partial)
		    return;

		const char *key = self.before.instanceKey.getString();
		if (self.edit == CallbackEdit::Change) {
		    BObolDatabaseSourceDisplayPatch patch;
		    patch.lineWidthValid = TRUE;
		    patch.lineWidth = changed_line_width;
		    self.failed =
			self.scene.setDatabaseSourceInstanceDisplayPatch(
			    key, patch) <= 0;
		    return;
		}

		self.failed = self.scene.removeDatabaseSourceInstance(key) <= 0;
		if (self.edit == CallbackEdit::Remove)
		    return;
		BObolDatabaseSourcePublishState replacement;
		replacement.sourceInstanceKey = key;
		replacement.sourcePath = self.before.path.getString();
		replacement.sourceRepresentationKey =
		    self.before.representationKey.getString();
		replacement.database = self.database;
		replacement.drawMode = self.before.drawMode;
		replacement.representationMode = self.before.representationMode;
		replacement.sourceRevisionValid = TRUE;
		replacement.sourceRevision = self.before.sourceRevision + 1;
		replacement.lineWidth = replacement_line_width;
		replacement.roleFlagsValid = TRUE;
		replacement.roleFlags =
		    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL;
		if (self.scene.publishDatabaseSourceInstance(replacement) < 0) {
		    self.failed = true;
		    return;
		}
		const SbVec3f replacement_points[2] = {
		    SbVec3f(0.0f, 0.0f, 0.0f),
		    SbVec3f(2.0f, 2.0f, 2.0f)
		};
		const int32_t replacement_commands[2] = {
		    SoBRLVListShape::MOVE,
		    SoBRLVListShape::DRAW
		};
		BObolExternalLineSet replacement_geometry;
		replacement_geometry.points = replacement_points;
		replacement_geometry.commands = replacement_commands;
		replacement_geometry.count = 2;
		replacement_geometry.sourceType = "callback";
		replacement_geometry.geometryKind = "wire";
		if (self.scene.publishDatabaseSourceInstanceExternalLineSet(
			key, replacement_geometry) <= 0) {
		    self.failed = true;
		    return;
		}
		BObolDatabaseSourceSummary replacement_summary;
		SoBRLDatabaseSource *replacement_source =
		    self.scene.findDatabaseSourceInstance(key);
		if (!replacement_source ||
		    !replacement_source->getSummary(replacement_summary) ||
		    !replacement_summary.valid ||
		    self.scene.setDatabaseSourceInstanceRealizationState(
			key, SoBRLDatabaseSource::UNREALIZED,
			replacement_summary.sourceRevision,
			replacement_summary.inputsRevision,
			SoBRLDatabaseSource::STALE_SOURCE, NULL,
			replacement_summary.realizationRoleFlags) <= 0)
		    self.failed = true;
	    }
	} observer{*scene, *source, before, expected_policy, gedp->dbip, edit,
	    policy_epsilon};

	SoNodeSensor sensor(Observer::changed, &observer);
	sensor.setPriority(0);
	sensor.attach(source);
	const uint64_t frame = scene->getFrameRevision();
	if (!ged_bobol_publication_begin(&publication, gedp, view_ctx,
		GED_DRAW_MODE_WIRE))
	    FAIL("retained-current test should restart its GED publication");
	const int changed = ged_bobol_database_source_ensure_for_path(
	    &publication, path, gedp->dbip, GED_DRAW_MODE_WIRE, 61);
	ged_bobol_publication_end(&publication);
	sensor.detach();

	SoBRLDatabaseSource *current =
	    scene->findDatabaseSourceInstance(instance_key.getString());
	bool final_state_valid = changed > 0 && observer.called &&
	    !observer.partial && !observer.failed &&
	    scene->getFrameRevision() == frame + 1;
	if (edit == CallbackEdit::Change) {
	    final_state_valid = final_state_valid && current == source &&
		source->lineWidth.getValue() == changed_line_width &&
		source->realizationStatus.getValue() ==
		    SoBRLDatabaseSource::REALIZED &&
		(source->realizationRoleFlags.getValue() &
		    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL);
	} else if (edit == CallbackEdit::Replace) {
	    final_state_valid = final_state_valid && current &&
		current != source &&
		current->lineWidth.getValue() == replacement_line_width &&
		current->sourceRevision.getValue() == before.sourceRevision + 1 &&
		current->realizationStatus.getValue() ==
		    SoBRLDatabaseSource::UNREALIZED && current->stale.getValue() &&
		(current->realizationRoleFlags.getValue() &
		    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL) &&
		current->getRealizedShape();
	} else {
	    final_state_valid = final_state_valid && !current;
	}
	if (!final_state_valid) {
	    fprintf(stderr, "retained current edit=%d changed=%d called=%d "
		"partial=%d failed=%d current=%p original=%p width=%d "
		"status=%d stale=%d roles=%d frame=%" PRIu64 "\n",
		int(edit), changed, observer.called, observer.partial,
		observer.failed, (void *)current, (void *)source,
		current ? current->lineWidth.getValue() : -1,
		current ? current->realizationStatus.getValue() : -1,
		current ? current->stale.getValue() : -1,
		current ? current->realizationRoleFlags.getValue() : -1,
		scene->getFrameRevision());
	    FAIL("retained-current publication must commit status and ownership before observers");
	}
	(void)scene->removeDatabaseSourceInstance(instance_key.getString());
	std::printf("PASS GED retained-current sequence edit%d: complete state and current callback target\n",
	    int(edit));
    }

    frontier_clear(gedp);
    if (!ged_view_lod_policy_apply(view_ctx, &original_policy))
	FAIL("retained-current test should restore its view policy");
    return 0;
}

static int
exercise_source_summary_copy_publication(struct ged *gedp,
	BObolSceneController *primary_scene)
{
    if (!gedp || !gedp->dbip || !primary_scene)
	FAIL("source-copy test needs GED and its primary scene");

    enum class CallbackEdit { Change, Replace, Remove };
    constexpr const char *path = "box.s";
    constexpr uint64_t primary_source_revision = 67;
    constexpr int changed_line_width = 33;
    constexpr int replacement_line_width = 35;

    frontier_clear(gedp);
    (void)primary_scene->clearDatabaseSources();
    struct ged_bobol_publication_context publication;
    if (!ged_bobol_publication_begin(&publication, gedp,
	    ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_WIRE))
	FAIL("source-copy test should start its primary publication");
    const int created = ged_bobol_database_source_ensure_for_path(
	&publication, path, gedp->dbip, GED_DRAW_MODE_WIRE,
	primary_source_revision);
    ged_bobol_publication_end(&publication);
    SoBRLDatabaseSource *primary = source_for_representation(primary_scene,
	path, SoBRLDatabaseSource::REPRESENTATION_WIRE);
    if (created <= 0 || !primary)
	FAIL("source-copy test should create its primary source");
    const char *primary_key = primary->instanceKey.getValue().getString();
    if (!primary_key || !primary_key[0] ||
	primary_scene->setDatabaseSourceInstanceMaterialPolicy(primary_key,
	    SoBRLDatabaseSource::MATERIAL_DATABASE) < 0)
	FAIL("source-copy test should set its primary material policy");
    BObolDatabaseSourceSummary expected;
    if (!primary->getSummary(expected) || !expected.valid ||
	expected.materialPolicy != SoBRLDatabaseSource::MATERIAL_DATABASE)
	FAIL("source-copy test needs a database material policy");

    for (CallbackEdit edit : {CallbackEdit::Change,
	CallbackEdit::Replace, CallbackEdit::Remove}) {
	SoSeparator *root = new SoSeparator;
	root->ref();
	{
	    BObolSceneController destination(root);
	    struct Observer {
		BObolSceneController &scene;
		const BObolDatabaseSourceSummary &expected;
		struct db_i *database;
		CallbackEdit edit;
		bool called = false;
		bool partial = false;
		bool failed = false;
		static void changed(void *data, SoSensor *)
		{
		    auto &self = *static_cast<Observer *>(data);
		    if (self.called)
			return;
		    self.called = true;
		    const char *key = self.expected.instanceKey.getString();
		    SoBRLDatabaseSource *source =
			self.scene.findDatabaseSourceInstance(key);
		    BObolDatabaseSourceSummary observed;
		    self.partial = !source || !source->getSummary(observed) ||
			!observed.valid ||
			observed.path != self.expected.path ||
			observed.representationKey !=
			    self.expected.representationKey ||
			observed.representationMode !=
			    self.expected.representationMode ||
			observed.sourceRevision !=
			    self.expected.sourceRevision ||
			observed.inputsRevision !=
			    self.expected.inputsRevision ||
			observed.materialPolicy !=
			    self.expected.materialPolicy ||
			observed.realizationRoleFlags !=
			    self.expected.realizationRoleFlags ||
			self.scene.findDatabaseSourceInstance(key) != source;
		    if (self.partial)
			return;

		    if (self.edit == CallbackEdit::Change) {
			BObolDatabaseSourceDisplayPatch patch;
			patch.lineWidthValid = TRUE;
			patch.lineWidth = changed_line_width;
			self.failed =
			    self.scene.setDatabaseSourceInstanceDisplayPatch(
				key, patch) <= 0;
			return;
		    }

		    self.failed = self.scene.removeDatabaseSourceInstance(key) <= 0;
		    if (self.failed || self.edit == CallbackEdit::Remove)
			return;
		    BObolDatabaseSourcePublishState replacement;
		    replacement.sourceInstanceKey = key;
		    replacement.sourcePath = self.expected.path.getString();
		    replacement.sourceRepresentationKey =
			self.expected.representationKey.getString();
		    replacement.database = self.database;
		    replacement.drawMode = self.expected.drawMode;
		    replacement.representationMode =
			self.expected.representationMode;
		    replacement.sourceRevisionValid = TRUE;
		    replacement.sourceRevision = self.expected.sourceRevision + 1;
		    replacement.lineWidth = replacement_line_width;
		    replacement.materialPolicyValid = TRUE;
		    replacement.materialPolicy =
			SoBRLDatabaseSource::MATERIAL_INHERIT;
		    replacement.roleFlagsValid = TRUE;
		    replacement.roleFlags =
			SoBRLDatabaseSource::REALIZATION_ROLE_NONE;
		    self.failed =
			self.scene.publishDatabaseSourceInstance(replacement) < 0;
		}
	    } observer{destination, expected, gedp->dbip, edit};

	    SoNodeSensor sensor(Observer::changed, &observer);
	    sensor.setPriority(0);
	    sensor.attach(root);
	    struct ged_scene_reducer_request request =
		ged_scene_reducer_request_make(GED_SCENE_REDUCER_DRAW, path);
	    request.mode = GED_DRAW_MODE_WIRE;
	    request.redraw = 1;
	    const uint64_t frame = destination.getFrameRevision();
	    const int changed = ged_draw_obol_scene_sync_transaction(gedp,
		&request, NULL, &destination);
	    sensor.detach();

	    SoBRLDatabaseSource *current =
		destination.findDatabaseSourceInstance(
		    expected.instanceKey.getString());
	    bool final_state_valid = changed > 0 && observer.called &&
		!observer.partial && !observer.failed &&
		destination.getFrameRevision() == frame + 1;
	    if (edit == CallbackEdit::Change) {
		final_state_valid = final_state_valid && current &&
		    current->lineWidth.getValue() == changed_line_width &&
		    current->materialPolicy.getValue() ==
			SoBRLDatabaseSource::MATERIAL_DATABASE;
	    } else if (edit == CallbackEdit::Replace) {
		final_state_valid = final_state_valid && current &&
		    current->lineWidth.getValue() == replacement_line_width &&
		    current->sourceRevision.getValue() ==
			expected.sourceRevision + 1 &&
		    current->materialPolicy.getValue() ==
			SoBRLDatabaseSource::MATERIAL_INHERIT;
	    } else {
		final_state_valid = final_state_valid && !current;
	    }
	    if (!final_state_valid) {
		fprintf(stderr, "source copy edit=%d changed=%d called=%d "
		    "partial=%d failed=%d current=%p policy=%d width=%d "
		    "status=%d stale=%d frame=%" PRIu64 "/%" PRIu64 "\n",
		    int(edit), changed, observer.called, observer.partial,
		    observer.failed, (void *)current,
		    current ? current->materialPolicy.getValue() : -1,
		    current ? current->lineWidth.getValue() : -1,
		    current ? current->realizationStatus.getValue() : -1,
		    current ? current->stale.getValue() : -1,
		    destination.getFrameRevision(), frame + 1);
		FAIL("source copy must publish material policy with the complete source");
	    }
	}
	root->unref();
	std::printf("PASS GED source-copy publication edit%d: complete state and current callback target\n",
	    int(edit));
    }

    (void)primary_scene->removeDatabaseSourceInstance(
	expected.instanceKey.getString());
    frontier_clear(gedp);
    return 0;
}

static int
exercise_presentation_publication(struct ged *gedp,
	BObolSceneController *scene)
{
    if (!gedp || !gedp->dbip || !scene)
	FAIL("presentation publication test needs GED and Obol scene state");

    constexpr const char *path = "box.s";
    constexpr int highlighted_width = 39;
    constexpr int cleared_width = 41;
    frontier_clear(gedp);
    (void)scene->clearDatabaseSources();
    struct ged_bobol_publication_context publication;
    if (!ged_bobol_publication_begin(&publication, gedp,
	    ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_WIRE))
	FAIL("presentation publication should start a GED publication");
    const int created = ged_bobol_database_source_ensure_for_path(
	&publication, path, gedp->dbip, GED_DRAW_MODE_WIRE, 73);
    ged_bobol_publication_end(&publication);
    SoBRLDatabaseSource *source = source_for_representation(scene, path,
	SoBRLDatabaseSource::REPRESENTATION_WIRE);
    if (created <= 0 || !source || !source->realizeDatabaseWireframe() ||
	!source->hasCompactInstanceIndex() || source->getCompactInstanceCount() <= 0)
	FAIL("presentation publication should create a compact source");
    const SbString key = source->instanceKey.getValue();
    SbModernUtils::SoNodeRef source_owner(source);

    struct Observer {
	BObolSceneController &scene;
	SoBRLDatabaseSource &source;
	SbString key;
	SbBool expected;
	int editedWidth;
	uint64_t frame;
	bool called = false;
	bool partial = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<Observer *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    BObolDatabaseSourceSummary summary;
	    self.partial = !self.source.getSummary(summary) || !summary.valid ||
		summary.highlighted != self.expected ||
		self.scene.findDatabaseSourceInstance(self.key.getString()) !=
		    &self.source ||
		self.scene.getFrameRevision() != self.frame + 1;
	    for (int i = 0; i < self.source.getCompactInstanceCount(); ++i) {
		BObolCompactInstanceHandle handle;
		BObolCompactInstanceSummary compact;
		self.partial = self.partial ||
		    !self.source.getCompactInstanceHandle(i, handle) ||
		    !self.source.getCompactInstanceSummary(handle, compact) ||
		    compact.highlighted != self.expected;
	    }
	    if (self.partial)
		return;
	    BObolDatabaseSourceDisplayPatch edit;
	    edit.lineWidthValid = TRUE;
	    edit.lineWidth = self.editedWidth;
	    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(
		self.key.getString(), edit) <= 0;
	}
    };

    const auto apply = [&](enum ged_scene_reducer_operation operation,
	SbBool expected, int editedWidth) {
	Observer observer{*scene, *source, key, expected, editedWidth,
	    scene->getFrameRevision()};
	SoNodeSensor sensor(Observer::changed, &observer);
	sensor.setPriority(0);
	sensor.attach(source);
	struct ged_scene_reducer_request request = operation ==
		GED_SCENE_REDUCER_HIGHLIGHT ?
	    ged_scene_reducer_request_make_value(operation, path,
		expected ? 1.0 : 0.0) :
	    ged_scene_reducer_request_make(operation, NULL);
	request.view = ged_draw_active_view_ctx(gedp);
	request.mode = GED_DRAW_MODE_WIRE;
	const int changed = ged_draw_obol_scene_sync_transaction(gedp,
	    &request, NULL, scene);
	sensor.detach();
	return changed > 0 && observer.called && !observer.partial &&
	    !observer.failed && source->highlighted.getValue() == expected &&
	    source->lineWidth.getValue() == editedWidth;
    };

    if (!apply(GED_SCENE_REDUCER_HIGHLIGHT, TRUE, highlighted_width))
	FAIL("GED highlight must publish source and compact state together");
    if (!apply(GED_SCENE_REDUCER_HIGHLIGHTS_CLEAR, FALSE, cleared_width))
	FAIL("GED highlight clear must publish source and compact state together");

    (void)scene->removeDatabaseSourceInstance(key.getString());
    frontier_clear(gedp);
    std::puts("PASS GED presentation publication: highlight/clear callbacks see one source record");
    return 0;
}

static int
exercise_deferred_selection_acceptance(struct ged *gedp,
	BObolSceneController *scene)
{
    if (!gedp || !gedp->dbip || !scene)
	FAIL("deferred-selection test needs GED and Obol scene state");

    constexpr const char *draw_path = "progressive_root.c";
    constexpr const char *retarget_path = "nested_parent.c";
    constexpr int replacement_line_width = 43;
    struct ged_view_context *view_ctx = ged_draw_active_view_ctx(gedp);
    if (!view_ctx || !ged_view_context_display_endpoint_ensure(view_ctx))
	FAIL("deferred-selection test needs an attached display endpoint");
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	FAIL("deferred-selection test needs an Obol view controller");

    enum class CallbackEdit { Retarget, Replace, Remove };
    struct Observer {
	BObolSceneController &scene;
	struct db_i *database;
	const char *drawPath;
	const char *retargetPath;
	CallbackEdit edit;
	SbString instanceKey;
	SoBRLDatabaseSource *source = nullptr;
	uint32_t revision = 0;
	bool called = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<Observer *>(data);
	    if (self.called)
		return;
	    SoBRLDatabaseSource *current = source_for_representation(
		&self.scene, self.drawPath,
		SoBRLDatabaseSource::REPRESENTATION_WIRE);
	    if (!current ||
		frontier_selected_occurrences_for_path(current,
		    self.drawPath) <= 0)
		return;

	    BObolDatabaseSourceSummary summary;
	    self.called = true;
	    self.source = current;
	    if (!current->getSummary(summary) || !summary.valid) {
		self.failed = true;
		return;
	    }
	    self.instanceKey = summary.instanceKey;
	    self.revision = summary.sourceRevision + 1;
	    if (self.edit == CallbackEdit::Retarget) {
		self.failed = current->retargetDatabaseSourceInstance(
		    self.instanceKey.getString(), self.retargetPath,
		    self.revision) <= 0;
		return;
	    }
	    self.failed = self.scene.removeDatabaseSourceInstance(
		self.instanceKey.getString()) <= 0;
	    if (self.failed || self.edit == CallbackEdit::Remove)
		return;

	    BObolDatabaseSourcePublishState replacement;
	    replacement.sourceInstanceKey = self.instanceKey.getString();
	    replacement.sourcePath = summary.path.getString();
	    replacement.sourceRepresentationKey =
		summary.representationKey.getString();
	    replacement.database = self.database;
	    replacement.drawMode = summary.drawMode;
	    replacement.representationMode = summary.representationMode;
	    replacement.sourceRevisionValid = TRUE;
	    replacement.sourceRevision = self.revision;
	    replacement.lineWidth = replacement_line_width;
	    self.failed = self.scene.publishDatabaseSourceInstance(
		replacement) <= 0;
	}
    };

    for (CallbackEdit edit : {CallbackEdit::Retarget,
	CallbackEdit::Replace, CallbackEdit::Remove}) {
	frontier_clear(gedp);
	(void)scene->clearDatabaseSources();
	(void)ged_selection_clear(gedp, NULL);
	if (!ged_selection_batch_begin(gedp) ||
	    !ged_selection_select_path(gedp, NULL, draw_path, 1))
	    FAIL("deferred-selection test should stage semantic selection");

	Observer observer{*scene, gedp->dbip, draw_path, retarget_path, edit,
	    SbString(), nullptr, 0, false, false};
	SoNodeSensor sensor(Observer::changed, &observer);
	sensor.setPriority(0);
	sensor.attach(scene->getSceneRoot());
	struct ged_draw_appearance_settings appearance =
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
	appearance.defer_leaf_expansion = 1;
	struct ged_scene_reducer_request draw =
	    ged_scene_reducer_request_make(GED_SCENE_REDUCER_DRAW, draw_path);
	draw.view = view_ctx;
	draw.mode = GED_DRAW_MODE_WIRE;
	draw.appearance = &appearance;
	const int drawn = ged_scene_reduce(gedp, &draw, NULL);
	sensor.detach();

	BObolProgressiveOptions options;
	options.maxProviderItems = 1;
	(void)controller->advanceProgressiveWork(&options, NULL);
	BObolLodConvergenceStatus convergence;
	controller->getLodConvergenceStatus(convergence);
	SoBRLDatabaseSource *current =
	    scene->findDatabaseSourceInstance(observer.instanceKey.getString());
	BObolDatabaseSourceSummary summary;
	const bool callback_edit_survived =
	    edit == CallbackEdit::Remove ? !current :
	    current && current->getSummary(summary) && summary.valid &&
	    summary.sourceRevision == observer.revision &&
	    (edit == CallbackEdit::Retarget ?
		(current == observer.source &&
		 BU_STR_EQUAL(summary.path.getString(), retarget_path)) :
		(current != observer.source &&
		 summary.lineWidth == replacement_line_width));

	(void)ged_selection_clear(gedp, NULL);
	if (!ged_selection_batch_end(gedp))
	    FAIL("deferred-selection test should finish its selection batch");
	frontier_clear(gedp);
	(void)scene->clearDatabaseSources();

	if (drawn <= 0 || !observer.called || observer.failed ||
	    !callback_edit_survived || convergence.sourcePreparationPending) {
	    fprintf(stderr, "deferred selection edit=%d drawn=%d called=%d "
		"failed=%d survived=%d pending=%d path=%s revision=%u/%u\n",
		int(edit), drawn, observer.called, observer.failed,
		callback_edit_survived, convergence.sourcePreparationPending,
		summary.valid ? summary.path.getString() : "",
		summary.valid ? summary.sourceRevision : 0,
		observer.revision);
	    FAIL("deferred selection callback must invalidate the captured worker target");
	}
    }

    (void)ged_view_context_obol_endpoint_set(view_ctx, NULL, 0);
    std::puts("PASS GED deferred-selection acceptance: callback retarget/replacement/removal skip superseded work");
    return 0;
}

static int
exercise_source_display_publication(struct ged *gedp,
	BObolSceneController *scene)
{
    if (!gedp || !gedp->dbip || !scene)
	FAIL("source-display test needs GED and Obol scene state");

    constexpr const char *path = "box.s";
    constexpr int requested_line_style = 1;
    constexpr int requested_line_width = 7;
    constexpr int callback_line_width = 47;
    constexpr float requested_transparency = 0.375f;
    constexpr uint32_t requested_material_revision = 211;
    const unsigned char requested_color[3] = {20, 40, 60};
    const unsigned char requested_material_color[3] = {80, 100, 120};

    frontier_clear(gedp);
    (void)scene->clearDatabaseSources();
    const char *draw_box[2] = {"draw", path};
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK)
	FAIL("source-display public-ref test should draw its source");
    record_source_state shape_record = {};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	&shape_record);
    SoBRLDatabaseSource *public_source = source_for_path(scene, path);
    if (!shape_record.found || !public_source)
	FAIL("source-display public-ref test should resolve its source reference");
    SbModernUtils::SoNodeRef public_source_owner(public_source);

    struct PublicObserver {
	BObolSceneController &scene;
	SoBRLDatabaseSource &source;
	SbString key;
	uint64_t frame;
	bool called = false;
	bool partial = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<PublicObserver *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    BObolDatabaseSourceSummary summary;
	    self.partial = !self.source.getSummary(summary) || !summary.valid ||
		summary.visible ||
		self.scene.findDatabaseSourceInstance(self.key.getString()) !=
		    &self.source ||
		self.scene.getFrameRevision() != self.frame + 1;
	    if (self.partial)
		return;
	    BObolDatabaseSourceDisplayPatch edit;
	    edit.visibleValid = TRUE;
	    edit.visible = TRUE;
	    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(
		self.key.getString(), edit) <= 0;
	}
    } public_observer{*scene, *public_source,
	public_source->instanceKey.getValue(), scene->getFrameRevision()};
    SoNodeSensor public_sensor(PublicObserver::changed, &public_observer);
    public_sensor.setPriority(0);
    public_sensor.attach(public_source);
    const int public_applied = ged_draw_shape_ref_set_visible(gedp,
	shape_record.ref, 0);
    public_sensor.detach();
    if (public_applied <= 0 || !public_observer.called ||
	public_observer.partial || public_observer.failed ||
	!public_source->visible.getValue())
	FAIL("GED public shape display must retain its callback edit");

    struct ged_draw_group_record_summary group_record;
    memset(&group_record, 0, sizeof(group_record));
    if (!ged_draw_group_ref_record_summary(gedp, shape_record.group,
	    &group_record) || !group_record.path || !group_record.path[0])
	FAIL("source-display public-ref test should resolve its group reference");
    const std::string group_path = group_record.path;
    SoGroup *group_node = scene->findGroup(group_path.c_str());
    if (!group_node || !group_node->isOfType(SoBRLSceneGroup::getClassTypeId()))
	FAIL("source-display public-ref test should retain its scene group");
    auto *public_group = static_cast<SoBRLSceneGroup *>(group_node);
    SbModernUtils::SoNodeRef public_group_owner(public_group);
    struct PublicGroupObserver {
	BObolSceneController &scene;
	SoBRLSceneGroup &group;
	const std::string &path;
	uint64_t frame;
	bool called = false;
	bool partial = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<PublicGroupObserver *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    self.partial = self.group.visible.getValue() ||
		self.scene.findGroup(self.path.c_str()) != &self.group ||
		self.scene.getFrameRevision() != self.frame + 1;
	    if (self.partial)
		return;
	    self.failed = self.scene.setGroupDisplayState(self.path.c_str(),
		TRUE, self.group.selected.getValue(),
		self.group.highlighted.getValue(),
		self.group.lineStyle.getValue(), self.group.lineWidth.getValue(),
		self.group.transparency.getValue(),
		self.group.colorOverride.getValue(), self.group.color.getValue(),
		self.group.materialColorValid.getValue(),
		self.group.materialColor.getValue(),
		self.group.materialRevision.getValue()) <= 0;
	}
    } public_group_observer{*scene, *public_group, group_path,
	scene->getFrameRevision()};
    SoNodeSensor public_group_sensor(PublicGroupObserver::changed,
	&public_group_observer);
    public_group_sensor.setPriority(0);
    public_group_sensor.attach(public_group);
    const int public_group_applied = ged_draw_group_ref_set_visible(gedp,
	shape_record.group, 0);
    public_group_sensor.detach();
    if (public_group_applied <= 0 || !public_group_observer.called ||
	public_group_observer.partial || public_group_observer.failed ||
	!public_group->visible.getValue())
	FAIL("GED public group display must retain its callback edit");

    frontier_clear(gedp);
    (void)scene->clearDatabaseSources();
    struct ged_bobol_publication_context publication;
    if (!ged_bobol_publication_begin(&publication, gedp,
	    ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_WIRE))
	FAIL("source-display test should start a GED publication");
    const int created = ged_bobol_database_source_ensure_for_path(
	&publication, path, gedp->dbip, GED_DRAW_MODE_WIRE, 79);
    ged_bobol_publication_end(&publication);
    SoBRLDatabaseSource *source = source_for_representation(scene, path,
	SoBRLDatabaseSource::REPRESENTATION_WIRE);
    if (created <= 0 || !source)
	FAIL("source-display test should create its source");
    const SbString key = source->instanceKey.getValue();
    SbModernUtils::SoNodeRef source_owner(source);

    struct Observer {
	BObolSceneController &scene;
	SoBRLDatabaseSource &source;
	SbString key;
	uint64_t frame;
	bool called = false;
	bool partial = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<Observer *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    BObolDatabaseSourceSummary summary;
	    self.partial = !self.source.getSummary(summary) || !summary.valid ||
		summary.drawMode != SoBRLDatabaseSource::SHADED ||
		summary.representationMode !=
		    SoBRLDatabaseSource::REPRESENTATION_SHADED ||
		summary.visible || !summary.selected || !summary.highlighted ||
		summary.lineStyle != requested_line_style ||
		summary.lineWidth != requested_line_width ||
		fabsf(summary.transparency - requested_transparency) > 1.0e-6f ||
		!summary.colorOverride ||
		fabsf(summary.color[0] - 20.0f / 255.0f) > 1.0e-6f ||
		fabsf(summary.color[1] - 40.0f / 255.0f) > 1.0e-6f ||
		fabsf(summary.color[2] - 60.0f / 255.0f) > 1.0e-6f ||
		!summary.materialColorValid ||
		fabsf(summary.materialColor[0] - 80.0f / 255.0f) > 1.0e-6f ||
		fabsf(summary.materialColor[1] - 100.0f / 255.0f) > 1.0e-6f ||
		fabsf(summary.materialColor[2] - 120.0f / 255.0f) > 1.0e-6f ||
		summary.materialRevision != requested_material_revision ||
		self.scene.findDatabaseSourceInstance(self.key.getString()) !=
		    &self.source ||
		self.scene.getFrameRevision() != self.frame + 1;
	    if (self.partial)
		return;
	    BObolDatabaseSourceDisplayPatch edit;
	    edit.lineWidthValid = TRUE;
	    edit.lineWidth = callback_line_width;
	    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(
		self.key.getString(), edit) <= 0;
	}
    } observer{*scene, *source, key, scene->getFrameRevision()};

    SoNodeSensor sensor(Observer::changed, &observer);
    sensor.setPriority(0);
    sensor.attach(source);
    const int applied = ged_draw_obol_database_source_update_display_for_path(
	gedp, path,
	1, 0,
	1, 1,
	1, 1,
	1, GED_DRAW_MODE_SHADED,
	1, requested_line_style,
	1, requested_line_width,
	1, requested_transparency,
	1, requested_color,
	1, requested_material_color,
	1, requested_material_revision);
    sensor.detach();

    BObolDatabaseSourceSummary final_summary;
    if (applied <= 0 || !observer.called || observer.partial ||
	observer.failed || !source->getSummary(final_summary) ||
	!final_summary.valid ||
	final_summary.lineWidth != callback_line_width) {
	fprintf(stderr, "source display applied=%d called=%d partial=%d "
	    "failed=%d width=%d mode=%d representation=%d frame=%" PRIu64
	    "/%" PRIu64 "\n", applied, observer.called, observer.partial,
	    observer.failed, final_summary.valid ? final_summary.lineWidth : -1,
	    final_summary.valid ? final_summary.drawMode : -1,
	    final_summary.valid ? final_summary.representationMode : -1,
	    scene->getFrameRevision(), observer.frame + 1);
	FAIL("GED source display must publish one complete record and retain callback edits");
    }

    (void)scene->removeDatabaseSourceInstance(key.getString());

    enum class PendingEdit { Change, Replace, Remove };
    constexpr int requested_pending_width = 9;
    constexpr int changed_pending_width = 53;
    constexpr int replacement_pending_width = 59;
    for (PendingEdit edit : {PendingEdit::Change, PendingEdit::Replace,
	PendingEdit::Remove}) {
	(void)scene->clearDatabaseSources();
	struct ged_bobol_publication_context wire_publication;
	if (!ged_bobol_publication_begin(&wire_publication, gedp,
		ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_WIRE) ||
	    ged_bobol_database_source_ensure_for_path(&wire_publication, path,
		gedp->dbip, GED_DRAW_MODE_WIRE, 83) <= 0) {
	    ged_bobol_publication_end(&wire_publication);
	    FAIL("source-display multi-target test should create its wire source");
	}
	ged_bobol_publication_end(&wire_publication);
	struct ged_bobol_publication_context shaded_publication;
	if (!ged_bobol_publication_begin(&shaded_publication, gedp,
		ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_SHADED) ||
	    ged_bobol_database_source_ensure_for_path(&shaded_publication, path,
		gedp->dbip, GED_DRAW_MODE_SHADED, 83) <= 0) {
	    ged_bobol_publication_end(&shaded_publication);
	    FAIL("source-display multi-target test should create its shaded source");
	}
	ged_bobol_publication_end(&shaded_publication);

	SoBRLDatabaseSource *wire = source_for_representation(scene, path,
	    SoBRLDatabaseSource::REPRESENTATION_WIRE);
	SoBRLDatabaseSource *shaded = source_for_representation(scene, path,
	    SoBRLDatabaseSource::REPRESENTATION_SHADED);
	if (!wire || !shaded || wire == shaded)
	    FAIL("source-display multi-target test needs distinct source instances");
	SbModernUtils::SoNodeRef wire_owner(wire);
	SbModernUtils::SoNodeRef shaded_owner(shaded);

	struct PendingObserver {
	    BObolSceneController &scene;
	    SoBRLDatabaseSource &wire;
	    SoBRLDatabaseSource &shaded;
	    struct db_i *database;
	    PendingEdit edit;
	    SbString publishedKey;
	    SbString pendingKey;
	    uint32_t replacementRevision = 0;
	    bool called = false;
	    bool failed = false;
	    static void changed(void *data, SoSensor *sensor)
	    {
		auto &self = *static_cast<PendingObserver *>(data);
		if (self.called)
		    return;
		self.called = true;
		auto *node_sensor = static_cast<SoNodeSensor *>(sensor);
		auto *published = static_cast<SoBRLDatabaseSource *>(
		    node_sensor->getAttachedNode());
		SoBRLDatabaseSource *pending = published == &self.wire ?
		    &self.shaded : &self.wire;
		self.publishedKey = published->instanceKey.getValue();
		self.pendingKey = pending->instanceKey.getValue();

		if (self.edit == PendingEdit::Change) {
		    BObolDatabaseSourceDisplayPatch patch;
		    patch.lineWidthValid = TRUE;
		    patch.lineWidth = changed_pending_width;
		    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(
			self.pendingKey.getString(), patch) <= 0;
		    return;
		}
		if (self.edit == PendingEdit::Remove) {
		    self.failed = self.scene.removeDatabaseSourceInstance(
			self.pendingKey.getString()) <= 0;
		    return;
		}

		BObolDatabaseSourceSummary summary;
		if (!pending->getSummary(summary) || !summary.valid) {
		    self.failed = true;
		    return;
		}
		BObolDatabaseSourcePublishState replacement;
		replacement.sourceInstanceKey = self.pendingKey.getString();
		replacement.sourcePath = summary.path.getString();
		replacement.sourceRepresentationKey =
		    summary.representationKey.getString();
		replacement.database = self.database;
		replacement.drawMode = summary.drawMode;
		replacement.representationMode = summary.representationMode;
		replacement.sourceRevisionValid = TRUE;
		self.replacementRevision = summary.sourceRevision + 1;
		replacement.sourceRevision = self.replacementRevision;
		replacement.lineWidth = replacement_pending_width;
		replacement.materialPolicyValid = TRUE;
		replacement.materialPolicy = summary.materialPolicy;
		self.failed = self.scene.removeDatabaseSourceInstance(
		    self.pendingKey.getString()) <= 0 ||
		    self.scene.publishDatabaseSourceInstance(replacement) <= 0;
	    }
	} pending_observer{*scene, *wire, *shaded, gedp->dbip, edit,
	    SbString(), SbString()};

	SoNodeSensor wire_sensor(PendingObserver::changed, &pending_observer);
	SoNodeSensor shaded_sensor(PendingObserver::changed, &pending_observer);
	wire_sensor.setPriority(0);
	shaded_sensor.setPriority(0);
	wire_sensor.attach(wire);
	shaded_sensor.attach(shaded);
	const int multi_applied =
	    ged_draw_obol_database_source_update_display_for_path(
		gedp, path,
		0, 0, 0, 0, 0, 0, 0, GED_DRAW_MODE_WIRE,
		0, 0, 1, requested_pending_width,
		0, 0.0, 0, NULL, 0, NULL, 0, 0);
	wire_sensor.detach();
	shaded_sensor.detach();

	SoBRLDatabaseSource *published = scene->findDatabaseSourceInstance(
	    pending_observer.publishedKey.getString());
	SoBRLDatabaseSource *pending = scene->findDatabaseSourceInstance(
	    pending_observer.pendingKey.getString());
	bool current_target_valid = multi_applied > 0 && pending_observer.called &&
	    !pending_observer.failed && published &&
	    published->lineWidth.getValue() == requested_pending_width;
	if (edit == PendingEdit::Change) {
	    current_target_valid = current_target_valid && pending &&
		pending->lineWidth.getValue() == changed_pending_width;
	} else if (edit == PendingEdit::Replace) {
	    current_target_valid = current_target_valid && pending &&
		pending != wire && pending != shaded &&
		pending->lineWidth.getValue() == replacement_pending_width &&
		pending->sourceRevision.getValue() ==
		    pending_observer.replacementRevision;
	} else {
	    current_target_valid = current_target_valid && !pending;
	}
	if (!current_target_valid) {
	    fprintf(stderr, "source display pending edit=%d applied=%d called=%d "
		"failed=%d published=%p pending=%p widths=%d/%d revision=%u/%u\n",
		int(edit), multi_applied, pending_observer.called,
		pending_observer.failed,
		(void *)published, (void *)pending,
		published ? published->lineWidth.getValue() : -1,
		pending ? pending->lineWidth.getValue() : -1,
		pending ? pending->sourceRevision.getValue() : 0,
		pending_observer.replacementRevision);
	    FAIL("GED source display must preserve callback changes to pending targets");
	}
    }

    (void)scene->clearDatabaseSources();
    frontier_clear(gedp);
    std::puts("PASS GED source-display publication: complete records and current callback targets");
    return 0;
}

static int
exercise_source_record_publication(struct ged *gedp,
	BObolSceneController *scene)
{
    if (!gedp || !gedp->dbip || !scene)
	FAIL("source-record test needs GED and Obol scene state");

    enum class CallbackEdit { Change, Replace, Remove };
    constexpr const char *path = "box.s";
    constexpr int changed_line_width = 35;
    constexpr int replacement_line_width = 37;
    for (CallbackEdit edit : {CallbackEdit::Change,
	CallbackEdit::Replace, CallbackEdit::Remove}) {
	frontier_clear(gedp);
	(void)scene->clearDatabaseSources();
	struct ged_bobol_publication_context publication;
	if (!ged_bobol_publication_begin(&publication, gedp,
		ged_draw_active_view_ctx(gedp), GED_DRAW_MODE_WIRE))
	    FAIL("source-record test should start a GED publication");
	const int created = ged_bobol_database_source_ensure_for_path(
	    &publication, path, gedp->dbip, GED_DRAW_MODE_WIRE, 71);
	ged_bobol_publication_end(&publication);
	SoBRLDatabaseSource *source = source_for_representation(scene, path,
	    SoBRLDatabaseSource::REPRESENTATION_WIRE);
	if (created <= 0 || !source)
	    FAIL("source-record test should create its source");

	SbVec3f points[2] = {
	    SbVec3f(0.0f, 0.0f, 0.0f),
	    SbVec3f(1.0f, 1.0f, 1.0f)
	};
	int32_t commands[2] = {
	    SoBRLVListShape::MOVE,
	    SoBRLVListShape::DRAW
	};
	BObolExternalLineSet proxy;
	proxy.points = points;
	proxy.commands = commands;
	proxy.count = 2;
	proxy.sourceType = "proxy";
	proxy.geometryKind = "aabb";
	const SbString instance_key = source->instanceKey.getValue();
	if (scene->publishDatabaseSourceInstanceExternalLineSet(
		instance_key.getString(), proxy) <= 0)
	    FAIL("source-record test should publish its retained geometry");

	struct ged_draw_obol_database_source_record record;
	if (!ged_draw_obol_database_source_record_for_path(gedp, path,
		&record))
	    FAIL("source-record test should read its initial record");
	record.draw_mode = GED_DRAW_MODE_HIDDEN_LINE;
	record.source_revision += 17;
	record.inputs_revision += 23;
	record.realization_status =
	    GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_CURRENT;
	record.realized_source_revision = record.source_revision;
	record.realized_inputs_revision = record.inputs_revision;
	record.stale_reason = GED_DRAW_STALE_NONE;
	record.material_policy =
	    GED_DRAW_OBOL_DATABASE_SOURCE_MATERIAL_INHERIT;
	record.realization_role_flags =
	    SoBRLDatabaseSource::REALIZATION_ROLE_MESH;
	record.realization_view_dependent = 1;
	record.realization_csg_lod_enabled = 0;
	record.realization_mesh_lod_enabled = 1;
	record.realization_view_scale = 13.0;
	record.realization_lod_scale = 2.5;
	record.realization_bot_threshold = 91;
	record.realization_curve_scale = 4.0;
	record.realization_point_scale = 5.0;
	SbModernUtils::SoNodeRef source_owner(source);

	struct Observer {
	    BObolSceneController &scene;
	    SoBRLDatabaseSource &source;
	    const struct ged_draw_obol_database_source_record &record;
	    struct db_i *database;
	    CallbackEdit edit;
	    uint64_t frame;
	    bool called = false;
	    bool partial = false;
	    bool failed = false;
	    static void changed(void *data, SoSensor *)
	    {
		auto &self = *static_cast<Observer *>(data);
		if (self.called)
		    return;
		self.called = true;
		BObolDatabaseSourceSummary current;
		self.partial = !self.source.getSummary(current) ||
		    !current.valid || current.drawMode !=
			SoBRLDatabaseSource::WIREFRAME ||
		    current.representationMode !=
			SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE ||
		    current.sourceRevision !=
			(uint32_t)self.record.source_revision ||
		    current.inputsRevision !=
			(uint32_t)self.record.inputs_revision ||
		    current.materialPolicy !=
			SoBRLDatabaseSource::MATERIAL_INHERIT ||
		    current.realizationStatus !=
			SoBRLDatabaseSource::REALIZED || current.stale ||
		    current.staleReason != SoBRLDatabaseSource::STALE_NONE ||
		    current.realizedSourceRevision != current.sourceRevision ||
		    current.realizedInputsRevision != current.inputsRevision ||
		    current.realizationRoleFlags !=
			SoBRLDatabaseSource::REALIZATION_ROLE_MESH ||
		    !current.realizationViewDependent ||
		    current.realizationCsgLodEnabled ||
		    !current.realizationMeshLodEnabled ||
		    std::abs(current.realizationViewScale - 13.0f) > 0.001f ||
		    std::abs(current.realizationLodScale - 2.5f) > 0.001f ||
		    current.realizationBotThreshold != 91 ||
		    std::abs(current.realizationCurveScale - 4.0f) > 0.001f ||
		    std::abs(current.realizationPointScale - 5.0f) > 0.001f ||
		    self.scene.findDatabaseSourceInstance(
			current.instanceKey.getString()) != &self.source ||
		    self.scene.getFrameRevision() != self.frame + 1;
		if (self.partial)
		    return;

		const char *key = current.instanceKey.getString();
		if (self.edit == CallbackEdit::Change) {
		    BObolDatabaseSourceDisplayPatch patch;
		    patch.lineWidthValid = TRUE;
		    patch.lineWidth = changed_line_width;
		    self.failed =
			self.scene.setDatabaseSourceInstanceDisplayPatch(
			    key, patch) <= 0;
		    return;
		}
		self.failed = self.scene.removeDatabaseSourceInstance(key) <= 0;
		if (self.edit == CallbackEdit::Remove)
		    return;
		BObolDatabaseSourcePublishState replacement;
		replacement.sourceInstanceKey = key;
		replacement.sourcePath = current.path.getString();
		replacement.sourceRepresentationKey =
		    current.representationKey.getString();
		replacement.database = self.database;
		replacement.drawMode = current.drawMode;
		replacement.representationMode = current.representationMode;
		replacement.sourceRevisionValid = TRUE;
		replacement.sourceRevision = current.sourceRevision + 1;
		replacement.lineWidth = replacement_line_width;
		replacement.roleFlagsValid = TRUE;
		replacement.roleFlags =
		    SoBRLDatabaseSource::REALIZATION_ROLE_NONE;
		self.failed =
		    self.scene.publishDatabaseSourceInstance(replacement) < 0 ||
		    self.failed;
	    }
	} observer{*scene, *source, record, gedp->dbip, edit,
	    scene->getFrameRevision()};

	SoNodeSensor sensor(Observer::changed, &observer);
	sensor.setPriority(0);
	sensor.attach(source);
	const int applied = ged_draw_obol_database_source_apply_record_for_path(
	    gedp, path, &record);
	sensor.detach();

	SoBRLDatabaseSource *current =
	    scene->findDatabaseSourceInstance(instance_key.getString());
	bool final_state_valid = applied > 0 && observer.called &&
	    !observer.partial && !observer.failed;
	if (edit == CallbackEdit::Change) {
	    final_state_valid = final_state_valid && current == source &&
		source->lineWidth.getValue() == changed_line_width &&
		source->realizationStatus.getValue() ==
		    SoBRLDatabaseSource::REALIZED && !source->stale.getValue();
	} else if (edit == CallbackEdit::Replace) {
	    final_state_valid = final_state_valid && current &&
		current != source &&
		current->lineWidth.getValue() == replacement_line_width &&
		current->realizationStatus.getValue() ==
		    SoBRLDatabaseSource::UNREALIZED && current->stale.getValue() &&
		current->realizationRoleFlags.getValue() ==
		    SoBRLDatabaseSource::REALIZATION_ROLE_NONE;
	} else {
	    final_state_valid = final_state_valid && !current;
	}
	if (!final_state_valid) {
	    fprintf(stderr, "source record edit=%d applied=%d called=%d "
		"partial=%d failed=%d current=%p original=%p width=%d "
		"status=%d stale=%d roles=%d frame=%" PRIu64 "\n",
		int(edit), applied, observer.called, observer.partial,
		observer.failed, (void *)current, (void *)source,
		current ? current->lineWidth.getValue() : -1,
		current ? current->realizationStatus.getValue() : -1,
		current ? current->stale.getValue() : -1,
		current ? current->realizationRoleFlags.getValue() : -1,
		scene->getFrameRevision());
	    FAIL("source-record publication must commit one complete record before observers");
	}
	(void)scene->removeDatabaseSourceInstance(instance_key.getString());
	std::printf("PASS GED source-record publication edit%d: complete state and current callback target\n",
	    int(edit));
    }

    frontier_clear(gedp);
    return 0;
}

static int
exercise_database_rename_publication(struct ged *gedp,
	BObolSceneController *scene)
{
    constexpr const char *old_path = "rename_source.s";
    constexpr const char *new_path = "rename_callback.s";
    constexpr const char *new_key = "/rename_callback.s";
    constexpr int initial_line_width = 11;
    constexpr int callback_line_width = 29;

    const char *draw_source[2] = {"draw", old_path};
    if (ged_exec_draw(gedp, 2, draw_source) != BRLCAD_OK)
	FAIL("database rename publication setup draw should succeed");
    SoBRLDatabaseSource *source = source_for_path(scene, old_path);
    if (!source || scene->setDatabaseSourceState(old_path, FALSE, 0, 7,
	    TRUE, FALSE, FALSE, 0, initial_line_width, 0.0f, FALSE,
	    SbColor(1.0f, 1.0f, 1.0f), FALSE,
	    SbColor(1.0f, 1.0f, 1.0f), 0) <= 0)
	FAIL("database rename publication setup state should succeed");

    struct Observer {
	BObolSceneController &scene;
	SoBRLDatabaseSource *source;
	uint64_t frame;
	bool called = false;
	bool partial = false;
	bool failed = false;
	static void changed(void *data, SoSensor *)
	{
	    auto &self = *static_cast<Observer *>(data);
	    if (self.called)
		return;
	    self.called = true;
	    BObolDatabaseSourceSummary summary;
	    BObolRealizedShapeSummary shape;
	    SoBRLDatabaseSource *current =
		self.scene.findDatabaseSourceInstance(new_key);
	    self.partial = current != self.source ||
		self.scene.findDatabaseSourceInstance("/rename_source.s") ||
		self.scene.findGroup(old_path) ||
		!self.scene.findGroup(new_path) ||
		!current->getSummary(summary) || !summary.valid ||
		!path_equal(summary.path.getString(), new_path) ||
		!BU_STR_EQUAL(summary.instanceKey.getString(), new_key) ||
		!BU_STR_EQUAL(summary.representationKey.getString(), new_key) ||
		summary.lineWidth != initial_line_width ||
		!current->getRealizedShapeSummary(0, shape) ||
		!path_equal(shape.path.getString(), new_path) ||
		!path_equal(shape.ownerSourcePath.getString(), new_path) ||
		!BU_STR_EQUAL(shape.ownerSourceInstanceKey.getString(), new_key) ||
		self.scene.getFrameRevision() <= self.frame;
	    if (self.partial)
		return;
	    BObolDatabaseSourceDisplayPatch patch;
	    patch.lineWidthValid = TRUE;
	    patch.lineWidth = callback_line_width;
	    self.failed = self.scene.setDatabaseSourceInstanceDisplayPatch(
		new_key, patch) <= 0;
	}
    } observer{*scene, source, scene->getFrameRevision()};

    SoNodeSensor sensor(Observer::changed, &observer);
    sensor.setPriority(0);
    sensor.attach(scene->getSceneRoot());
    const char *rename_source[3] = {"move", old_path, new_path};
    const int renamed = ged_exec(gedp, 3, rename_source);
    sensor.detach();
    SoBRLDatabaseSource *renamed_source = source_for_path(scene, new_path);
    BObolDatabaseSourceSummary renamed_summary;
    if (renamed != BRLCAD_OK || !observer.called || observer.partial ||
	observer.failed || renamed_source != source ||
	!renamed_source->getSummary(renamed_summary) ||
	renamed_summary.lineWidth != callback_line_width)
	FAIL("database rename must publish complete state and retain callback edits");

    const char *restore_source[3] = {"move", new_path, old_path};
    if (ged_exec(gedp, 3, restore_source) != BRLCAD_OK ||
	!source_for_path(scene, old_path) || source_for_path(scene, new_path))
	FAIL("database rename publication restore should succeed");
    const char *erase_source[2] = {"erase", old_path};
    if (ged_exec_erase(gedp, 2, erase_source) != BRLCAD_OK)
	FAIL("database rename publication cleanup should succeed");

    constexpr const char *old_component = "nested_child.c";
    constexpr const char *new_component = "nested_child_callback.c";
    constexpr const char *old_nested =
	"nested_parent.c/nested_child.c/nested_leaf.s";
    constexpr const char *new_nested =
	"nested_parent.c/nested_child_callback.c/nested_leaf.s";
    const char *draw_nested[2] = {"draw", old_nested};
    if (ged_exec_draw(gedp, 2, draw_nested) != BRLCAD_OK ||
	!source_for_path(scene, old_nested))
	FAIL("nested database rename setup draw should succeed");
    const char *rename_nested[3] = {"mvall", old_component, new_component};
    if (ged_exec(gedp, 3, rename_nested) != BRLCAD_OK)
	FAIL("nested database component rename should succeed");
    SoBRLDatabaseSource *nested = source_for_path(scene, new_nested);
    BObolDatabaseSourceSummary nested_summary;
    if (source_for_path(scene, old_nested) ||
	!nested || !nested->getSummary(nested_summary) ||
	!path_equal(nested_summary.path.getString(), new_nested) ||
	!path_equal(nested_summary.instanceKey.getString(), new_nested) ||
	!path_equal(nested_summary.representationKey.getString(), new_nested) ||
	(nested->hasCompactInstanceIndex() &&
	 (nested->getCompactInstanceCountForPath(old_nested, TRUE) != 0 ||
	  nested->getCompactInstanceCountForPath(new_nested, TRUE) <= 0)))
	FAIL("nested database rename must retarget the component without collapsing the source path");

    const char *restore_nested[3] = {"mvall", new_component, old_component};
    if (ged_exec(gedp, 3, restore_nested) != BRLCAD_OK ||
	!source_for_path(scene, old_nested) || source_for_path(scene, new_nested))
	FAIL("nested database component rename restore should succeed");
    const char *erase_nested[2] = {"erase", old_nested};
    if (ged_exec_erase(gedp, 2, erase_nested) != BRLCAD_OK)
	FAIL("nested database rename cleanup should succeed");

    std::puts("PASS GED database rename publication: complete callback state, retained edit and nested component identity");
    return 0;
}

struct faceplate_composition_observer {
    faceplate_composition_observer(BObolViewController *view_controller,
	struct ged_view_context *view_context, bool rectangle, bool scale) :
	controller(view_controller), view_ctx(view_context),
	presentation_revision(view_controller ?
	    view_controller->features().presentationRevision() : 0),
	expect_rectangle(rectangle), expect_scale(scale)
    {
    }

    BObolViewController *controller;
    struct ged_view_context *view_ctx;
    uint64_t presentation_revision;
    size_t graph_callbacks = 0;
    size_t frame_callbacks = 0;
    bool partial = false;
    bool expect_rectangle;
    bool expect_scale;

    bool completeFeature(const char *name) const
    {
	BObolFeatureRecord record;
	const BObolFeatureHandle handle = controller->features().find(name);
	return handle.isValid() && controller->features().record(handle, record) &&
	    record.overlay.isOverlay && record.overlay.ownerToken == view_ctx &&
	    record.overlay.overlayClass == BObolOverlayClass::Faceplate &&
	    record.overlay.role == BObolOverlayRole::Screen &&
	    record.overlay.order == BObolOverlayOrder::Screen;
    }

    bool complete() const
    {
	if (!controller)
	    return false;
	const BObolFeatureHandle rectangle_handle =
	    controller->features().find("_faceplate/interactive_rect");
	const BObolFeatureHandle scale_handle =
	    controller->features().find("_faceplate/scale");
	const BObolFeatureHandle scale_labels_handle =
	    controller->features().find("_faceplate/scale_labels");
	return completeFeature("_faceplate/center_dot") &&
	    rectangle_handle.isValid() == expect_rectangle &&
	    (!expect_rectangle ||
	     completeFeature("_faceplate/interactive_rect")) &&
	    scale_handle.isValid() == expect_scale &&
	    scale_labels_handle.isValid() == expect_scale &&
	    (!expect_scale || (completeFeature("_faceplate/scale") &&
	     completeFeature("_faceplate/scale_labels"))) &&
	    controller->features().presentationRevision() ==
		presentation_revision + 1 && controller->isRenderRequested();
    }

    static void graph_changed(void *data, SoSensor *)
    {
	auto &observer = *static_cast<faceplate_composition_observer *>(data);
	++observer.graph_callbacks;
	if (!observer.complete()) {
	    fprintf(stderr,
		"faceplate composition partial graph callback: revision=%" PRIu64 " center=%d rectangle=%d\n",
		observer.controller->features().presentationRevision(),
		observer.controller->features().exists(
		    "_faceplate/center_dot") ? 1 : 0,
		observer.controller->features().exists(
		    "_faceplate/interactive_rect") ? 1 : 0);
	    observer.partial = true;
	}
    }

    static void frame_requested(void *data, const char *reason)
    {
	auto &observer = *static_cast<faceplate_composition_observer *>(data);
	if (!reason || !BU_STR_EQUAL(reason, "view-feature-store"))
	    return;
	++observer.frame_callbacks;
	if (!observer.complete()) {
	    fprintf(stderr,
		"faceplate composition partial frame callback: reason=%s revision=%" PRIu64 " center=%d rectangle=%d\n",
		reason ? reason : "", observer.controller->features().presentationRevision(),
		observer.controller->features().exists(
		    "_faceplate/center_dot") ? 1 : 0,
		observer.controller->features().exists(
		    "_faceplate/interactive_rect") ? 1 : 0);
	    observer.partial = true;
	}
    }
};

struct edit_preview_publication_observer {
    edit_preview_publication_observer(BObolViewController *view_controller,
	struct ged_view_context *view_context, const char *feature_name,
	const char *identity_value, const char *intent_value,
	const char *role_value, const point_t *ged_points,
	const int *ged_commands, size_t count, uint32_t source_revision,
	uint32_t inputs_revision, const unsigned char color[3],
	bool require_nonempty_geometry = false) :
	controller(view_controller), view_ctx(view_context), name(feature_name),
	identity(identity_value), intent(intent_value), role(role_value),
	sourceRevision(source_revision), inputsRevision(inputs_revision),
	expectColor(color != NULL),
	requireNonemptyGeometry(require_nonempty_geometry),
	presentationRevision(view_controller ?
	    view_controller->features().presentationRevision() : 0),
	root(changed, this)
    {
	for (size_t i = 0; i < count; ++i) {
	    points.emplace_back(static_cast<float>(ged_points[i][X]),
		static_cast<float>(ged_points[i][Y]),
		static_cast<float>(ged_points[i][Z]));
	    const int command = ged_commands ? ged_commands[i] : -1;
	    commands.push_back(command == GED_DRAW_VIEW_LINE_POINT_DRAW ?
		static_cast<int32_t>(BObolLineCommand::Point) :
		command == GED_DRAW_VIEW_LINE_MOVE ?
		static_cast<int32_t>(BObolLineCommand::Move) :
		static_cast<int32_t>(i ? BObolLineCommand::Draw :
		    BObolLineCommand::Move));
	}
	if (color) {
	    expectedColor = SbColor(static_cast<float>(color[0]) / 255.0f,
		static_cast<float>(color[1]) / 255.0f,
		static_cast<float>(color[2]) / 255.0f);
	}
	root.setPriority(0);
	root.attach(controller->getSceneRoot());
	controller->setFrameRequestCallback(frame_requested, this);
    }

    ~edit_preview_publication_observer()
    {
	controller->clearFrameRequestCallback(this);
	root.detach();
    }

    bool complete() const
    {
	BObolFeatureRecord record;
	const BObolFeatureHandle handle = controller->features().find(name,
	    BOBOL_FEATURE_SCOPE_LOCAL);
	if (!handle.isValid() || !controller->features().record(handle, record))
	    return false;
	const bool geometry_complete = requireNonemptyGeometry ?
	    !record.points.empty() &&
	    record.commands.size() == record.points.size() :
	    record.points == points && record.commands == commands;
	return record.kind == BObolFeatureKind::EditPreview &&
	    record.scope == BObolFeatureScope::Local &&
	    record.identity == identity && record.editIntentId == intent &&
	    record.editIntentRole == role && geometry_complete &&
	    record.sourceRevision == sourceRevision &&
	    record.inputsRevision == inputsRevision &&
	    record.style.hasColor == expectColor &&
	    (!expectColor || record.style.color == expectedColor) &&
	    record.overlay.isOverlay && record.overlay.ownerToken == view_ctx &&
	    record.overlay.overlayClass == BObolOverlayClass::EditHandle &&
	    record.overlay.lifecycle == BObolOverlayLifecycle::PerTool &&
	    record.overlay.order == BObolOverlayOrder::PostTransparent &&
	    record.overlay.sourcePath == identity &&
	    controller->features().presentationRevision() ==
		presentationRevision + 1 && controller->isRenderRequested();
    }

    static void changed(void *data, SoSensor *)
    {
	auto &observer =
	    *static_cast<edit_preview_publication_observer *>(data);
	++observer.graphCallbacks;
	observer.partial = observer.partial || !observer.complete();
    }

    static void frame_requested(void *data, const char *reason)
    {
	if (!reason || !BU_STR_EQUAL(reason, "view-feature-store"))
	    return;
	auto &observer =
	    *static_cast<edit_preview_publication_observer *>(data);
	++observer.frameCallbacks;
	observer.partial = observer.partial || !observer.complete();
    }

    BObolViewController *controller;
    struct ged_view_context *view_ctx;
    const char *name;
    const char *identity;
    const char *intent;
    const char *role;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    uint32_t sourceRevision;
    uint32_t inputsRevision;
    bool expectColor;
    bool requireNonemptyGeometry;
    SbColor expectedColor;
    uint64_t presentationRevision;
    SoNodeSensor root;
    size_t graphCallbacks = 0;
    size_t frameCallbacks = 0;
    bool partial = false;
};

static int
exercise_edit_transaction_publication(struct ged *gedp,
	struct ged_view_context *view_ctx, BObolViewController *controller)
{
    if (!gedp || !gedp->dbip || !view_ctx || !controller)
	FAIL("edit transaction publication requires an attached view controller");

    static const char ensure_name[] = "edit-transaction::ensure";
    static const char ensure_source[] = "edit-transaction::ensure-source.s";
    controller->clearRenderRequest();
    ged_view_edit_ref ensured = GED_VIEW_EDIT_REF_NULL;
    {
	edit_preview_publication_observer observer(controller, view_ctx,
	    ensure_name, ensure_source, ensure_name, "edit-handle", NULL,
	    NULL, 0, 1, 1, NULL);
	ensured = ged_view_edit_overlay_ensure(view_ctx, ensure_name,
	    ensure_source);
	if (ged_view_edit_ref_is_null(ensured) || !observer.complete() ||
	    observer.partial || observer.graphCallbacks != 1 ||
	    observer.frameCallbacks != 1)
	    FAIL("edit overlay ensure exposed a partial retained publication");
    }
    if (!ged_view_edit_remove_ref(view_ctx, ensured))
	FAIL("edit overlay ensure fixture should remove its retained preview");

    static const char feature_name[] = "edit-transaction::preview";
    static const char initial_source[] = "edit-transaction::initial.s";
    static const char initial_intent[] = "edit-transaction::translate";
    static const char initial_role[] = "translate-handle";
    static const uint32_t initial_source_revision = 41;
    static const uint32_t initial_inputs_revision = 42;
    point_t initial_points[3] = {
	{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 2.0, 0.0}
    };
    int initial_commands[3] = {
	GED_DRAW_VIEW_LINE_MOVE, GED_DRAW_VIEW_LINE_DRAW,
	GED_DRAW_VIEW_LINE_DRAW
    };
    const unsigned char initial_color[3] = {23, 67, 149};
    struct ged_view_edit_transaction transaction =
	ged_view_edit_transaction_default();
    transaction.event = GED_VIEW_EDIT_PREVIEW_BEGIN;
    transaction.feature_name = feature_name;
    transaction.source_path = initial_source;
    transaction.edit_intent_id = initial_intent;
    transaction.edit_intent_role = initial_role;
    transaction.points = initial_points;
    transaction.commands = initial_commands;
    transaction.point_count = 3;
    transaction.source_revision = initial_source_revision;
    transaction.inputs_revision = initial_inputs_revision;
    transaction.color_valid = 1;
    VMOVE(transaction.color, initial_color);

    controller->clearRenderRequest();
    ged_view_edit_ref feature = GED_VIEW_EDIT_REF_NULL;
    {
	edit_preview_publication_observer observer(controller, view_ctx,
	    feature_name, initial_source, initial_intent, initial_role,
	    initial_points, initial_commands, 3, initial_source_revision,
	    initial_inputs_revision, initial_color);
	if (!ged_view_edit_transaction_apply(view_ctx, &transaction, &feature) ||
	    ged_view_edit_ref_is_null(feature) || !observer.complete() ||
	    observer.partial || observer.graphCallbacks != 1 ||
	    observer.frameCallbacks != 1)
	    FAIL("initial edit transaction exposed a partial retained publication");
    }

    static const char replacement_source[] =
	"edit-transaction::replacement.s";
    static const char replacement_intent[] = "edit-transaction::rotate";
    static const char replacement_role[] = "rotate-handle";
    static const uint32_t replacement_source_revision = 51;
    static const uint32_t replacement_inputs_revision = 52;
    point_t replacement_points[4] = {
	{0.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
	{3.0, 3.0, 0.0}, {0.0, 3.0, 0.0}
    };
    int replacement_commands[4] = {
	GED_DRAW_VIEW_LINE_MOVE, GED_DRAW_VIEW_LINE_DRAW,
	GED_DRAW_VIEW_LINE_DRAW, GED_DRAW_VIEW_LINE_DRAW
    };
    const unsigned char replacement_color[3] = {181, 83, 29};
    transaction.event = GED_VIEW_EDIT_PREVIEW_UPDATE;
    transaction.feature = feature;
    transaction.source_path = replacement_source;
    transaction.edit_intent_id = replacement_intent;
    transaction.edit_intent_role = replacement_role;
    transaction.points = replacement_points;
    transaction.commands = replacement_commands;
    transaction.point_count = 4;
    transaction.source_revision = replacement_source_revision;
    transaction.inputs_revision = replacement_inputs_revision;
    VMOVE(transaction.color, replacement_color);

    controller->clearRenderRequest();
    {
	edit_preview_publication_observer observer(controller, view_ctx,
	    feature_name, replacement_source, replacement_intent,
	    replacement_role, replacement_points, replacement_commands, 4,
	    replacement_source_revision, replacement_inputs_revision,
	    replacement_color);
	ged_view_edit_ref replacement = GED_VIEW_EDIT_REF_NULL;
	if (!ged_view_edit_transaction_apply(view_ctx, &transaction,
		&replacement) || ged_view_edit_ref_is_null(replacement) ||
	    replacement.id != feature.id || !observer.complete() ||
	    observer.partial || observer.graphCallbacks != 1 ||
	    observer.frameCallbacks != 1)
	    FAIL("replacement edit transaction exposed a partial retained publication");
	feature = replacement;
    }

    const uint64_t no_op_revision =
	controller->features().presentationRevision();
    controller->clearRenderRequest();
    transaction.feature = feature;
    transaction.source_path = replacement_source;
    transaction.edit_intent_id = NULL;
    transaction.edit_intent_role = NULL;
    transaction.points = NULL;
    transaction.commands = NULL;
    transaction.point_count = 0;
    transaction.source_revision = 0;
    transaction.inputs_revision = 0;
    transaction.color_valid = 0;
    if (!ged_view_edit_transaction_apply(view_ctx, &transaction, NULL) ||
	controller->features().presentationRevision() != no_op_revision ||
	controller->isRenderRequested())
	FAIL("unchanged edit event should be an exact publication no-op");

    transaction.event = GED_VIEW_EDIT_PREVIEW_COMMIT;
    if (!ged_view_edit_transaction_apply(view_ctx, &transaction, NULL) ||
	controller->features().exists(feature_name,
	    BOBOL_FEATURE_SCOPE_LOCAL))
	FAIL("terminal edit transaction should retire its retained preview");

    static const char primitive_name[] = "edit-transaction::primitive";
    static const char primitive_source[] = "box.s";
    static const char primitive_intent[] = "edit-transaction::scale";
    static const char primitive_role[] = "scale-handle";
    static const uint32_t primitive_source_revision = 61;
    static const uint32_t primitive_inputs_revision = 62;
    struct directory *box_dp = db_lookup(gedp->dbip, primitive_source,
	LOOKUP_QUIET);
    struct rt_db_internal box_internal;
    RT_DB_INTERNAL_INIT(&box_internal);
    if (!box_dp || rt_db_get_internal(&box_internal, box_dp, gedp->dbip,
	    NULL) < 0)
	FAIL("primitive edit transaction should load its source geometry");

    mat_t primitive_matrix;
    MAT_IDN(primitive_matrix);
    MAT_DELTAS(primitive_matrix, 6.0, 0.0, 0.0);
    struct ged_view_edit_transaction primitive_transaction =
	ged_view_edit_transaction_default();
    primitive_transaction.event = GED_VIEW_EDIT_PREVIEW_BEGIN;
    primitive_transaction.feature_name = primitive_name;
    primitive_transaction.source_path = primitive_source;
    primitive_transaction.edit_intent_id = primitive_intent;
    primitive_transaction.edit_intent_role = primitive_role;
    primitive_transaction.dbip = gedp->dbip;
    primitive_transaction.internal = &box_internal;
    primitive_transaction.matrix = primitive_matrix;
    primitive_transaction.source_revision = primitive_source_revision;
    primitive_transaction.inputs_revision = primitive_inputs_revision;
    controller->clearRenderRequest();
    ged_view_edit_ref primitive = GED_VIEW_EDIT_REF_NULL;
    {
	edit_preview_publication_observer observer(controller, view_ctx,
	    primitive_name, primitive_source, primitive_intent,
	    primitive_role, NULL, NULL, 0, primitive_source_revision,
	    primitive_inputs_revision, NULL, true);
	if (!ged_view_edit_transaction_apply(view_ctx, &primitive_transaction,
		&primitive) || ged_view_edit_ref_is_null(primitive) ||
	    !observer.complete() || observer.partial ||
	    observer.graphCallbacks != 1 || observer.frameCallbacks != 1) {
	    rt_db_free_internal(&box_internal);
	    FAIL("primitive edit transaction exposed a partial retained publication");
	}
    }
    rt_db_free_internal(&box_internal);

    BObolFeatureRecord primitive_record;
    const BObolFeatureHandle primitive_handle =
	controller->features().find(primitive_name, BOBOL_FEATURE_SCOPE_LOCAL);
    if (!primitive_handle.isValid() ||
	!controller->features().record(primitive_handle, primitive_record))
	FAIL("primitive edit transaction should retain its complete preview");
    bool transformed = false;
    for (const SbVec3f &point : primitive_record.points) {
	if (point[X] > 4.5f) {
	    transformed = true;
	    break;
	}
    }
    if (!transformed)
	FAIL("primitive edit transaction should prepare transformed geometry");

    primitive_transaction.event = GED_VIEW_EDIT_PREVIEW_CANCEL;
    primitive_transaction.internal = NULL;
    if (!ged_view_edit_transaction_apply(view_ctx, &primitive_transaction,
	    NULL) || controller->features().exists(primitive_name,
		BOBOL_FEATURE_SCOPE_LOCAL))
	FAIL("terminal primitive edit transaction should retire its preview");

    std::puts("PASS GED edit transaction publication: atomic point/primitive insertion, replacement, exact no-op and terminal retirement");
    return 0;
}

static int
exercise_faceplate_composition_publication(struct ged *gedp,
	struct ged_view_context *view_ctx, BObolViewController *controller)
{
    struct bv *view = bv_context_view(
	reinterpret_cast<struct bv_context *>(view_ctx));
    if (!gedp || !view_ctx || !controller || !view)
	FAIL("faceplate composition fixture requires an attached view controller");

    struct bv_other_state center = BV_OTHER_STATE_INIT;
    center.gos_draw = 1;
    VSET(center.gos_line_color, 255, 240, 32);
    struct bv_interactive_rect_state rectangle =
	BV_INTERACTIVE_RECT_STATE_INIT;
    rectangle.draw = 1;
    rectangle.line_width = 2;
    rectangle.x = -0.5;
    rectangle.y = -0.25;
    rectangle.width = 1.0;
    rectangle.height = 0.5;
    VSET(rectangle.color, 32, 192, 255);
    if (!bv_center_dot_state_set(view, &center) ||
	!bv_interactive_rect_state_set(view, &rectangle))
	FAIL("faceplate composition fixture could not set view state");

    controller->clearRenderRequest();
    faceplate_composition_observer observer(controller, view_ctx, true, false);
    SoNodeSensor scene_sensor(
	faceplate_composition_observer::graph_changed, &observer);
    SoNodeSensor overlay_sensor(
	faceplate_composition_observer::graph_changed, &observer);
    scene_sensor.setPriority(0);
    overlay_sensor.setPriority(0);
    scene_sensor.attach(controller->getSceneRoot());
    overlay_sensor.attach(controller->getFramebufferOverlayRoot());
    controller->setFrameRequestCallback(
	faceplate_composition_observer::frame_requested, &observer);
    const int result = ged_view_faceplate_sync(gedp, view_ctx);
    controller->clearFrameRequestCallback(&observer);
    overlay_sensor.detach();
    scene_sensor.detach();

    if (result != BRLCAD_OK || !observer.complete() || observer.partial ||
	!observer.graph_callbacks || !observer.frame_callbacks) {
	fprintf(stderr,
	    "faceplate composition: result=%d complete=%d partial=%d graph=%zu frame=%zu revision=%" PRIu64 "/%" PRIu64 " center=%d rectangle=%d\n",
	    result, observer.complete() ? 1 : 0, observer.partial ? 1 : 0,
	    observer.graph_callbacks, observer.frame_callbacks,
	    observer.presentation_revision,
	    controller->features().presentationRevision(),
	    controller->features().exists("_faceplate/center_dot") ? 1 : 0,
	    controller->features().exists("_faceplate/interactive_rect") ? 1 : 0);
	FAIL("GED faceplate refresh exposed a partial multi-record publication");
    }

    SoNode *initial_center = controller->features().node(
	controller->features().find("_faceplate/center_dot"));
    VSET(center.gos_line_color, 255, 64, 32);
    rectangle.draw = 0;
    struct bv_other_state scale = BV_OTHER_STATE_INIT;
    scale.gos_draw = 1;
    VSET(scale.gos_line_color, 255, 255, 0);
    scale.gos_font_size = 18;
    if (!bv_center_dot_state_set(view, &center) ||
	!bv_interactive_rect_state_set(view, &rectangle) ||
	!bv_scale_overlay_state_set(view, &scale))
	FAIL("faceplate composition fixture could not set replacement state");

    controller->clearRenderRequest();
    faceplate_composition_observer replacement_observer(
	controller, view_ctx, false, true);
    SoNodeSensor replacement_scene_sensor(
	faceplate_composition_observer::graph_changed, &replacement_observer);
    SoNodeSensor replacement_overlay_sensor(
	faceplate_composition_observer::graph_changed, &replacement_observer);
    replacement_scene_sensor.setPriority(0);
    replacement_overlay_sensor.setPriority(0);
    replacement_scene_sensor.attach(controller->getSceneRoot());
    replacement_overlay_sensor.attach(controller->getFramebufferOverlayRoot());
    controller->setFrameRequestCallback(
	faceplate_composition_observer::frame_requested, &replacement_observer);
    const int replacement_result = ged_view_faceplate_sync(gedp, view_ctx);
    controller->clearFrameRequestCallback(&replacement_observer);
    replacement_overlay_sensor.detach();
    replacement_scene_sensor.detach();
    if (replacement_result != BRLCAD_OK ||
	!replacement_observer.complete() || replacement_observer.partial ||
	!replacement_observer.graph_callbacks ||
	!replacement_observer.frame_callbacks ||
	controller->features().node(controller->features().find(
	    "_faceplate/center_dot")) == initial_center)
	FAIL("GED faceplate replacement exposed a partial multi-record publication");

    std::puts("PASS GED faceplate composition publication: atomic insertion, replacement/removal, overlay metadata, graph and frame observers");
    return 0;
}

int
main(int argc, char **argv)
{
    if (exercise_typed_pick_result())
	return 1;

    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);
    /* The first deferred draw exercises failure after starting one worker. */
    bu_setenv("BOBOL_SOURCE_REALIZATION_WORKERS", "2", 1);
    char lcache[MAXPATHLEN] = {0};
    char cache_leaf[64] = {0};
    snprintf(cache_leaf, sizeof(cache_leaf), "ged_obol_draw_sync_cache_%d",
	bu_pid());
    bu_dir(lcache, MAXPATHLEN, BU_DIR_CURR, cache_leaf, NULL);
    bu_dirclear(lcache);
    bu_mkdir(lcache);
    bu_setenv("BU_DIR_CACHE", lcache, 1);
    bobol_init(NULL);

    const char *dbpath = "ged_obol_draw_sync_tmp.g";
    bu_file_delete(dbpath);
    if (!make_obol_sync_db(dbpath))
	FAIL("failed to create GED Obol draw-sync test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open GED Obol draw-sync test database");

    struct ged_view_context *initial_view_ctx = ged_draw_active_view_ctx(gedp);
    if (ged_draw_obol_scene_controller(gedp) ||
	    ged_draw_obol_controller(gedp) ||
	    ged_draw_ensure_root_attached(gedp) ||
	    (initial_view_ctx && ged_view_context_scene_attached(initial_view_ctx)))
	FAIL("GED without an Obol owner should expose an explicitly empty scene");

    if (!initial_view_ctx ||
	!ged_view_context_display_endpoint_ensure(initial_view_ctx))
	FAIL("GED should create an owned display endpoint for its active view");
    bobol_display_endpoint_t *initial_endpoint =
	ged_view_context_obol_endpoint_get(initial_view_ctx);
    BObolViewController *initial_view_controller = initial_endpoint ?
	static_cast<BObolViewController *>(
	    bobol_display_endpoint_controller(initial_endpoint)) : NULL;
    BObolSceneController *owned_scene =
	ged_draw_obol_scene_controller(gedp);
    BObolViewController *owned_controller = ged_draw_obol_controller(gedp);
    if (!initial_endpoint || !initial_view_controller || !owned_scene ||
	!owned_controller || !ged_draw_obol_scene_controller_owned(gedp) ||
	owned_controller->getSceneController() != owned_scene)
	FAIL("GED endpoint ensure should create one owned shared scene");

    if (argc > 1 && BU_STR_EQUAL(argv[1],
	    "faceplate-composition-publication")) {
	const int result = exercise_faceplate_composition_publication(gedp,
	    initial_view_ctx, initial_view_controller);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }
    if (argc > 1 && BU_STR_EQUAL(argv[1],
	    "edit-transaction-publication")) {
	const int result = exercise_edit_transaction_publication(
	    gedp, initial_view_ctx, initial_view_controller);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }
    if (argc > 1 && BU_STR_EQUAL(argv[1],
	    "production-flow-lifecycle")) {
	const char *lod_enable[3] = {"view", "lod", "1"};
	const char *mesh_enable[4] = {"view", "lod", "mesh", "1"};
	const char *bot_threshold[4] = {
	    "view", "lod", "bot_threshold", "0"
	};
	if (ged_exec_view(gedp, 3, lod_enable) != BRLCAD_OK ||
	    ged_exec_view(gedp, 4, mesh_enable) != BRLCAD_OK ||
	    ged_exec_view(gedp, 4, bot_threshold) != BRLCAD_OK)
	    FAIL("production lifecycle should enable automatic mesh LoD");
	bu_setenv("BOBOL_LOD_TASK_DELAY_MS", "75", 1);
	int result = exercise_delayed_mesh_camera_close(gedp, initial_view_ctx);
	if (!result) {
	    BObolViewController *reopened_controller =
		ged_bobol_view_controller(initial_view_ctx);
	    result = exercise_progressive_autoview_lifecycle(gedp,
		reopened_controller, initial_view_ctx, GED_DRAW_MODE_WIRE,
		"progressive_root.c", "progressive_root_async.c", 4, true);
	}
	bu_setenv("BOBOL_LOD_TASK_DELAY_MS", "0", 1);
	(void)ged_view_context_obol_endpoint_set(initial_view_ctx, NULL, 0);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }
    /* This is a broad scene/transaction lifecycle test, not a default-policy
     * performance test.  Leaving AUTO active made hundreds of unrelated draw
     * assertions start background realization and wait for transient states
     * that dedicated LoD tests already cover.  Keep general draws eager here;
     * the two progressive fixtures below explicitly request deferred leaf
     * expansion and therefore still exercise their intended provider paths. */
    const char *mesh_lod_off[4] = {"view", "lod", "mesh", "0"};
    const char *csg_lod_off[4] = {"view", "lod", "csg", "0"};
    if (ged_exec_view(gedp, 4, mesh_lod_off) != BRLCAD_OK ||
	ged_exec_view(gedp, 4, csg_lod_off) != BRLCAD_OK)
	FAIL("GED lifecycle test should isolate unrelated draws from automatic LoD");
    /* The bridge test below expects eager headless realization.  Detaching the
     * initial endpoint leaves the shared per-GED scene alive while removing
     * the per-view progressive provider. */
    if (!ged_view_context_obol_endpoint_set(initial_view_ctx, NULL, 0))
	FAIL("GED should detach an ensured endpoint without releasing its scene");

    if (argc > 1 && BU_STR_EQUAL(argv[1], "presentation-publication")) {
	const int result = exercise_presentation_publication(gedp, owned_scene);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }

    if (argc > 1 && BU_STR_EQUAL(argv[1], "source-display-publication")) {
	const int result = exercise_source_display_publication(gedp, owned_scene);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }

    if (argc > 1 && BU_STR_EQUAL(argv[1], "deferred-selection-acceptance")) {
	const int result = exercise_deferred_selection_acceptance(gedp,
	    owned_scene);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }

    if (argc > 1 && BU_STR_EQUAL(argv[1], "subtract-style-publication")) {
	const int result = exercise_progressive_occurrence_and_boolean_identity(
	    gedp, owned_scene);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }

    if (argc > 1 && BU_STR_EQUAL(argv[1],
	    "deferred-result-acceptance")) {
	const int result = exercise_progressive_occurrence_and_boolean_identity(
	    gedp, owned_scene);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }

    if (argc > 1 && BU_STR_EQUAL(argv[1],
	    "deferred-source-replacement")) {
	const int result = exercise_deferred_source_replacement(gedp,
	    owned_scene);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }

    if (argc > 1 && BU_STR_EQUAL(argv[1], "database-rename-publication")) {
	const int result = exercise_database_rename_publication(gedp,
	    owned_scene);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }

    if (argc > 1 && BU_STR_EQUAL(argv[1], "redraw-erase-publication")) {
	const int result = exercise_redraw_erase_publication(gedp,
	    owned_scene);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(lcache);
	return result;
    }

    if (exercise_base_source_promotion(gedp, owned_scene))
	return 1;

    if (exercise_deferred_proxy_transition(gedp, owned_scene))
	return 1;

    if (exercise_retained_current_publication(gedp, owned_scene))
	return 1;

    if (exercise_source_summary_copy_publication(gedp, owned_scene))
	return 1;

    if (exercise_presentation_publication(gedp, owned_scene))
	return 1;

    if (exercise_source_display_publication(gedp, owned_scene))
	return 1;

    if (exercise_source_record_publication(gedp, owned_scene))
	return 1;

    if (exercise_evaluated_region_publication(gedp, owned_scene))
	return 1;

    if (exercise_metadata_effects(gedp, owned_scene))
	return 1;

    if (exercise_redraw_appearance(gedp, owned_scene))
	return 1;

    if (exercise_redraw_erase_publication(gedp, owned_scene))
	return 1;

    if (exercise_draw_frontier(gedp, owned_scene))
	return 1;

    const char *draw_box[2] = {"draw", "box.s"};
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK)
	FAIL("real GED draw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(owned_scene, "box.s"))
	FAIL("endpoint-owned Obol scene should mirror GED draw command");
    SoBRLDatabaseSource *box_source = source_for_path(owned_scene, "box.s");
    if (!box_source ||
	    box_source->realizationStatus.getValue() !=
	    SoBRLDatabaseSource::REALIZED ||
	    !box_source->hasRealizedWireGeometry() ||
	    box_source->getRealizedShapeCount() != 0)
	FAIL("mirrored GED draw should realize Obol wire geometry");

    const char *draw_ball[2] = {"draw", "ball.s"};
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("second real GED draw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene should retain both endpoint-routed draws");
    const int owned_source_count_before_binunif =
	owned_scene->getDatabaseSourceCount();
    const char *draw_binunif[2] = {"draw", "payload.binunif"};
    if (ged_exec_draw(gedp, 2, draw_binunif) != BRLCAD_OK)
	FAIL("real GED draw of non-drawable binunif should report command success");
    if (owned_scene->getDatabaseSourceCount() !=
	    owned_source_count_before_binunif ||
	    source_for_path(owned_scene, "payload.binunif"))
	FAIL("Obol draw bridge should not publish non-drawable binunif sources");
    if (exercise_duplicate_occurrence_pick_identity(gedp, owned_scene))
	return 1;
    if (exercise_progressive_occurrence_and_boolean_identity(gedp,
	    owned_scene))
	return 1;
    if (exercise_deferred_source_replacement(gedp, owned_scene))
	return 1;

    struct ged_view_context *feature_view_ctx = ged_view_active_ctx(gedp);
    if (!feature_view_ctx ||
	!ged_view_context_display_endpoint_ensure(feature_view_ctx) ||
	    !ged_draw_obol_view_context_feature_store_active(feature_view_ctx))
	FAIL("GED active view endpoint should expose its Obol feature store");
    bobol_display_endpoint_t *feature_endpoint =
	ged_view_context_obol_endpoint_get(feature_view_ctx);
    BObolViewController *feature_view_controller = feature_endpoint ?
	static_cast<BObolViewController *>(
	    bobol_display_endpoint_controller(feature_endpoint)) : NULL;
    if (!feature_view_controller)
	FAIL("GED active view endpoint should carry an Obol controller");
    point_t feature_points[2] = {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}};
    int feature_cmds[2] = {
	GED_DRAW_VIEW_LINE_MOVE,
	GED_DRAW_VIEW_LINE_DRAW
    };
    struct ged_view_feature_style feature_style =
	ged_view_feature_style_default();
    feature_style.visible = 1;
    feature_style.color_valid = 1;
    feature_style.color[0] = 20;
    feature_style.color[1] = 40;
    feature_style.color[2] = 60;
    feature_style.line_width = 3;
    if (!test_line_replace(feature_view_ctx, "cap2::line", 0,
	    feature_points, feature_cmds, 2, &feature_style,
	    GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION,
	    GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT,
	    GED_VIEW_FEATURE_OVERLAY_ORDER_MODEL) ||
	    !owned_controller->features().exists("cap2::line"))
	FAIL("GED feature line replacement should publish into the owned Obol feature store");
    struct ged_view_feature_summary feature_summary =
	ged_view_feature_summary_default();
    if (!ged_view_feature_get_summary(feature_view_ctx,
	    "cap2::line", &feature_summary) ||
	    !feature_summary.exists ||
	    feature_summary.geometry_command_count != 2 ||
	    feature_summary.color[0] != 20 ||
	    feature_summary.color[1] != 40 ||
	    feature_summary.color[2] != 60)
	FAIL("GED feature summary should read the owned Obol feature store");
    point_t *copied_points = NULL;
    size_t copied_point_count = 0;
    if (!ged_view_feature_points_copy(feature_view_ctx,
	    "cap2::line", &copied_points, &copied_point_count) ||
	    copied_point_count != 2 ||
	    !NEAR_EQUAL(copied_points[1][X], 2.0, SMALL_FASTF))
	FAIL("GED feature point readback should copy owned Obol line geometry");
    bu_free(copied_points, "GED Obol feature test copied points");
    int copied_cmd = 0;
    if (!ged_view_feature_line_command_at(feature_view_ctx,
	    "cap2::line", 1, &copied_cmd) ||
	    copied_cmd != GED_DRAW_VIEW_LINE_DRAW)
	FAIL("GED feature command readback should copy owned Obol line commands");
    point_t appended_points[3] = {
	{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}
    };
    int appended_commands[3] = {
	GED_DRAW_VIEW_LINE_MOVE, GED_DRAW_VIEW_LINE_DRAW,
	GED_DRAW_VIEW_LINE_DRAW
    };
    if (!test_line_replace(feature_view_ctx, "cap2::line", 0,
	    appended_points, appended_commands, 3, &feature_style,
	    GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION,
	    GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT,
	    GED_VIEW_FEATURE_OVERLAY_ORDER_MODEL) ||
	    !ged_view_feature_points_copy(feature_view_ctx,
		"cap2::line", &copied_points, &copied_point_count) ||
	    copied_point_count != 3)
	FAIL("GED line append should mutate owned Obol line geometry");
    bu_free(copied_points, "GED Obol feature test appended points");
    if (!ged_view_feature_visible_set(feature_view_ctx,
	    "cap2::line", 0) ||
	    ged_view_feature_visible(feature_view_ctx,
		"cap2::line") != 0)
	FAIL("GED feature visibility should mutate owned Obol feature style");
    struct ged_view_feature_style changed_line_style =
	ged_view_feature_style_default();
    changed_line_style.color_valid = 1;
    VSET(changed_line_style.color, 90, 80, 70);
    changed_line_style.line_width = 5;
    if (!ged_view_feature_style_apply(feature_view_ctx, "cap2::line",
	    &changed_line_style, 0))
	FAIL("GED line style setters should mutate owned Obol feature style");
    struct ged_view_feature_style line_style = ged_view_feature_style_default();
    if (!ged_view_feature_style_get(feature_view_ctx,
	    "cap2::line", &line_style) || !line_style.color_valid ||
	    line_style.color[0] != 90 || line_style.color[1] != 80 ||
	    line_style.color[2] != 70 || line_style.line_width != 5)
	FAIL("GED line style readback should read owned Obol feature style");

    point_t tcl_line_points[2] = {{0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}};
    int tcl_line_commands[2] = {
	GED_DRAW_VIEW_LINE_MOVE, GED_DRAW_VIEW_LINE_DRAW
    };
    struct ged_view_feature_style tcl_feature_style =
	ged_view_feature_style_default();
    tcl_feature_style.visible = 1;
    tcl_feature_style.color_valid = 1;
    VSET(tcl_feature_style.color, 101, 102, 103);
    tcl_feature_style.line_width = 4;
    if (!test_line_replace(feature_view_ctx, "cap2::tcl-line", 1,
	    tcl_line_points, tcl_line_commands, 2, &tcl_feature_style,
	    GED_VIEW_FEATURE_OVERLAY_CLASS_TCL_OVERLAY,
	    GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND,
	    GED_VIEW_FEATURE_OVERLAY_ORDER_POST_TRANSPARENT) ||
	    !feature_overlay_matches(feature_view_controller, "cap2::tcl-line",
		BObolOverlayClass::TclOverlay,
		BObolOverlayLifecycle::PerCommand,
		BObolOverlayOrder::PostTransparent))
	FAIL("GED Tcl line replacement should publish typed Obol overlay metadata");
    struct ged_view_feature_summary tcl_line_summary =
	ged_view_feature_summary_default();
    if (!ged_view_feature_get_summary(feature_view_ctx,
	    "cap2::tcl-line", &tcl_line_summary) ||
	    !tcl_line_summary.exists ||
	    !tcl_line_summary.is_overlay ||
	    tcl_line_summary.geometry_command_count != 2)
	FAIL("GED Tcl line summary should read typed Obol overlay feature state");

    point_t annotation_point = {5.0, 5.0, 0.0};
    int annotation_command = GED_DRAW_VIEW_LINE_MOVE;
    if (!test_line_replace(feature_view_ctx, "cap2::annotation", 1,
	    &annotation_point, &annotation_command, 1, NULL,
	    GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION,
	    GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT,
	    GED_VIEW_FEATURE_OVERLAY_ORDER_MODEL) ||
	    !feature_overlay_matches(feature_view_controller, "cap2::annotation",
		BObolOverlayClass::UserAnnotation,
		BObolOverlayLifecycle::Persistent,
		BObolOverlayOrder::Model))
	FAIL("GED model annotation creation should publish typed Obol overlay metadata");

    point_t polygon_points[2] = {{1.0, 0.0, 0.0}, {1.0, 2.0, 0.0}};
    int polygon_cmds[2] = {
	GED_DRAW_VIEW_LINE_MOVE,
	GED_DRAW_VIEW_LINE_DRAW
    };
    if (!test_line_replace(feature_view_ctx, "cap2::polygon-overlay", 1,
	    polygon_points, polygon_cmds, 2, &feature_style,
	    GED_VIEW_FEATURE_OVERLAY_CLASS_POLYGON_EDIT,
	    GED_VIEW_FEATURE_LIFECYCLE_PER_TOOL,
	    GED_VIEW_FEATURE_OVERLAY_ORDER_POST_TRANSPARENT) ||
	    !feature_overlay_matches(feature_view_controller,
		"cap2::polygon-overlay",
		BObolOverlayClass::PolygonEdit,
		BObolOverlayLifecycle::PerTool,
		BObolOverlayOrder::PostTransparent))
	FAIL("GED Tcl polygon replacement should publish typed Obol overlay metadata");

    struct ged_view_feature_label label = {};
    label.text = "cap2 label";
    VSET(label.point, 1.0, 2.0, 3.0);
    label.color_valid = 1;
    label.color[0] = 7;
    label.color[1] = 8;
    label.color[2] = 9;
    label.font_size = 18.0;
    if (!test_labels_replace(feature_view_ctx, "cap2::label", 0,
	    &label, 1, NULL) ||
	    !owned_controller->features().exists("cap2::label") ||
	    ged_view_feature_label_count(feature_view_ctx,
		"cap2::label") != 1)
	FAIL("GED label replacement should publish into the owned Obol feature store");
    BObolFeatureHandle label_handle =
	owned_controller->features().find("cap2::label");
    BObolFeatureRecord label_record;
    if (!label_handle.isValid() ||
	    !owned_controller->features().record(label_handle,
		label_record) ||
	    label_record.labels.size() != 1 ||
	    fabs(label_record.labels[0].fontSize - 18.0f) > 0.001f)
	FAIL("GED label replacement should preserve explicit Obol font size");
    struct bu_vls label_text = BU_VLS_INIT_ZERO;
    point_t label_point = VINIT_ZERO;
    unsigned char label_rgb[3] = {0, 0, 0};
    if (!ged_view_feature_label_copy(feature_view_ctx, "cap2::label",
	    0, &label_text, label_point, label_rgb) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&label_text), "cap2 label") ||
	    !NEAR_EQUAL(label_point[X], 1.0, SMALL_FASTF) ||
	    label_rgb[0] != 7 ||
	    label_rgb[1] != 8 ||
	    label_rgb[2] != 9) {
	bu_vls_free(&label_text);
	FAIL("GED label copy should read owned Obol label data");
    }
    bu_vls_free(&label_text);
    point_t moved_label_point = {4.0, 5.0, 6.0};
    VMOVE(label.point, moved_label_point);
    if (!test_labels_replace(feature_view_ctx, "cap2::label", 0,
	    &label, 1, NULL) ||
	    !ged_view_feature_label_copy(feature_view_ctx, "cap2::label",
		0, NULL, label_point, NULL) ||
	    !NEAR_EQUAL(label_point[X], 4.0, SMALL_FASTF) ||
	    !NEAR_EQUAL(label_point[Y], 5.0, SMALL_FASTF) ||
	    !NEAR_EQUAL(label_point[Z], 6.0, SMALL_FASTF))
	FAIL("GED label point mutation should update owned Obol label data");

    point_t created_label_point = {1.0, 1.0, 1.0};
    point_t created_label_target = {0.0, 0.0, 0.0};
    struct ged_view_feature_label created_label = {};
    created_label.text = "created";
    VMOVE(created_label.point, created_label_point);
    created_label.line_flag = 1;
    VMOVE(created_label.target, created_label_target);
    struct ged_view_feature_style created_label_style =
	ged_view_feature_style_default();
    created_label_style.visible = 1;
    created_label_style.color_valid = 1;
    VSET(created_label_style.color, 255, 255, 0);
    if (!test_labels_replace(feature_view_ctx, "cap2::created-label", 0,
	    &created_label, 1, &created_label_style) ||
	    !ged_view_feature_label_copy(feature_view_ctx,
		"cap2::created-label", 0, NULL, NULL, label_rgb) ||
	    label_rgb[0] != 255 ||
	    label_rgb[1] != 255 ||
	    label_rgb[2] != 0)
	FAIL("GED label create should preserve the legacy yellow feature color in Obol");

    point_t arrow_points[2] = {{0.0, 0.0, 0.0}, {0.0, 3.0, 0.0}};
    struct ged_view_feature_style arrow_style = feature_style;
    arrow_style.arrow = 1;
    arrow_style.arrow_tip_length = 0.25;
    arrow_style.arrow_tip_width = 0.5;
    if (!test_arrows_replace(feature_view_ctx, "cap2::arrow",
	    arrow_points, 2, &arrow_style) ||
	    !feature_view_controller->features().exists("cap2::arrow") ||
	    !feature_overlay_matches(feature_view_controller, "cap2::arrow",
		BObolOverlayClass::TclOverlay,
		BObolOverlayLifecycle::PerCommand,
		BObolOverlayOrder::PostTransparent))
	FAIL("GED arrow replacement should publish into the owned Obol feature store");
    fastf_t tip_length = 0.0;
    fastf_t tip_width = 0.0;
    struct ged_view_feature_style arrow_style_read =
	ged_view_feature_style_default();
    if (!ged_view_feature_style_get(feature_view_ctx, "cap2::arrow",
	    &arrow_style_read))
	FAIL("GED arrow style readback should read owned Obol feature style");
    tip_length = arrow_style_read.arrow_tip_length;
    tip_width = arrow_style_read.arrow_tip_width;
    if (
	    fabs(tip_length - 0.25) > 0.001 ||
	    fabs(tip_width - 0.5) > 0.001)
	FAIL("GED arrow tip readback should read owned Obol feature style");

    point_t axes_center = {1.0, 1.0, 1.0};
    struct ged_view_feature_style axes_style = ged_view_feature_style_default();
    axes_style.visible = 1;
    axes_style.color_valid = 1;
    VSET(axes_style.color, 11, 22, 33);
    axes_style.line_width = 2;
    if (!test_axes_replace(feature_view_ctx, "cap2::axes", 0,
	    &axes_center, 1, 4.0, &axes_style,
	    GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION,
	    GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT,
	    GED_VIEW_FEATURE_OVERLAY_ORDER_MODEL) ||
	    !owned_controller->features().exists("cap2::axes"))
	FAIL("GED axes creation should publish into the owned Obol feature store");
    point_t axes_readback = VINIT_ZERO;
    fastf_t axes_size = 0.0;
    struct ged_view_feature_style axes_style_read =
	ged_view_feature_style_default();
    if (!ged_view_feature_axes_copy(feature_view_ctx, "cap2::axes", 0,
	    axes_readback, &axes_size) ||
	    !ged_view_feature_style_get(feature_view_ctx, "cap2::axes",
		&axes_style_read) ||
	    !NEAR_EQUAL(axes_readback[X], 1.0, SMALL_FASTF) ||
	    !NEAR_EQUAL(axes_size, 4.0, SMALL_FASTF) ||
	    axes_style_read.line_width != 2 ||
	    axes_style_read.color[0] != 11 ||
	    axes_style_read.color[1] != 22 ||
	    axes_style_read.color[2] != 33)
	FAIL("GED axes readback should read owned Obol axes state");

    point_t tcl_axes_centers[1] = {{2.0, 2.0, 0.0}};
    if (!test_axes_replace(feature_view_ctx, "cap2::tcl-axes", 1,
	    tcl_axes_centers, 1, 2.5, &feature_style,
	    GED_VIEW_FEATURE_OVERLAY_CLASS_TCL_OVERLAY,
	    GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND,
	    GED_VIEW_FEATURE_OVERLAY_ORDER_POST_TRANSPARENT) ||
	    !feature_overlay_matches(feature_view_controller, "cap2::tcl-axes",
		BObolOverlayClass::TclOverlay,
		BObolOverlayLifecycle::PerCommand,
		BObolOverlayOrder::PostTransparent))
	FAIL("GED Tcl axes replacement should publish typed Obol overlay metadata");

    point_t face_points[4] = {
	{0.0, 0.0, 0.0},
	{1.0, 0.0, 0.0},
	{1.0, 1.0, 0.0},
	{0.0, 1.0, 0.0}
    };
    int face_indices[6] = {0, 1, 2, 0, 2, 3};
    if (!test_mesh_replace(feature_view_ctx, "cap2::mesh", face_points,
	    4, face_indices, 6, &feature_style) ||
	    !owned_controller->features().exists("cap2::mesh"))
	FAIL("GED indexed-face replacement should publish into the owned Obol feature store");

    BObolFeatureHandle patch_mesh_handle =
	owned_controller->features().find("cap2::mesh");
    SoNode *patch_mesh_node = owned_controller->features().node(
	patch_mesh_handle);
    int patch_index = 2;
    point_t patch_point = {1.0, 1.0, 2.0};
    if (!patch_mesh_node ||
	!patch_mesh_node->isOfType(SoBRLMeshShape::getClassTypeId()) ||
	!ged_view_feature_indexed_face_points_update(feature_view_ctx,
	    "cap2::mesh", &patch_index, &patch_point, 1) ||
	owned_controller->features().node(patch_mesh_handle) != patch_mesh_node ||
	static_cast<SoBRLMeshShape *>(patch_mesh_node)->point[2] !=
	    SbVec3f(1.0f, 1.0f, 2.0f))
	FAIL("GED indexed-face point patch should update retained geometry in place");
    patch_index = 99;
    VSET(patch_point, 9.0, 9.0, 9.0);
    if (ged_view_feature_indexed_face_points_update(feature_view_ctx,
	    "cap2::mesh", &patch_index, &patch_point, 1) ||
	static_cast<SoBRLMeshShape *>(patch_mesh_node)->point[2] !=
	    SbVec3f(1.0f, 1.0f, 2.0f))
	FAIL("GED invalid indexed-face point patch should be atomic");

    struct bg_line_layer_builder *diagnostic_builder =
	bg_line_layer_builder_create();
    if (!diagnostic_builder)
	FAIL("diagnostic line-layer builder should allocate");
    point_t diagnostic_a = {0.0, 0.0, 0.0};
    point_t diagnostic_b = {0.5, 0.5, 0.0};
    if (!bg_line_layer_builder_add(diagnostic_builder, 255, 0, 0,
		diagnostic_a, BG_GEOMETRY_LINE_MOVE) ||
	    !bg_line_layer_builder_add(diagnostic_builder, 255, 0, 0,
		diagnostic_b, BG_GEOMETRY_LINE_DRAW)) {
	bg_line_layer_builder_free(diagnostic_builder);
	FAIL("diagnostic line-layer builder should accept test geometry");
    }
    if (!test_line_builder_replace(feature_view_ctx, "cap2::diagnostic",
	    diagnostic_builder)) {
	bg_line_layer_builder_free(diagnostic_builder);
	FAIL("GED diagnostic line-layer replacement should publish into the owned Obol feature store");
    }
    bg_line_layer_builder_free(diagnostic_builder);
    BObolFeatureRecord diagnostic_record;
    BObolFeatureHandle diagnostic_handle =
	owned_controller->features().find("cap2::diagnostic");
    if (!diagnostic_handle.isValid() ||
	    !owned_controller->features().record(diagnostic_handle,
		diagnostic_record) ||
	    diagnostic_record.kind != BObolFeatureKind::LineLayer ||
	    !feature_overlay_matches(owned_controller, "cap2::diagnostic",
		BObolOverlayClass::Diagnostic,
		BObolOverlayLifecycle::PerCommand,
		BObolOverlayOrder::PostTransparent))
	FAIL("GED diagnostic line-layer replacement should stamp typed Obol diagnostic metadata");

    struct command_result_callback_state command_callback_state = {};
    struct ged_view_feature_batch_desc feature_batch_desc =
	ged_view_feature_batch_desc_default();
    feature_batch_desc.owner_id = "rtcheck";
    feature_batch_desc.owner_role = "command-result";
    feature_batch_desc.event_cb = command_event_cb;
    feature_batch_desc.event_cb_data = &command_callback_state;
    struct ged_view_feature_batch *aborted_scene =
	ged_view_feature_batch_begin(feature_view_ctx, &feature_batch_desc);
    if (!aborted_scene ||
	!ged_view_feature_batch_line_set_replace(aborted_scene,
	    "rtcheck::aborted", feature_points, feature_cmds, 2, NULL) ||
	owned_controller->features().exists("rtcheck::aborted"))
	FAIL("GED feature batches must not publish staged geometry before commit");
    ged_view_feature_batch_abort(aborted_scene);
    if (owned_controller->features().exists("rtcheck::aborted"))
	FAIL("GED feature-batch abort must discard staged geometry");

    struct ged_view_feature_batch *feature_batch =
	ged_view_feature_batch_begin(feature_view_ctx, &feature_batch_desc);
    if (!feature_batch)
	FAIL("GED feature-batch begin should create an Obol-backed publication context");
    struct ged_view_feature_style command_style =
	ged_view_feature_style_default();
    command_style.visible = 1;
    command_style.selectable = 0;
    struct ged_view_feature_line_layer command_layer =
	ged_view_feature_line_layer_default();
    command_layer.name = "rtcheck::overlaps/yellow";
    command_layer.points = feature_points;
    command_layer.commands = feature_cmds;
    command_layer.point_count = 2;
    struct ged_view_feature_metadata command_metadata[2] = {
	{"result.kind", "overlap"},
	{"result.count", "1"}
    };
    struct ged_view_feature_metadata command_primitive_metadata[2] = {
	{"overlap.objects", "box.s cone.s"},
	{"overlap.depth", "0.25"}
    };
    if (!ged_view_feature_batch_line_layers_replace(feature_batch,
	    "rtcheck::overlaps", &command_layer, 1, &command_style) ||
	    !ged_view_feature_batch_metadata_replace(feature_batch,
		"rtcheck::overlaps", command_metadata, 2) ||
	    !ged_view_feature_batch_primitive_metadata_replace(
		feature_batch, "rtcheck::overlaps", 0,
		command_primitive_metadata, 2))
	FAIL("GED feature-batch line-layer replacement should stage");
    if (owned_controller->features().exists("rtcheck::overlaps"))
	FAIL("GED feature-batch geometry and metadata must remain staged until commit");
    if (!ged_view_feature_batch_commit(feature_batch))
	FAIL("GED feature-batch line-layer replacement should commit");
    BObolFeatureHandle command_handle =
	owned_controller->features().find("rtcheck::overlaps",
		BOBOL_FEATURE_SCOPE_SHARED);
    BObolFeatureRecord command_record;
    if (!command_handle.isValid() ||
	    !owned_controller->features().record(command_handle,
		command_record) ||
	    command_record.scope != BObolFeatureScope::Shared ||
	    !BU_STR_EQUAL(command_record.owner.ownerId.getString(),
		"rtcheck") ||
	    !BU_STR_EQUAL(command_record.owner.ownerRole.getString(),
		"command-result") ||
	    !command_record.style.hasSelectable ||
	    command_record.style.selectable ||
	    command_record.metadata.size() != 2 ||
	    command_record.primitiveMetadata.size() != 1 ||
	    command_record.primitiveMetadata[0].primitiveIndex != 0 ||
	    command_record.primitiveMetadata[0].metadata.size() != 2 ||
	    !BU_STR_EQUAL(command_record.metadata[0].key.getString(),
		"result.kind") ||
	    !BU_STR_EQUAL(command_record.metadata[0].value.getString(),
		"overlap") ||
	    !BU_STR_EQUAL(command_record.metadata[1].key.getString(),
		"result.count") ||
	    !BU_STR_EQUAL(command_record.metadata[1].value.getString(), "1") ||
	    !feature_overlay_matches(owned_controller, "rtcheck::overlaps",
		BObolOverlayClass::CommandResult,
		BObolOverlayLifecycle::PerCommand,
		BObolOverlayOrder::PostTransparent))
	FAIL("GED feature-batch result should be shared, owned, selectable-aware command content");
    if (command_callback_state.accepted_count < 2 ||
	    command_callback_state.updated_count < 2 ||
	    !command_callback_state.saw_line_layers_update ||
	    !command_callback_state.saw_metadata_update ||
	    !command_callback_state.saw_primitive_metadata_update ||
	    command_callback_state.line_layers_feature_id != command_handle.id ||
	    command_callback_state.metadata_feature_id != command_handle.id ||
	    command_callback_state.primitive_metadata_feature_id !=
		command_handle.id)
	FAIL("GED feature-batch result callback should report line-layer and metadata updates with feature handles");
    struct bu_vls primitive_key = BU_VLS_INIT_ZERO;
    struct bu_vls primitive_value = BU_VLS_INIT_ZERO;
    if (ged_view_feature_primitive_metadata_count(
		feature_view_ctx, "rtcheck::overlaps", 0) != 2 ||
	    !ged_view_feature_primitive_metadata_copy(
		feature_view_ctx, "rtcheck::overlaps", 0, 0,
		&primitive_key, &primitive_value) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&primitive_key), "overlap.objects") ||
	    !BU_STR_EQUAL(bu_vls_cstr(&primitive_value), "box.s cone.s")) {
	bu_vls_free(&primitive_key);
	bu_vls_free(&primitive_value);
	FAIL("GED feature-batch primitive metadata should read back through neutral view APIs");
    }
    bu_vls_free(&primitive_key);
    bu_vls_free(&primitive_value);

    struct bu_vls resolved_feature = BU_VLS_INIT_ZERO;
    int resolved_primitive = -1;
    int primitive_index = -1;
    if (!ged_view_feature_pick_primitive_resolve(
		feature_view_ctx, "rtcheck::overlaps/yellow", 0, 1, 1,
		&resolved_feature, &resolved_primitive) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&resolved_feature),
		"rtcheck::overlaps") ||
	    resolved_primitive != 0 ||
	    ged_view_feature_selection_count(
		feature_view_ctx, "rtcheck::overlaps") != 1 ||
	    ged_view_feature_highlight_count(
		feature_view_ctx, "rtcheck::overlaps") != 1 ||
	    !ged_view_feature_selection_at(
		feature_view_ctx, "rtcheck::overlaps", 0, &primitive_index) ||
	    primitive_index != 0 ||
	    !ged_view_feature_highlight_at(
		feature_view_ctx, "rtcheck::overlaps", 0, &primitive_index) ||
	    primitive_index != 0) {
	bu_vls_free(&resolved_feature);
	FAIL("GED feature-batch child primitive picks should resolve and set parent primitive state");
    }
    bu_vls_free(&resolved_feature);
    if (ged_view_feature_pick_primitive_resolve(
	    feature_view_ctx, "rtcheck::overlaps/yellow", 1, 0, 0, NULL,
	    &resolved_primitive))
	FAIL("GED feature-batch child primitive resolver should reject out-of-range picks");
    struct ged_view_feature_summary command_summary =
	ged_view_feature_summary_default();
    if (!ged_view_feature_get_summary(feature_view_ctx,
	    "rtcheck::overlaps", &command_summary) ||
	    command_summary.primitive_metadata_count != 1 ||
	    command_summary.selected_primitive_count != 1 ||
	    command_summary.highlighted_primitive_count != 1)
	FAIL("GED feature-batch result summary should report primitive state");
    if (!owned_controller->features().record(command_handle,
		command_record) ||
	    command_record.selectedPrimitives.size() != 1 ||
	    command_record.highlightedPrimitives.size() != 1)
	FAIL("GED feature-batch result record should preserve primitive state");
    SoBRLVListShape *command_vlist = first_feature_vlist(
	    owned_controller->features().node(command_handle));
    if (!command_vlist ||
	    command_vlist->selectedPrimitive.getNum() != 1 ||
	    command_vlist->selectedPrimitive[0] != 0 ||
	    command_vlist->highlightedPrimitive.getNum() != 1 ||
	    command_vlist->highlightedPrimitive[0] != 0)
	FAIL("GED feature-batch result primitive state should reach realized Coin VLIST");

    feature_batch = ged_view_feature_batch_begin(feature_view_ctx,
	    &feature_batch_desc);
    if (!feature_batch ||
	    ged_view_feature_batch_remove_prefix(feature_batch,
		"rtcheck::") != 1 ||
	    !ged_view_feature_batch_commit(feature_batch) ||
	    owned_controller->features().exists("rtcheck::overlaps"))
	FAIL("GED feature-batch remove-prefix should remove owned shared command results");
    if (command_callback_state.removed_count < 1 ||
	    !command_callback_state.saw_remove_prefix)
	FAIL("GED feature-batch result callback should report owner-scoped feature removal");

    feature_batch_desc.generation = 10;
    struct ged_view_feature_batch *stale_scene =
	ged_view_feature_batch_begin(feature_view_ctx, &feature_batch_desc);
    if (!stale_scene)
	FAIL("GED feature-batch stale generation test should create old scene");

    point_t latest_points[2] = {
	{0.0, 0.0, 0.0},
	{0.0, 2.0, 0.0}
    };
    struct ged_view_feature_line_layer latest_layer =
	ged_view_feature_line_layer_default();
    latest_layer.name = "rtcheck::generation/latest";
    latest_layer.points = latest_points;
    latest_layer.commands = feature_cmds;
    latest_layer.point_count = 2;
    feature_batch_desc.generation = 11;
    feature_batch = ged_view_feature_batch_begin(feature_view_ctx,
	    &feature_batch_desc);
    if (!feature_batch ||
	    !ged_view_feature_batch_line_layers_replace(feature_batch,
		"rtcheck::generation", &latest_layer, 1, &command_style) ||
	    !ged_view_feature_batch_commit(feature_batch))
	FAIL("GED feature-batch latest generation should publish");

    point_t stale_points[2] = {
	{0.0, 0.0, 0.0},
	{0.0, 3.0, 0.0}
    };
    struct ged_view_feature_line_layer stale_layer =
	ged_view_feature_line_layer_default();
    stale_layer.name = "rtcheck::generation/stale";
    stale_layer.points = stale_points;
    stale_layer.commands = feature_cmds;
    stale_layer.point_count = 2;
    if (ged_view_feature_batch_line_layers_replace(stale_scene,
	    "rtcheck::generation", &stale_layer, 1, &command_style) ||
	    ged_view_feature_batch_commit(stale_scene))
	FAIL("GED feature-batch stale generation should be rejected");
    if (command_callback_state.failed_count < 2 ||
	    !command_callback_state.saw_stale_failure ||
	    !command_callback_state.saw_commit_failure)
	FAIL("GED feature-batch result callback should report stale generation rejection");

    command_handle = owned_controller->features().find("rtcheck::generation",
	    BOBOL_FEATURE_SCOPE_SHARED);
    if (!command_handle.isValid() ||
	    !owned_controller->features().record(command_handle,
		command_record) ||
	    command_record.owner.generation != 11 ||
	    command_record.points.size() != 2 ||
	    !NEAR_EQUAL(command_record.points[1][1], 2.0f, SMALL_FASTF))
	FAIL("GED feature-batch stale generation should not replace latest result");
    feature_batch = ged_view_feature_batch_begin(feature_view_ctx,
	    &feature_batch_desc);
    if (!feature_batch ||
	    ged_view_feature_batch_remove_prefix(feature_batch,
		"rtcheck::generation") != 1 ||
	    !ged_view_feature_batch_commit(feature_batch) ||
	    owned_controller->features().exists("rtcheck::generation"))
	FAIL("GED feature-batch generation cleanup should remove latest result");

    feature_batch_desc.generation = 0;
    feature_batch_desc.event_cb = NULL;
    feature_batch_desc.event_cb_data = NULL;

    FILE *nirt_plot = tmpfile();
    if (!nirt_plot)
	FAIL("NIRT/qray feature-batch uplot test should create a temporary plot stream");
    int old_plot_mode = pl_getOutputMode();
    pl_setOutputMode(PL_OUTPUT_MODE_BINARY);
    point_t nirt_a = VINIT_ZERO;
    point_t nirt_b = {1.0, 0.0, 0.0};
    pl_color(nirt_plot, 0, 255, 255);
    pdv_3line(nirt_plot, nirt_a, nirt_b);
    pl_setOutputMode(old_plot_mode);
    rewind(nirt_plot);
    if (_ged_view_feature_batch_publish_uplot(gedp, nirt_plot,
	    "query_ray", 1.0, PL_OUTPUT_MODE_BINARY, "nirt",
	    "command-result", "query_ray", "query-ray", 0) != BRLCAD_OK) {
	fclose(nirt_plot);
	FAIL("NIRT/qray uplot import should publish through feature-batch ownership");
    }
    fclose(nirt_plot);
    BObolFeatureHandle nirt_handle =
	owned_controller->features().find("query_ray",
		BOBOL_FEATURE_SCOPE_SHARED);
    BObolFeatureRecord nirt_record;
    if (!nirt_handle.isValid() ||
	    !owned_controller->features().record(nirt_handle, nirt_record) ||
	    nirt_record.kind != BObolFeatureKind::LineLayer ||
	    nirt_record.scope != BObolFeatureScope::Shared ||
	    !BU_STR_EQUAL(nirt_record.owner.ownerId.getString(), "nirt") ||
	    !BU_STR_EQUAL(nirt_record.owner.ownerRole.getString(),
		"command-result") ||
	    nirt_record.layers.size() != 1 ||
	    nirt_record.points.size() != 2 ||
	    nirt_record.metadata.size() < 8 ||
	    !BU_STR_EQUAL(nirt_record.metadata[0].key.getString(),
		"result.feature") ||
	    !BU_STR_EQUAL(nirt_record.metadata[0].value.getString(),
		"query_ray") ||
	    !BU_STR_EQUAL(nirt_record.metadata[1].key.getString(),
		"result.format") ||
	    !BU_STR_EQUAL(nirt_record.metadata[1].value.getString(),
		"uplot-line-layers") ||
	    !BU_STR_EQUAL(nirt_record.metadata[5].key.getString(),
		"result.kind") ||
	    !BU_STR_EQUAL(nirt_record.metadata[5].value.getString(),
		"query-ray") ||
	    !BU_STR_EQUAL(nirt_record.metadata[6].key.getString(),
		"result.schema") ||
	    !BU_STR_EQUAL(nirt_record.metadata[6].value.getString(),
		"brlcad.nirt.query-ray.v1") ||
	    !BU_STR_EQUAL(nirt_record.metadata[7].key.getString(),
		"result.severity") ||
	    !BU_STR_EQUAL(nirt_record.metadata[7].value.getString(), "mixed") ||
	    nirt_record.primitiveMetadata.size() != 1 ||
	    nirt_record.primitiveMetadata[0].primitiveIndex != 0 ||
	    nirt_record.primitiveMetadata[0].metadata.size() != 10 ||
	    !BU_STR_EQUAL(nirt_record.primitiveMetadata[0].metadata[0].key.getString(),
		"result.schema") ||
	    !BU_STR_EQUAL(nirt_record.primitiveMetadata[0].metadata[0].value.getString(),
		"brlcad.nirt.query-ray.v1") ||
	    !BU_STR_EQUAL(nirt_record.primitiveMetadata[0].metadata[2].value.getString(),
		"partition") ||
	    !BU_STR_EQUAL(nirt_record.primitiveMetadata[0].metadata[3].value.getString(),
		"info") ||
	    !BU_STR_EQUAL(nirt_record.primitiveMetadata[0].metadata[5].key.getString(),
		"segment.start_mm") ||
	    !BU_STR_EQUAL(nirt_record.primitiveMetadata[0].metadata[7].key.getString(),
		"hit.entry_mm") ||
	    !BU_STR_EQUAL(nirt_record.primitiveMetadata[0].metadata[9].key.getString(),
		"nirt.partition.parity") ||
	    !BU_STR_EQUAL(nirt_record.primitiveMetadata[0].metadata[9].value.getString(),
		"odd") ||
	    !feature_overlay_matches(owned_controller, "query_ray",
		BObolOverlayClass::CommandResult,
		BObolOverlayLifecycle::PerCommand,
		BObolOverlayOrder::PostTransparent))
	FAIL("NIRT/qray feature-batch uplot result should be shared owned command content");
    feature_batch_desc.owner_id = "nirt";
    feature_batch = ged_view_feature_batch_begin(feature_view_ctx,
	    &feature_batch_desc);
    if (!feature_batch ||
	    ged_view_feature_batch_remove_prefix(feature_batch,
		"query_ray") != 1 ||
	    !ged_view_feature_batch_commit(feature_batch) ||
	    owned_controller->features().exists("query_ray"))
	FAIL("NIRT/qray feature-batch cleanup should remove owned shared command results");

    FILE *rtcheck_plot = tmpfile();
    if (!rtcheck_plot)
	FAIL("rtcheck feature-batch schema test should create a temporary plot stream");
    old_plot_mode = pl_getOutputMode();
    pl_setOutputMode(PL_OUTPUT_MODE_BINARY);
    pl_color(rtcheck_plot, 255, 255, 0);
    pdv_3line(rtcheck_plot, nirt_a, nirt_b);
    pl_setOutputMode(old_plot_mode);
    rewind(rtcheck_plot);
    if (_ged_view_feature_batch_publish_uplot(gedp, rtcheck_plot,
	    "rtcheck::schema", 1.0, PL_OUTPUT_MODE_BINARY, "rtcheck",
	    "command-result", "rtcheck::schema", "overlap", 17) != BRLCAD_OK) {
	fclose(rtcheck_plot);
	FAIL("rtcheck uplot import should publish versioned overlap metadata");
    }
    fclose(rtcheck_plot);
    BObolFeatureHandle rtcheck_schema_handle =
	owned_controller->features().find("rtcheck::schema",
		BOBOL_FEATURE_SCOPE_SHARED);
    BObolFeatureRecord rtcheck_schema_record;
    if (!rtcheck_schema_handle.isValid() ||
	!owned_controller->features().record(rtcheck_schema_handle,
	    rtcheck_schema_record) ||
	rtcheck_schema_record.owner.generation != 17 ||
	rtcheck_schema_record.metadata.size() < 9 ||
	!BU_STR_EQUAL(rtcheck_schema_record.metadata[7].value.getString(),
	    "brlcad.rtcheck.overlap.v1") ||
	!BU_STR_EQUAL(rtcheck_schema_record.metadata[8].value.getString(),
	    "error") ||
	rtcheck_schema_record.primitiveMetadata.size() != 1 ||
	!BU_STR_EQUAL(rtcheck_schema_record.primitiveMetadata[0].metadata[2].value.getString(),
	    "overlap") ||
	!BU_STR_EQUAL(rtcheck_schema_record.primitiveMetadata[0].metadata[3].value.getString(),
	    "error") ||
	!BU_STR_EQUAL(rtcheck_schema_record.primitiveMetadata[0].metadata[7].key.getString(),
	    "hit.entry_mm"))
	FAIL("rtcheck result should retain overlap schema, generation, severity, and hit metadata");
    feature_batch_desc.owner_id = "rtcheck";
    feature_batch_desc.generation = 17;
    feature_batch = ged_view_feature_batch_begin(feature_view_ctx,
	&feature_batch_desc);
    if (!feature_batch ||
	ged_view_feature_batch_remove_prefix(feature_batch,
	    "rtcheck::schema") != 1 ||
	!ged_view_feature_batch_commit(feature_batch))
	FAIL("rtcheck schema result should use owner-scoped cleanup");
    feature_batch_desc.generation = 0;

    struct bg_line_layer_builder *builder_publish =
	bg_line_layer_builder_create();
    if (!builder_publish)
	FAIL("feature-batch builder publish test should allocate a builder");
    point_t builder_a = {0.0, 0.0, 0.0};
    point_t builder_b = {0.0, 1.0, 0.0};
    if (!bg_line_layer_builder_add(builder_publish, 255, 255, 0,
		builder_a, BG_GEOMETRY_LINE_MOVE) ||
	    !bg_line_layer_builder_add(builder_publish, 255, 255, 0,
		builder_b, BG_GEOMETRY_LINE_DRAW)) {
	bg_line_layer_builder_free(builder_publish);
	FAIL("feature-batch builder publish test should accept line geometry");
    }
    if (_ged_view_feature_batch_publish_line_layer_builder(gedp,
	    "nmg::_helper_test", builder_publish, "nmg",
	    "command-result", "nmg::_helper_test", "nmg-test", 0) != BRLCAD_OK) {
	bg_line_layer_builder_free(builder_publish);
	FAIL("line-layer builder helper should publish through feature-batch ownership");
    }
    bg_line_layer_builder_free(builder_publish);
    BObolFeatureHandle builder_handle =
	owned_controller->features().find("nmg::_helper_test",
		BOBOL_FEATURE_SCOPE_SHARED);
    BObolFeatureRecord builder_record;
    if (!builder_handle.isValid() ||
	    !owned_controller->features().record(builder_handle,
		builder_record) ||
	    builder_record.kind != BObolFeatureKind::LineLayer ||
	    builder_record.scope != BObolFeatureScope::Shared ||
	    !BU_STR_EQUAL(builder_record.owner.ownerId.getString(), "nmg") ||
	    !BU_STR_EQUAL(builder_record.owner.ownerRole.getString(),
		"command-result") ||
	    builder_record.points.size() != 2 ||
	    builder_record.metadata.size() < 6 ||
	    !BU_STR_EQUAL(builder_record.metadata[1].key.getString(),
		"result.format") ||
	    !BU_STR_EQUAL(builder_record.metadata[1].value.getString(),
		"line-layer-builder") ||
	    !BU_STR_EQUAL(builder_record.metadata[5].key.getString(),
		"result.kind") ||
	    !BU_STR_EQUAL(builder_record.metadata[5].value.getString(),
		"nmg-test") ||
	    !feature_overlay_matches(owned_controller, "nmg::_helper_test",
		BObolOverlayClass::CommandResult,
		BObolOverlayLifecycle::PerCommand,
		BObolOverlayOrder::PostTransparent))
	FAIL("line-layer builder helper should preserve shared owned command-result metadata");
    feature_batch_desc.owner_id = "nmg";
    feature_batch = ged_view_feature_batch_begin(feature_view_ctx,
	    &feature_batch_desc);
    if (!feature_batch ||
	    ged_view_feature_batch_remove_prefix(feature_batch,
		"nmg::_helper_test") != 1 ||
	    !ged_view_feature_batch_commit(feature_batch) ||
	    owned_controller->features().exists("nmg::_helper_test"))
	FAIL("line-layer builder helper result should clean up by owner-scoped prefix");

    struct bg_line_layer_builder *builder_generation =
	bg_line_layer_builder_create();
    if (!builder_generation)
	FAIL("feature-batch builder generation test should create builder");
    point_t builder_latest_b = {0.0, 2.0, 0.0};
    if (!bg_line_layer_builder_add(builder_generation, 255, 255, 0,
		builder_a, BG_GEOMETRY_LINE_MOVE) ||
	    !bg_line_layer_builder_add(builder_generation, 255, 255, 0,
		builder_latest_b, BG_GEOMETRY_LINE_DRAW)) {
	bg_line_layer_builder_free(builder_generation);
	FAIL("feature-batch builder generation test should accept latest geometry");
    }
    if (_ged_view_feature_batch_publish_line_layer_builder(gedp,
	    "nmg::_helper_generation", builder_generation, "nmg",
	    "command-result", "nmg::_helper_generation", "nmg-generation",
	    22) != BRLCAD_OK) {
	bg_line_layer_builder_free(builder_generation);
	FAIL("line-layer builder helper should publish latest generation");
    }
    bg_line_layer_builder_free(builder_generation);

    builder_generation = bg_line_layer_builder_create();
    if (!builder_generation)
	FAIL("feature-batch builder stale generation test should create builder");
    point_t builder_stale_b = {0.0, 3.0, 0.0};
    if (!bg_line_layer_builder_add(builder_generation, 255, 255, 0,
		builder_a, BG_GEOMETRY_LINE_MOVE) ||
	    !bg_line_layer_builder_add(builder_generation, 255, 255, 0,
		builder_stale_b, BG_GEOMETRY_LINE_DRAW)) {
	bg_line_layer_builder_free(builder_generation);
	FAIL("feature-batch builder stale generation test should accept stale geometry");
    }
    if (_ged_view_feature_batch_publish_line_layer_builder(gedp,
	    "nmg::_helper_generation", builder_generation, "nmg",
	    "command-result", "nmg::_helper_generation", "nmg-generation",
	    21) == BRLCAD_OK) {
	bg_line_layer_builder_free(builder_generation);
	FAIL("line-layer builder helper should reject stale generation without fallback");
    }
    bg_line_layer_builder_free(builder_generation);

    builder_handle = owned_controller->features().find(
	    "nmg::_helper_generation", BOBOL_FEATURE_SCOPE_SHARED);
    if (!builder_handle.isValid() ||
	    !owned_controller->features().record(builder_handle,
		builder_record) ||
	    builder_record.owner.generation != 22 ||
	    builder_record.points.size() != 2 ||
	    !NEAR_EQUAL(builder_record.points[1][1], 2.0f, SMALL_FASTF) ||
	    builder_record.metadata.size() < 7 ||
	    !BU_STR_EQUAL(builder_record.metadata[5].value.getString(),
		"nmg-generation") ||
	    !BU_STR_EQUAL(builder_record.metadata[6].key.getString(),
		"result.generation") ||
	    !BU_STR_EQUAL(builder_record.metadata[6].value.getString(), "22"))
	FAIL("line-layer builder helper stale generation should not replace latest result");
    feature_batch_desc.owner_id = "nmg";
    feature_batch_desc.generation = 22;
    feature_batch = ged_view_feature_batch_begin(feature_view_ctx,
	    &feature_batch_desc);
    if (!feature_batch ||
	    ged_view_feature_batch_remove_prefix(feature_batch,
		"nmg::_helper_generation") != 1 ||
	    !ged_view_feature_batch_commit(feature_batch) ||
	    owned_controller->features().exists("nmg::_helper_generation"))
	FAIL("line-layer builder helper generation result should clean up by owner-scoped prefix");
    feature_batch_desc.generation = 0;

    ged_view_edit_ref ged_preview_ref =
	ged_view_edit_overlay_ensure(feature_view_ctx,
		"cap2::ged-preview", "cap2::ged-source.s");
    if (ged_view_edit_ref_is_null(ged_preview_ref))
	FAIL("GED feature overlay ensure should return an Obol feature ref");
    BObolFeatureHandle ged_preview_handle =
	feature_view_controller->features().find("cap2::ged-preview");
    BObolFeatureSummary ged_preview_summary;
    if (!ged_preview_handle.isValid() ||
	    !feature_view_controller->features().summary("cap2::ged-preview",
		ged_preview_summary) ||
	    !ged_preview_summary.exists ||
	    ged_preview_summary.kind != BObolFeatureKind::EditPreview ||
	    !ged_preview_summary.overlay.isOverlay ||
	    ged_preview_summary.overlay.overlayClass !=
		BObolOverlayClass::EditHandle ||
	    ged_preview_summary.overlay.lifecycle !=
		BObolOverlayLifecycle::PerTool ||
	    ged_preview_summary.overlay.order !=
		BObolOverlayOrder::PostTransparent ||
	    ged_preview_summary.overlay.ownerToken != feature_view_ctx ||
	    !BU_STR_EQUAL(ged_preview_summary.overlay.sourcePath.getString(),
		"cap2::ged-source.s"))
	FAIL("GED feature overlay API should publish typed Obol edit-preview metadata");
    point_t ged_preview_points[3] = {
	{0.0, 0.0, 0.0},
	{1.0, 0.0, 0.0},
	{1.0, 1.0, 0.0}
    };
    int ged_preview_cmds[3] = {
	GED_DRAW_VIEW_LINE_MOVE,
	GED_DRAW_VIEW_LINE_DRAW,
	GED_DRAW_VIEW_LINE_DRAW
    };
    if (!ged_view_edit_points_replace(feature_view_ctx, ged_preview_ref,
		GED_VIEW_EDIT_GEOMETRY_TRANSIENT_PREVIEW, ged_preview_points,
		ged_preview_cmds, 3))
	FAIL("GED feature points replacement should update Obol edit-preview geometry");
    BObolFeatureRecord ged_preview_record;
	if (!feature_view_controller->features().record(ged_preview_handle,
		ged_preview_record) ||
	    ged_preview_record.kind != BObolFeatureKind::EditPreview ||
	    ged_preview_record.points.size() != 3 ||
	    ged_preview_record.commands.size() != 3)
	FAIL("GED feature points replacement should preserve Obol edit-preview records");
    struct directory *preview_box_dp = db_lookup(gedp->dbip, "box.s", LOOKUP_QUIET);
    struct rt_db_internal preview_box_intern;
    RT_DB_INTERNAL_INIT(&preview_box_intern);
    if (!preview_box_dp ||
	    rt_db_get_internal(&preview_box_intern, preview_box_dp, gedp->dbip,
		NULL) < 0)
	FAIL("GED feature primitive wireframe helper should load test primitive");
    mat_t preview_xform;
    MAT_IDN(preview_xform);
    MAT_DELTAS(preview_xform, 3.0, 0.0, 0.0);
    if (!ged_view_edit_primitive_wireframe_replace(feature_view_ctx,
		ged_preview_ref,
		gedp->dbip, &preview_box_intern, preview_xform, NULL, NULL)) {
	rt_db_free_internal(&preview_box_intern);
	FAIL("GED feature primitive wireframe helper should publish transformed primitive geometry");
    }
    rt_db_free_internal(&preview_box_intern);
	if (!feature_view_controller->features().record(ged_preview_handle,
		ged_preview_record) ||
	    ged_preview_record.kind != BObolFeatureKind::EditPreview ||
	    ged_preview_record.points.size() <= 3 ||
	    ged_preview_record.commands.size() != ged_preview_record.points.size())
	FAIL("GED feature primitive wireframe helper should replace edit-preview geometry");
    int preview_xform_seen = 0;
    for (size_t i = 0; i < ged_preview_record.points.size(); i++) {
	if (ged_preview_record.points[i][0] > 1.5f) {
	    preview_xform_seen = 1;
	    break;
	}
    }
    if (!preview_xform_seen)
	FAIL("GED feature primitive wireframe helper should apply the supplied transform");
    if (!ged_view_edit_preview_replace(feature_view_ctx,
		ged_preview_ref,
		"cap2::ged-preview-explicit.s",
		"edit::cap2::ged-preview",
		"move-handle",
		ged_preview_points,
		ged_preview_cmds,
		3,
		31,
		32) ||
	    !feature_view_controller->features().record(ged_preview_handle,
		ged_preview_record) ||
	    ged_preview_record.kind != BObolFeatureKind::EditPreview ||
	    ged_preview_record.identity != "cap2::ged-preview-explicit.s" ||
	    ged_preview_record.editIntentId != "edit::cap2::ged-preview" ||
	    ged_preview_record.editIntentRole != "move-handle" ||
	    ged_preview_record.sourceRevision != 31 ||
	    ged_preview_record.inputsRevision != 32)
	FAIL("GED feature edit preview replace should preserve explicit identity, intent, and revision metadata");
    if (!ged_view_edit_preview_replace(feature_view_ctx,
		ged_preview_ref,
		NULL,
		NULL,
		NULL,
		ged_preview_points,
		ged_preview_cmds,
		3,
		0,
		0) ||
	    !feature_view_controller->features().record(ged_preview_handle,
		ged_preview_record) ||
	    ged_preview_record.sourceRevision != 32 ||
	    ged_preview_record.inputsRevision != 33)
	FAIL("GED feature edit preview replace should advance preview revisions when explicit values are omitted");
    ged_view_edit_visible_set(feature_view_ctx, ged_preview_ref, 0);
    ged_view_edit_color_set(feature_view_ctx, ged_preview_ref, 12, 34, 56);
    BObolFeatureStyle ged_preview_style;
	if (!feature_view_controller->features().style(ged_preview_handle,
		ged_preview_style) ||
	    !ged_preview_style.hasVisible ||
	    ged_preview_style.visible ||
	    !ged_preview_style.hasColor)
	FAIL("GED feature style mutations should update Obol feature style");
    if (!ged_view_edit_preview_publish_event(feature_view_ctx,
		ged_preview_ref, GED_VIEW_EDIT_PREVIEW_UPDATE,
		"cap2::ged-source.s") ||
	    !feature_view_controller->features().exists("cap2::ged-preview"))
	FAIL("GED feature preview events should route through the Obol feature API");
    if (!ged_view_edit_preview_publish_event(feature_view_ctx,
		ged_preview_ref, GED_VIEW_EDIT_PREVIEW_COMMIT,
		"cap2::ged-source.s") ||
	    feature_view_controller->features().exists("cap2::ged-preview"))
	FAIL("GED feature commit preview events should retire transient Obol edit previews");
    ged_preview_ref =
	ged_view_edit_overlay_ensure(feature_view_ctx,
		"cap2::ged-preview", "cap2::ged-source.s");
	ged_preview_handle =
	    feature_view_controller->features().find("cap2::ged-preview");
    if (ged_view_edit_ref_is_null(ged_preview_ref) ||
	    !ged_preview_handle.isValid() ||
	    !ged_view_edit_points_replace(feature_view_ctx, ged_preview_ref,
		GED_VIEW_EDIT_GEOMETRY_TRANSIENT_PREVIEW, ged_preview_points,
		ged_preview_cmds, 3))
	FAIL("GED feature overlay ensure should recreate edit preview state after commit teardown");
    if (!ged_view_edit_geometry_clear(feature_view_ctx, ged_preview_ref) ||
	    !feature_view_controller->features().record(ged_preview_handle,
		ged_preview_record) ||
	    !ged_preview_record.points.empty())
	FAIL("GED feature clear geometry should clear the Obol feature record");

    ged_view_edit_ref ged_label_ref =
	ged_view_edit_label_ensure(feature_view_ctx, "cap2::ged-label");
    struct ged_view_feature_label ged_label;
    memset(&ged_label, 0, sizeof(ged_label));
    ged_label.text = "ged label";
    VSET(ged_label.point, 3.0, 4.0, 5.0);
    ged_label.color_valid = 1;
    ged_label.color[0] = 210;
    ged_label.color[1] = 211;
    ged_label.color[2] = 212;
    ged_label.font_size = 16.0;
    if (ged_view_edit_ref_is_null(ged_label_ref) ||
	    !ged_view_edit_labels_replace(feature_view_ctx, ged_label_ref,
		&ged_label, 1))
	FAIL("GED feature label replacement should route into Obol labels");
    BObolFeatureHandle ged_label_handle =
	feature_view_controller->features().find("cap2::ged-label");
    BObolFeatureRecord ged_label_record;
    if (!ged_label_handle.isValid() ||
	    !feature_view_controller->features().record(ged_label_handle,
		ged_label_record) ||
	    ged_label_record.kind != BObolFeatureKind::Labels ||
	    ged_label_record.labels.size() != 1 ||
	    !BU_STR_EQUAL(ged_label_record.labels[0].text.getString(),
		"ged label") ||
	    fabs(ged_label_record.labels[0].fontSize - 16.0f) > 0.001f ||
	    !ged_label_record.overlay.isOverlay ||
	    ged_label_record.overlay.ownerToken != feature_view_ctx)
	FAIL("GED feature label API should publish typed Obol label records");
    if (!ged_view_edit_remove_ref(feature_view_ctx, ged_label_ref) ||
	    feature_view_controller->features().exists("cap2::ged-label"))
	FAIL("GED feature reference removal should delete the exact Obol-backed feature record");
    ged_label_ref = ged_view_edit_label_ensure(feature_view_ctx,
	    "cap2::ged-label");
    if (ged_view_edit_ref_is_null(ged_label_ref) ||
	    !ged_view_feature_remove(feature_view_ctx, "cap2::ged-label") ||
	    feature_view_controller->features().exists("cap2::ged-label"))
	FAIL("GED feature name removal should delete Obol-backed feature records");

    if (ged_view_feature_remove_prefix(feature_view_ctx,
	    "cap2::") < 11 ||
	    owned_controller->features().exists("cap2::line") ||
	    feature_view_controller->features().exists("cap2::tcl-line") ||
	    feature_view_controller->features().exists("cap2::annotation") ||
	    feature_view_controller->features().exists("cap2::polygon-overlay") ||
	    owned_controller->features().exists("cap2::label") ||
	    feature_view_controller->features().exists("cap2::arrow") ||
	    owned_controller->features().exists("cap2::axes") ||
	    feature_view_controller->features().exists("cap2::tcl-axes") ||
	    owned_controller->features().exists("cap2::mesh") ||
	    owned_controller->features().exists("cap2::diagnostic") ||
	    feature_view_controller->features().exists("cap2::rt-preview"))
	FAIL("GED feature prefix removal should clear owned Obol feature store entries");
    if (!ged_view_context_obol_endpoint_set(feature_view_ctx, NULL, 0))
	FAIL("GED feature test should detach its temporary endpoint");
    feature_view_controller = NULL;

    size_t root_group_count = 77;
    if (!ged_draw_obol_group_descendant_group_count_for_path(gedp, "/",
	    &root_group_count) ||
	    root_group_count !=
		(size_t)owned_scene->getGroupDescendantGroupCount("/") ||
	    ged_draw_has_groups(gedp) !=
		(root_group_count > 0 ? 1 : 0))
	FAIL("GED public group presence should match owned Obol scene groups");
    const size_t original_root_group_count = root_group_count;
    if (!ged_draw_obol_group_ensure_for_path(gedp,
	    "__obol_root_group_presence.s",
	    "__obol_root_group_presence.s",
	    GED_DRAW_MODE_WIRE,
	    0))
	FAIL("GED Obol root group-presence sentinel should be ensured");
    root_group_count = 0;
    if (!ged_draw_obol_group_descendant_group_count_for_path(gedp, "/",
	    &root_group_count) ||
	    root_group_count != original_root_group_count + 1 ||
	    owned_scene->getGroupDescendantGroupCount("/") !=
		(int)(original_root_group_count + 1) ||
	    !ged_draw_has_groups(gedp))
	FAIL("GED public group presence should prefer owned Obol scene groups");
    group_source_state root_group_state = {0, GED_DRAW_GROUP_REF_NULL, NULL,
	"__obol_root_group_presence.s"};
    ged_draw_foreach_group_record(gedp, group_source_state_cb,
	    &root_group_state);
    if (!root_group_state.found ||
	    ged_draw_group_ref_is_null(root_group_state.ref))
	FAIL("GED source-root group traversal should enumerate owned Obol groups");
    group_index_state root_group_index = {0, GED_DRAW_GROUP_REF_NULL_INIT};
    if (ged_draw_group_ref_index_for_component(gedp,
	    "__obol_root_group_presence.s",
	    group_index_state_cb, &root_group_index) != 1 ||
	    root_group_index.count != 1 ||
	    ged_draw_group_ref_is_null(root_group_index.ref))
	FAIL("GED Obol group component index should enumerate owned Obol groups");
    if (!ged_draw_group_ref_set_mode(gedp, root_group_state.ref,
	    GED_DRAW_MODE_SHADED))
	FAIL("GED source-root group traversal should return a mutable Obol group ref");
    SoGroup *root_presence_group =
	owned_scene->findGroup("__obol_root_group_presence.s");
    if (!root_presence_group ||
	    !root_presence_group->isOfType(SoBRLSceneGroup::getClassTypeId()) ||
	    static_cast<SoBRLSceneGroup *>(root_presence_group)->
		drawMode.getValue() != BOBOL_LOD_DRAW_SHADED)
	FAIL("GED source-root group traversal ref should mutate the owned Obol group");
    if (owned_scene->removeGroup("__obol_root_group_presence.s") <= 0)
	FAIL("GED Obol root group-presence sentinel should be removable");
    root_group_count = 77;
    if (!ged_draw_obol_group_descendant_group_count_for_path(gedp, "/",
	    &root_group_count) ||
	    root_group_count != original_root_group_count ||
	    owned_scene->getGroupDescendantGroupCount("/") !=
		(int)original_root_group_count ||
	    ged_draw_has_groups(gedp) !=
		(original_root_group_count > 0 ? 1 : 0))
	FAIL("GED public group presence should clear with owned Obol scene groups");
    if (!ged_draw_obol_database_source_ensure_for_path(gedp,
	    "group_only.s", gedp->dbip, GED_DRAW_MODE_WIRE, 1002))
	FAIL("GED Obol root shape iterator sentinel should insert an owned source");
    record_source_state root_shape_record = {0, GED_DRAW_SHAPE_REF_NULL_INIT,
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "group_only.s", 0, 0, 0, 0, 0.0,
	0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &root_shape_record);
    if (!root_shape_record.found ||
	    root_shape_record.sourceRevision != 1002 ||
	    root_shape_record.drawMode != GED_DRAW_MODE_WIRE)
	FAIL("GED source-root shape traversal should enumerate owned Obol database sources");
    shape_index_state root_shape_index = {0, GED_DRAW_SHAPE_REF_NULL_INIT};
    if (ged_draw_shape_ref_index_for_component(gedp, "group_only.s",
	    shape_index_state_cb, &root_shape_index) != 1 ||
	    root_shape_index.count != 1 ||
	    ged_draw_shape_ref_is_null(root_shape_index.ref))
	FAIL("GED Obol shape component index should enumerate owned Obol database sources");
    shape_index_state root_shape_path_hash_index = {
	0, GED_DRAW_SHAPE_REF_NULL_INIT
    };
    ged_draw_index_stats_reset(gedp);
    if (!root_shape_record.pathHash ||
	    ged_draw_shape_ref_index_for_path_hash(gedp,
		root_shape_record.pathHash, shape_index_state_cb,
		&root_shape_path_hash_index) != 1 ||
	    root_shape_path_hash_index.count != 1 ||
	    ged_draw_shape_ref_is_null(root_shape_path_hash_index.ref))
	FAIL("GED Obol shape path-hash index should enumerate owned Obol database sources");
    struct ged_draw_index_stats root_shape_path_hash_stats;
    memset(&root_shape_path_hash_stats, 0,
	    sizeof(root_shape_path_hash_stats));
    ged_draw_index_stats_get(gedp, &root_shape_path_hash_stats);
    if (root_shape_path_hash_stats.path_queries ||
	    root_shape_path_hash_stats.path_candidates)
	FAIL("GED Obol shape path-hash index should avoid registry path-index queries");
    if (ged_draw_group_ref_is_null(root_shape_record.group) ||
	    !ged_draw_group_ref_set_mode(gedp, root_shape_record.group,
		GED_DRAW_MODE_SHADED))
	FAIL("GED source-root shape traversal should expose an owned Obol group ref");
    SoGroup *root_shape_group = owned_scene->findGroup("group_only.s");
    if (!root_shape_group ||
	    !root_shape_group->isOfType(SoBRLSceneGroup::getClassTypeId()) ||
	    static_cast<SoBRLSceneGroup *>(root_shape_group)->
		drawMode.getValue() != BOBOL_LOD_DRAW_SHADED)
	FAIL("GED source-root shape group ref should mutate the owned Obol group");
    if (!ged_draw_group_ref_set_mode(gedp, root_shape_record.group,
	    GED_DRAW_MODE_WIRE))
	FAIL("GED source-root shape group mode restore should succeed");
    if (owned_scene->removeDatabaseSource("group_only.s") <= 0)
	FAIL("GED Obol root shape iterator sentinel should be removable");
    if (owned_scene->removeGroup("group_only.s") <= 0)
	FAIL("GED Obol root shape iterator sentinel group should be removable");

    if (exercise_mode_specific_source_lifecycle(gedp, owned_scene,
	    "box.s", GED_DRAW_MODE_EVAL_WIRE,
	    SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE, 0, 0, 1,
	    "evaluated-wire"))
	return 1;
    if (exercise_evaluated_wire_shape_ref_realize_context(gedp, owned_scene,
	    "box.s"))
	return 1;
    if (exercise_mode_specific_source_lifecycle(gedp, owned_scene,
	    "box.s", GED_DRAW_MODE_HIDDEN_LINE,
	    SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE, 0, 1, 1,
	    "hidden-line"))
	return 1;
    if (exercise_mode_specific_source_lifecycle(gedp, owned_scene,
	    "box.s", GED_DRAW_MODE_EVAL_POINTS,
	    SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS, 0, 0, 1,
	    "evaluated-points"))
	return 1;
    if (exercise_deferred_mode_replacement(gedp, owned_scene, "box.s"))
	return 1;
    const char *draw_eval_points_option[4] = {"draw",
	"--evaluated-points", "--add-mode", "box.s"};
    if (ged_exec_draw(gedp, 4, draw_eval_points_option) != BRLCAD_OK)
	FAIL("GED evaluated-points option draw should succeed");
    if (source_representation_count(owned_scene, "box.s",
	    SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS) != 1)
	FAIL("GED evaluated-points option should create one mode source");
    (void)try_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, "box.s",
	    ged_draw_active_view_ctx(gedp), -1);
    if (!ged_draw_obol_database_source_ensure_for_path(gedp, "box.s",
	    gedp->dbip, GED_DRAW_MODE_WIRE, 0))
	FAIL("GED evaluated-points option test should restore shared wire source");
    const char *draw_annot_line[2] = {"draw", "annot_line.s"};
    if (ged_exec_draw(gedp, 2, draw_annot_line) != BRLCAD_OK)
	FAIL("GED annotation draw should succeed for owned Obol publication");
    record_source_state annot_record = {0, GED_DRAW_SHAPE_REF_NULL_INIT,
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "annot_line.s", 0, 0, 0, 0, 0.0,
	0};
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
    BObolRealizedShapeSummary annot_summary;
    SoBRLExportAction annot_export;
    if (annot_source)
	annot_export.apply(annot_source);
    if (owned_scene->getDatabaseSourceCount() != 3 ||
	    !annot_source || annot_source->getRealizedShapeCount() != 0 ||
	    !annot_source->getRealizedShapeSummary(0, annot_summary) ||
	    annot_summary.segmentCount != 1 ||
	    bu_strcmp(annot_summary.sourceType.getString(), "annotation") != 0 ||
	    bu_strcmp(annot_summary.geometryKind.getString(), "annotation") != 0 ||
	    annot_export.getLineCount() != 1 ||
	    fabs(annot_export.getLine(0).b[0] - 50.25f) > 0.001f ||
	    fabs(annot_export.getLine(0).b[1] - 0.5f) > 0.001f ||
	    fabs(annot_export.getLine(0).b[2]) > 0.001f)
	FAIL("GED annotation draw should publish line segments into the owned Obol source");
    const char *erase_annot_line[2] = {"erase", "annot_line.s"};
    if (ged_exec_erase(gedp, 2, erase_annot_line) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "annot_line.s"))
	FAIL("GED annotation erase should restore the owned Obol source baseline");
    const char *draw_submodel_owner[2] = {"draw", "submodel_owner.s"};
    if (ged_exec_draw(gedp, 2, draw_submodel_owner) != BRLCAD_OK)
	FAIL("GED submodel draw should succeed for owned Obol direct publication");
    SoBRLDatabaseSource *submodel_source =
	source_for_path(owned_scene, "submodel_owner.s");
    if (!submodel_source || owned_scene->getDatabaseSourceCount() != 3)
	FAIL("GED submodel draw should create an owned Obol source");
    BObolRealizedShapeSummary submodel_summary;
    if (submodel_source->getRealizedShapeCount() != 0 ||
	    !submodel_source->getRealizedShapeSummary(0, submodel_summary) ||
	    submodel_summary.pointCount == 0 ||
	    submodel_summary.commandCount != submodel_summary.pointCount ||
	    bu_strcmp(submodel_summary.recordRole.getString(), "database") != 0 ||
	    bu_strcmp(submodel_summary.sourceType.getString(), "submodel") != 0 ||
	    auxiliary_for_path_variant(submodel_source, "box.s")) {
	FAIL("GED submodel draw should realize direct primary owned Obol geometry without legacy auxiliary staging");
    }
    const char *erase_submodel_owner[2] = {"erase", "submodel_owner.s"};
    if (ged_exec_erase(gedp, 2, erase_submodel_owner) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "submodel_owner.s"))
	FAIL("GED submodel erase should restore the owned Obol source baseline");
    const char *draw_submodel_temp_owner[2] = {
	"draw", "submodel_temp_owner.s"
    };
    if (ged_exec_draw(gedp, 2, draw_submodel_temp_owner) != BRLCAD_OK)
	FAIL("GED submodel temp-source draw should succeed for owned Obol direct publication");
    SoBRLDatabaseSource *submodel_temp_source =
	source_for_path(owned_scene, "submodel_temp_owner.s");
    if (!submodel_temp_source || owned_scene->getDatabaseSourceCount() != 3 ||
	    source_for_path(owned_scene, "nested_leaf.s"))
	FAIL("GED submodel temp-source draw should not leak a temporary owned Obol leaf source");
    BObolRealizedShapeSummary submodel_temp_summary;
    if (submodel_temp_source->getRealizedShapeCount() != 0 ||
	    !submodel_temp_source->getRealizedShapeSummary(0,
		submodel_temp_summary) ||
	    submodel_temp_summary.pointCount == 0 ||
	    bu_strcmp(submodel_temp_summary.recordRole.getString(), "database") != 0 ||
	    auxiliary_for_path_variant(submodel_temp_source, "nested_leaf.s"))
	FAIL("GED submodel temp-source draw should realize direct primary owned Obol geometry");
    const char *erase_submodel_temp_owner[2] = {
	"erase", "submodel_temp_owner.s"
    };
    if (ged_exec_erase(gedp, 2, erase_submodel_temp_owner) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "submodel_temp_owner.s") ||
	    source_for_path(owned_scene, "nested_leaf.s"))
	FAIL("GED submodel temp-source erase should restore the owned Obol source baseline");
    const char *lod_mesh_on[4] = {"view", "lod", "mesh", "1"};
    if (ged_exec_view(gedp, 4, lod_mesh_on) != BRLCAD_OK)
	FAIL("GED view LoD mesh enable should succeed for owned BoT update test");
    const char *lod_bot_threshold_zero[4] = {
	"view", "lod", "bot_threshold", "0"
    };
    if (ged_exec_view(gedp, 4, lod_bot_threshold_zero) != BRLCAD_OK)
	FAIL("GED view LoD BoT threshold should be settable for owned BoT update test");
    const char *lod_mesh_cache[3] = {"view", "lod", "cache"};
    if (ged_exec_view(gedp, 3, lod_mesh_cache) != BRLCAD_OK)
	FAIL("GED view LoD cache should be buildable for owned BoT update test");
    const char *draw_mesh_owner[2] = {"draw", "mesh_owner.bot"};
    if (ged_exec_draw(gedp, 2, draw_mesh_owner) != BRLCAD_OK)
	FAIL("GED BoT mesh LoD draw should succeed for owned Obol mesh update");
    record_source_state mesh_record = {0, GED_DRAW_SHAPE_REF_NULL_INIT,
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "mesh_owner.bot", 0, 0, 0, 0,
	0.0, 0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &mesh_record);
    if (!mesh_record.found)
	FAIL("GED BoT mesh LoD draw should create a shape record");
    struct ged_view_context *mesh_lod_view_ctx = ged_draw_shape_ref_view_context(gedp,
	    mesh_record.ref);
    if (!mesh_lod_view_ctx)
	FAIL("GED BoT mesh LoD view context should be available");
    struct ged_view_context *mesh_lod_view_ctxs[1] = {mesh_lod_view_ctx};
    if (!ged_draw_shape_ref_lod_ensure(gedp, mesh_record.ref,
	    mesh_lod_view_ctx, mesh_lod_view_ctxs, 1))
	FAIL("GED BoT mesh LoD ensure should succeed for owned Obol mesh update");
    BObolViewController *mesh_view_controller =
	ged_bobol_view_controller(mesh_lod_view_ctx);
    SoBRLDatabaseSource *mesh_source = NULL;
    SoBRLMeshShape *mesh_shape = NULL;
    BObolProgressiveOptions mesh_progressive_options;
    BObolProgressiveStatus mesh_progressive_status;
    for (int attempt = 0; attempt < 2000; ++attempt) {
	mesh_source = source_for_path(owned_scene, "mesh_owner.bot");
	mesh_shape = mesh_source ? mesh_source->getRealizedMesh() : NULL;
	if (mesh_shape && mesh_shape->point.getNum() > 0 &&
	    mesh_shape->coordIndex.getNum() > 0)
	    break;
	if (mesh_view_controller)
	    (void)mesh_view_controller->advanceProgressiveWork(
		&mesh_progressive_options, &mesh_progressive_status);
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    SbVec3f mesh_lod_bmin;
    SbVec3f mesh_lod_bmax;
    if (!mesh_source || !mesh_shape ||
	    mesh_shape->point.getNum() == 0 ||
	    mesh_shape->coordIndex.getNum() == 0) {
	FAIL("GED BoT mesh LoD update should publish owned Obol mesh fields");
    }
    if (!mesh_shape->isLodBackedMesh())
	FAIL("GED BoT mesh LoD update should realize a LoD-backed Obol mesh");
    if (!mesh_source->getMeshLod())
	FAIL("GED Obol mesh LoD draw should store the runtime handle on the owned source");
    if (!mesh_source->getMeshLodBounds(mesh_lod_bmin, mesh_lod_bmax))
	FAIL("GED Obol mesh LoD draw should store runtime bounds on the owned source");

    struct BObolMeshLod *retained_mesh_lod = mesh_source->getMeshLod();
    mesh_source->markStale(SoBRLDatabaseSource::STALE_VIEW);
    SbVec3f retained_mesh_lod_bmin;
    SbVec3f retained_mesh_lod_bmax;
    if (mesh_source->getMeshLod() != retained_mesh_lod ||
	!mesh_source->getMeshLodBounds(retained_mesh_lod_bmin,
	    retained_mesh_lod_bmax) ||
	!retained_mesh_lod_bmin.equals(mesh_lod_bmin, 0.0f) ||
	!retained_mesh_lod_bmax.equals(mesh_lod_bmax, 0.0f))
	FAIL("view-only invalidation should retain the source-owned mesh LoD record");
    if (!ged_draw_shape_ref_lod_ensure(gedp, mesh_record.ref,
	    mesh_lod_view_ctx, mesh_lod_view_ctxs, 1) ||
	mesh_source->getMeshLod() != retained_mesh_lod)
	FAIL("GED BoT mesh LoD ensure should reuse the retained view-independent reader");

    const BObolSourceRealizationStamp valid_mesh_lod_stamp =
	mesh_source->captureRealizationStamp();
    const uint64_t mesh_lod_frame = owned_scene->getFrameRevision();
    struct BObolMeshLod *invalid_bounds_candidate =
	bobol_mesh_lod_clone_reader(retained_mesh_lod);
    SbVec3f invalid_mesh_lod_bmin = mesh_lod_bmin;
    invalid_mesh_lod_bmin[0] = std::numeric_limits<float>::quiet_NaN();
    if (!invalid_bounds_candidate ||
	owned_scene->adoptDatabaseSourceInstanceMeshLod(
	    valid_mesh_lod_stamp.sourceInstanceKey().getString(),
	    valid_mesh_lod_stamp, invalid_bounds_candidate,
	    invalid_mesh_lod_bmin, mesh_lod_bmax) ||
	mesh_source->getMeshLod() != retained_mesh_lod ||
	owned_scene->getFrameRevision() != mesh_lod_frame)
	FAIL("mesh LoD adoption should reject invalid bounds without changing ownership");
    bobol_mesh_lod_destroy(invalid_bounds_candidate);

    struct BObolMeshLod *replacement_mesh_lod =
	bobol_mesh_lod_clone_reader(retained_mesh_lod);
    SbVec3f replacement_mesh_lod_bmin;
    SbVec3f replacement_mesh_lod_bmax;
    if (!replacement_mesh_lod ||
	!owned_scene->adoptDatabaseSourceInstanceMeshLod(
	    valid_mesh_lod_stamp.sourceInstanceKey().getString(),
	    valid_mesh_lod_stamp, replacement_mesh_lod,
	    mesh_lod_bmin, mesh_lod_bmax) ||
	mesh_source->getMeshLod() != replacement_mesh_lod ||
	!mesh_source->getMeshLodBounds(replacement_mesh_lod_bmin,
	    replacement_mesh_lod_bmax) ||
	!replacement_mesh_lod_bmin.equals(mesh_lod_bmin, 0.0f) ||
	!replacement_mesh_lod_bmax.equals(mesh_lod_bmax, 0.0f) ||
	owned_scene->getFrameRevision() != mesh_lod_frame)
	FAIL("mesh LoD adoption should transfer one complete handle-and-bounds record");

    struct BObolMeshLod *stale_mesh_lod_candidate =
	bobol_mesh_lod_clone_reader(replacement_mesh_lod);
    const BObolSourceRealizationStamp stale_mesh_lod_stamp =
	mesh_source->captureRealizationStamp();
    BObolDatabaseSourceSummary mesh_source_summary;
    if (!stale_mesh_lod_candidate ||
	!mesh_source->getSummary(mesh_source_summary) ||
	!mesh_source_summary.valid ||
	owned_scene->setDatabaseSourceInstanceState(
	    mesh_source_summary.instanceKey.getString(), FALSE,
	    mesh_source_summary.sourceRevision,
	    mesh_source_summary.inputsRevision + 1,
	    mesh_source_summary.visible, mesh_source_summary.selected,
	    mesh_source_summary.highlighted, mesh_source_summary.lineStyle,
	    mesh_source_summary.lineWidth, mesh_source_summary.transparency,
	    mesh_source_summary.colorOverride, mesh_source_summary.color,
	    mesh_source_summary.materialColorValid,
	    mesh_source_summary.materialColor,
	    mesh_source_summary.materialRevision) <= 0)
	FAIL("mesh LoD source-input invalidation setup should succeed");
    if (mesh_source->getMeshLod() ||
	mesh_source->getMeshLodBounds(retained_mesh_lod_bmin,
	    retained_mesh_lod_bmax))
	FAIL("source-input invalidation should retire the mesh LoD record");
    if (owned_scene->adoptDatabaseSourceInstanceMeshLod(
	    stale_mesh_lod_stamp.sourceInstanceKey().getString(),
	    stale_mesh_lod_stamp, stale_mesh_lod_candidate,
	    mesh_lod_bmin, mesh_lod_bmax))
	FAIL("mesh LoD adoption should reject a result from stale source inputs");
    bobol_mesh_lod_destroy(stale_mesh_lod_candidate);
    if (!ged_draw_shape_ref_lod_ensure(gedp, mesh_record.ref,
	    mesh_lod_view_ctx, mesh_lod_view_ctxs, 1) ||
	!mesh_source->getMeshLod() ||
	!mesh_source->getMeshLodBounds(retained_mesh_lod_bmin,
	    retained_mesh_lod_bmax))
	FAIL("GED BoT mesh LoD ensure should rebuild a retired source-resource record");

    const char *erase_mesh_owner[2] = {"erase", "mesh_owner.bot"};
    if (ged_exec_erase(gedp, 2, erase_mesh_owner) != BRLCAD_OK ||
	owned_scene->getDatabaseSourceCount() != 2 ||
	source_for_path(owned_scene, "mesh_owner.bot"))
	FAIL("GED BoT mesh LoD exact-identity setup should restore the source baseline");

    const char *draw_mesh_owner_shaded[3] = {
	"draw", "-m1", "mesh_owner.bot"
    };
    if (ged_exec_draw(gedp, 3, draw_mesh_owner_shaded) != BRLCAD_OK)
	FAIL("GED BoT shaded sibling draw should succeed for exact mesh LoD ownership");
    const char *draw_mesh_owner_wire[4] = {
	"draw", "-m0", "--add-mode", "mesh_owner.bot"
    };
    if (ged_exec_draw(gedp, 4, draw_mesh_owner_wire) != BRLCAD_OK)
	FAIL("GED BoT wire sibling draw should succeed for exact mesh LoD ownership");
    record_source_mode_state shaded_mesh_record = {};
    shaded_mesh_record.recordState.ref = GED_DRAW_SHAPE_REF_NULL;
    shaded_mesh_record.recordState.group = GED_DRAW_GROUP_REF_NULL;
    shaded_mesh_record.recordState.matchPath = "mesh_owner.bot";
    shaded_mesh_record.matchDrawMode = GED_DRAW_MODE_SHADED_BOTS;
    ged_draw_foreach_shape_record(gedp, record_source_mode_state_cb,
	    &shaded_mesh_record);
    record_source_mode_state wire_mesh_record = {};
    wire_mesh_record.recordState.ref = GED_DRAW_SHAPE_REF_NULL;
    wire_mesh_record.recordState.group = GED_DRAW_GROUP_REF_NULL;
    wire_mesh_record.recordState.matchPath = "mesh_owner.bot";
    wire_mesh_record.matchDrawMode = GED_DRAW_MODE_WIRE;
    ged_draw_foreach_shape_record(gedp, record_source_mode_state_cb,
	    &wire_mesh_record);
    SoBRLDatabaseSource *shaded_mesh_source = source_for_representation(
	owned_scene, "mesh_owner.bot", GED_DRAW_MODE_SHADED_BOTS);
    SoBRLDatabaseSource *wire_mesh_source = source_for_representation(
	owned_scene, "mesh_owner.bot", GED_DRAW_MODE_WIRE);
    if (!shaded_mesh_record.recordState.found || !shaded_mesh_source)
	FAIL("GED BoT shaded mesh LoD source should have an exact shape record");
    if (!wire_mesh_record.recordState.found || !wire_mesh_source)
	FAIL("GED BoT wire mesh LoD source should have an exact shape record");
    if (shaded_mesh_source == wire_mesh_source)
	FAIL("GED BoT mesh LoD modes should retain distinct source identities");
    if (!ged_draw_shape_ref_lod_ensure(gedp,
	    shaded_mesh_record.recordState.ref, mesh_lod_view_ctx,
	    mesh_lod_view_ctxs, 1) || !shaded_mesh_source->getMeshLod())
	FAIL("GED BoT shaded sibling should install its own mesh LoD record");
    if (!ged_draw_shape_ref_lod_ensure(gedp,
	    wire_mesh_record.recordState.ref, mesh_lod_view_ctx,
	    mesh_lod_view_ctxs, 1) || !wire_mesh_source->getMeshLod())
	FAIL("GED BoT wire sibling should install its own mesh LoD record");

    struct BObolMeshLod *wire_mesh_lod = wire_mesh_source->getMeshLod();
    BObolDatabaseSourceSummary shaded_mesh_summary;
    if (!shaded_mesh_source->getSummary(shaded_mesh_summary) ||
	!shaded_mesh_summary.valid ||
	owned_scene->setDatabaseSourceInstanceState(
	    shaded_mesh_summary.instanceKey.getString(), FALSE,
	    shaded_mesh_summary.sourceRevision,
	    shaded_mesh_summary.inputsRevision + 1,
	    shaded_mesh_summary.visible, shaded_mesh_summary.selected,
	    shaded_mesh_summary.highlighted, shaded_mesh_summary.lineStyle,
	    shaded_mesh_summary.lineWidth, shaded_mesh_summary.transparency,
	    shaded_mesh_summary.colorOverride, shaded_mesh_summary.color,
	    shaded_mesh_summary.materialColorValid,
	    shaded_mesh_summary.materialColor,
	    shaded_mesh_summary.materialRevision) <= 0)
	FAIL("GED BoT shaded mesh LoD invalidation setup should succeed");
    if (shaded_mesh_source->getMeshLod() ||
	wire_mesh_source->getMeshLod() != wire_mesh_lod)
	FAIL("mode-scoped invalidation should retire only the exact mesh LoD record");
    if (!ged_draw_shape_ref_lod_ensure(gedp,
	    shaded_mesh_record.recordState.ref, mesh_lod_view_ctx,
	    mesh_lod_view_ctxs, 1) || !shaded_mesh_source->getMeshLod() ||
	wire_mesh_source->getMeshLod() != wire_mesh_lod)
	FAIL("mode-scoped mesh LoD ensure should rebuild only the exact source");
    std::puts("PASS GED mesh-LoD attachment: exact ownership, atomic bounds, stale rejection and invalidation");

    if (ged_exec_erase(gedp, 2, erase_mesh_owner) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "mesh_owner.bot"))
	FAIL("GED BoT mesh LoD erase should restore the owned Obol source baseline");
    const char *draw_brep_wire[2] = {"draw", "brep_owner.brep"};
    if (ged_exec_draw(gedp, 2, draw_brep_wire) != BRLCAD_OK)
	FAIL("GED BREP wireframe draw should succeed for owned Obol line-set update");
    record_source_state brep_wire_record = {0, GED_DRAW_SHAPE_REF_NULL_INIT,
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "brep_owner.brep", 0, 0, 0,
	0, 0.0, 0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &brep_wire_record);
    if (!brep_wire_record.found)
	FAIL("GED BREP wireframe draw should create a shape record");
    struct directory *brep_wire_dp = db_lookup(gedp->dbip,
	    "brep_owner.brep", LOOKUP_QUIET);
    if (!brep_wire_dp)
	FAIL("brep_owner.brep should be available for BREP wireframe publication");
    struct rt_db_internal brep_wire_intern;
    RT_DB_INTERNAL_INIT(&brep_wire_intern);
    mat_t brep_wire_identity;
    MAT_IDN(brep_wire_identity);
    if (rt_db_get_internal(&brep_wire_intern, brep_wire_dp, gedp->dbip,
	    brep_wire_identity) < 0)
	FAIL("brep_owner.brep internal lookup should succeed");
    struct bn_tol brep_wire_tol;
    BN_TOL_INIT_SET_TOL(&brep_wire_tol);
    int brep_wire_publish =
	ged_draw_obol_database_source_publish_primitive_wireframe_for_path(
		gedp, "brep_owner.brep", &brep_wire_intern, NULL,
		&brep_wire_tol);
    rt_db_free_internal(&brep_wire_intern);
    if (brep_wire_publish < 0)
	FAIL("GED BREP wireframe publication should succeed for owned Obol line-set update");
    SoBRLDatabaseSource *brep_wire_source =
	source_for_path(owned_scene, "brep_owner.brep");
    SoBRLVListShape *brep_wire_shape = brep_wire_source ?
	brep_wire_source->getRealizedShape() : NULL;
    if (!brep_wire_source || !brep_wire_shape ||
	    brep_wire_shape->point.getNum() == 0 ||
	    brep_wire_shape->command.getNum() == 0 ||
	    bu_strcmp(brep_wire_shape->sourceType.getValue().getString(),
		"line-set") != 0)
	FAIL("GED BREP wireframe draw should publish owned Obol line geometry");
    const char *erase_brep_wire[2] = {"erase", "brep_owner.brep"};
    if (ged_exec_erase(gedp, 2, erase_brep_wire) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "brep_owner.brep"))
	FAIL("GED BREP wireframe erase should restore the owned Obol source baseline");
    const char *lod_brep_threshold_one[4] = {
	"view", "lod", "bot_threshold", "1"
    };
    if (ged_exec_view(gedp, 4, lod_brep_threshold_one) != BRLCAD_OK)
	FAIL("GED view LoD threshold should enable adaptive BRep realization");
    struct ged_view_context *brep_draw_view_ctx = ged_view_active_ctx(gedp);
    if (!brep_draw_view_ctx ||
	    !ged_view_context_display_endpoint_ensure(brep_draw_view_ctx))
	FAIL("GED BREP mesh LoD draw should attach a progressive display endpoint");
    const char *draw_brep_owner[3] = {"draw", "-m1", "brep_owner.brep"};
    if (ged_exec_draw(gedp, 3, draw_brep_owner) != BRLCAD_OK)
	FAIL("GED BREP mesh LoD draw should succeed for owned Obol mesh update");
    record_source_state brep_record = {0, GED_DRAW_SHAPE_REF_NULL_INIT,
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "brep_owner.brep", 0, 0, 0, 0,
	0.0, 0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &brep_record);
    if (!brep_record.found)
	FAIL("GED BREP mesh LoD draw should create a shape record");
    struct ged_view_context *brep_lod_view_ctx = ged_draw_shape_ref_view_context(gedp,
	    brep_record.ref);
    if (!brep_lod_view_ctx)
	FAIL("GED BREP mesh LoD view context should be available");
    struct ged_view_context *brep_lod_view_ctxs[1] = {brep_lod_view_ctx};
    if (!ged_draw_shape_ref_lod_ensure(gedp, brep_record.ref,
	    brep_lod_view_ctx, brep_lod_view_ctxs, 1))
	FAIL("GED BREP mesh LoD ensure should succeed for owned Obol mesh update");

    /* Unlike the database BoT path, an adaptive BRep tessellation is a
     * background provider job.  Pump until its first useful retained mesh is
     * published before inspecting fields.  Full presentation-gated terminal
     * refinement is covered by the shared-library NIST image matrix; this
     * white-box test has no render host and must not pretend that it can
     * certify a presented frame. */
    BObolViewController *brep_view_controller =
	ged_bobol_view_controller(brep_lod_view_ctx);
    bv_dimensions_set(DRAW_TEST_BV(brep_lod_view_ctx), 512, 512);
    const char *brep_autoview[1] = {"autoview"};
    if (!brep_view_controller ||
	    ged_exec_autoview(gedp, 1, brep_autoview) != BRLCAD_OK ||
	    !brep_view_controller->syncCameraFromViewContext(
		brep_lod_view_ctx))
	FAIL("GED BREP mesh LoD update should establish a visible test camera");
    SoBRLDatabaseSource *brep_source = NULL;
    const BObolViewLodState::CadPayload *brep_payload = NULL;
    BObolProgressiveOptions brep_progressive_options;
    brep_progressive_options.forceTerminalLodRefinement = TRUE;
    BObolProgressiveStatus brep_progressive_status;
    for (int attempt = 0; brep_view_controller && attempt < 10000;
	    attempt++) {
	brep_source = source_for_path(owned_scene, "brep_owner.brep");
	BObolViewLodState *brep_view_state =
	    brep_view_controller->getViewLodState();
	std::vector<const BObolViewLodState::CadPayload *> brep_payloads;
	if (brep_source && brep_view_state)
	    brep_view_state->findCadPayloadsUnordered(brep_source,
		brep_payloads);
	brep_payload = brep_payloads.size() == 1 ? brep_payloads[0] : NULL;
	if (brep_payload && brep_payload->isValid() &&
		brep_payload->progressiveMesh &&
		brep_payload->counts.faceCount > 0)
	    break;
	(void)brep_view_controller->advanceProgressiveWork(
	    &brep_progressive_options, &brep_progressive_status);
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!brep_source || !brep_payload || !brep_payload->isValid() ||
	    !brep_payload->progressiveMesh ||
	    !brep_payload->progressiveMesh->isValid() ||
	    !brep_payload->preparedCadGeometry ||
	    brep_payload->counts.faceCount == 0 ||
	    brep_payload->counts.pointCount == 0 ||
	    brep_payload->activeCut <
		brep_payload->progressiveMesh->minimumCut() ||
	    brep_payload->bounds.isEmpty())
	FAIL("GED BREP mesh LoD update should publish a view-local retained PoP payload");
    if (brep_source->getRealizedMesh() || brep_source->getMeshLod())
	FAIL("GED compact BREP LoD should not copy view-local mesh state onto the shared source");
    const char *erase_brep_owner[2] = {"erase", "brep_owner.brep"};
    if (ged_exec_erase(gedp, 2, erase_brep_owner) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "brep_owner.brep"))
	FAIL("GED BREP mesh LoD erase should restore the owned Obol source baseline");
    if (!ged_view_context_obol_endpoint_set(brep_lod_view_ctx, NULL, 0))
	FAIL("GED BREP mesh LoD test should restore the detached lifecycle fixture");
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
    BObolDatabaseSourceSummary renamed_summary;
    if (source_for_path(owned_scene, "rename_source.s") ||
	    !renamed_source ||
	    owned_scene->getDatabaseSourceCount() != 3 ||
	    !renamed_source->getSummary(renamed_summary) ||
	    bu_strcmp(renamed_summary.path.getString(), "renamed_source.s") != 0 ||
	    renamed_summary.lineWidth != 11 ||
	    renamed_summary.inputsRevision != 7 ||
	    renamed_summary.sourceRevision == 9191)
	FAIL("GED rename transaction should rename the owned Obol source in place");
    BObolRealizedShapeSummary renamed_shape_summary;
    if (renamed_source->getRealizedShapeCount() != 0 ||
	    !renamed_source->getRealizedShapeSummary(0,
		renamed_shape_summary) ||
	    !path_equal(renamed_shape_summary.path.getString(),
		"renamed_source.s") ||
	    !path_equal(renamed_shape_summary.ownerSourcePath.getString(),
		"renamed_source.s"))
	FAIL("GED rename transaction should retarget owned Obol realized shape metadata");
    const char *erase_renamed_source[2] = {"erase", "renamed_source.s"};
    if (ged_exec_erase(gedp, 2, erase_renamed_source) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "renamed_source.s"))
	FAIL("GED renamed source erase should restore the owned Obol source baseline");
    if (owned_scene->findGroup("rename_source.s") ||
	    owned_scene->findGroup("renamed_source.s"))
	FAIL("GED source rename/erase should retarget and prune owned Obol source-owner groups");
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
    ged_draw_scene_handle group_only_scene_handle =
	ged_draw_registry_group_ref_scene_handle(gedp, group_only_ref);
    if (ged_draw_scene_handle_backend(group_only_scene_handle) !=
	    GED_DRAW_SCENE_BACKEND_OBOL)
	FAIL("GED group lookup/create should return an owned Obol group ref");
    if (!owned_scene->findGroup("group_only.s") ||
	    owned_scene->getGroupDatabaseSourceCount("group_only.s") != 0)
	FAIL("GED group lookup/create should ensure an empty owned Obol group");
    struct ged_draw_group_record group_only_record;
    memset(&group_only_record, 0, sizeof(group_only_record));
    if (!ged_draw_group_record_get(gedp, group_only_ref,
	    &group_only_record) ||
	    !path_equal(group_only_record.path, "group_only.s") ||
	    group_only_record.shape_count !=
		owned_scene->getGroupDatabaseSourceCount("group_only.s"))
	FAIL("GED group lookup/create should keep public records aligned with owned Obol empty groups");
    if (!path_equal(ged_draw_registry_group_ref_semantic_path(gedp,
	    group_only_ref), "group_only.s"))
	FAIL("GED group refs should retain a semantic registry path for Obol lookup");
    if (!owned_scene->ensureGroup("group_only.s/obol_child.s"))
	FAIL("owned Obol group child-count sentinel should be created");
    struct ged_draw_scene_tree_summary group_only_tree;
    memset(&group_only_tree, 0, sizeof(group_only_tree));
    if (!ged_draw_group_ref_tree_summary(gedp, group_only_ref,
	    &group_only_tree) ||
	    !group_only_tree.valid ||
	    !group_only_tree.is_group ||
	    group_only_tree.is_shape ||
	    !group_only_tree.has_parent ||
	    !path_equal(group_only_tree.name, "group_only.s") ||
	    !group_only_tree.fullpath ||
	    !path_equal(DB_FULL_PATH_CUR_DIR(group_only_tree.fullpath)->d_namep,
		"group_only.s") ||
	    group_only_tree.draw_tree_depth != 1 ||
	    group_only_tree.child_count != 1)
	FAIL("GED group tree summaries should prefer owned Obol group tree metadata");
    ged_scene_node_ref group_only_ctx = ged_scene_group_node(gedp,
	group_only_ref);
    struct ged_draw_scene_tree_summary group_only_context_tree;
    memset(&group_only_context_tree, 0, sizeof(group_only_context_tree));
    if (ged_scene_node_ref_is_null(group_only_ctx) ||
	    !ged_scene_node_tree_summary(gedp, group_only_ctx,
		&group_only_context_tree) ||
	    group_only_context_tree.child_count != 1)
	FAIL("GED group-ref context should resolve to an owned Obol group context");
    ged_scene_node_ref group_only_parent_ctx =
	ged_scene_node_parent(gedp, group_only_ctx);
    struct ged_draw_scene_tree_summary group_only_parent_tree;
    memset(&group_only_parent_tree, 0, sizeof(group_only_parent_tree));
    if (ged_scene_node_ref_is_null(group_only_parent_ctx) ||
	    !ged_scene_node_tree_summary(gedp, group_only_parent_ctx,
		&group_only_parent_tree) ||
	    !group_only_parent_tree.valid ||
	    !group_only_parent_tree.is_group ||
	    group_only_parent_tree.has_parent ||
	    !path_equal(ged_scene_node_name(gedp, group_only_parent_ctx),
		"/"))
	FAIL("GED group-ref context should expose an owned Obol parent context");
    ged_scene_node_ref group_only_child_ctx =
	ged_scene_node_child_at(gedp, group_only_ctx, 0);
    struct ged_draw_scene_tree_summary group_only_child_tree;
    memset(&group_only_child_tree, 0, sizeof(group_only_child_tree));
    if (ged_scene_node_ref_is_null(group_only_child_ctx) ||
	    !ged_scene_node_tree_summary(gedp, group_only_child_ctx,
		&group_only_child_tree) ||
	    !group_only_child_tree.valid ||
	    !group_only_child_tree.is_group ||
	    !path_equal(ged_scene_node_name(gedp, group_only_child_ctx),
		"obol_child.s") ||
	    !ged_scene_node_ref_equal(
		ged_scene_node_parent(gedp, group_only_child_ctx),
		group_only_ctx))
	FAIL("GED scene-context child traversal should return owned Obol child contexts");
    if (owned_scene->setGroupDrawIntent("group_only.s",
	    "ged-draw-group:group_only.s", BOBOL_LOD_DRAW_WIRE,
	    BOBOL_LOD_DRAW_WIRE, TRUE, 0) < 0 ||
	    !ged_scene_group_is_overlay(gedp, group_only_ctx))
	FAIL("GED group contexts should read owned Obol overlay state");
    if (owned_scene->setGroupDrawIntent("group_only.s",
	    "ged-draw-group:group_only.s", BOBOL_LOD_DRAW_WIRE,
	    BOBOL_LOD_DRAW_WIRE, FALSE, 0) < 0 ||
	    ged_scene_group_is_overlay(gedp, group_only_ctx))
	FAIL("GED group contexts should clear owned Obol overlay state");
    if (owned_scene->removeGroup("group_only.s/obol_child.s") <= 0)
	FAIL("owned Obol group child-count sentinel should be removable");
    if (owned_scene->setGroupDrawIntent("group_only.s",
	    "ged-draw-group:group_only.s", BOBOL_LOD_DRAW_WIRE,
	    BOBOL_LOD_DRAW_WIRE, TRUE, 0) < 0)
	FAIL("owned Obol group overlay erase sentinel update should succeed");
    ged_draw_erase_name(gedp, "group_only.s");
    memset(&group_only_tree, 0, sizeof(group_only_tree));
    if (!ged_draw_group_ref_tree_summary(gedp, group_only_ref,
	    &group_only_tree) ||
	    !group_only_tree.valid ||
	    !group_only_tree.is_group)
	FAIL("GED public name erase should preserve overlay groups from owned Obol state");
    if (owned_scene->setGroupDrawIntent("group_only.s",
	    "ged-draw-group:group_only.s", BOBOL_LOD_DRAW_WIRE,
	    BOBOL_LOD_DRAW_WIRE, FALSE, 0) < 0)
	FAIL("owned Obol group overlay erase sentinel restore should succeed");
    struct ged_draw_scene_tree_summary original_root_count_tree;
    memset(&original_root_count_tree, 0, sizeof(original_root_count_tree));
    if (!ged_scene_node_tree_summary(gedp, group_only_parent_ctx,
	    &original_root_count_tree) ||
	    !original_root_count_tree.valid ||
	    !original_root_count_tree.is_group)
	FAIL("GED root scene context should summarize owned Obol root children");
    const size_t original_root_child_count =
	(size_t)original_root_count_tree.child_count;
    if (!ged_draw_obol_group_ensure_for_path(gedp,
	    "__obol_root_count_only.s", "__obol_root_count_only.s",
	    GED_DRAW_MODE_WIRE, 0))
	FAIL("GED Obol root child-count sentinel group should be ensured");
    struct ged_draw_scene_tree_summary updated_root_count_tree;
    memset(&updated_root_count_tree, 0, sizeof(updated_root_count_tree));
    if (owned_scene->getGroupChildCount("/") !=
	    (int)(original_root_child_count + 1) ||
	    !ged_scene_node_tree_summary(gedp, group_only_parent_ctx,
		&updated_root_count_tree) ||
	    updated_root_count_tree.child_count !=
	    original_root_child_count + 1)
	FAIL("GED root scene context child count should prefer owned Obol root children");
    if (owned_scene->removeGroup("__obol_root_count_only.s") <= 0)
	FAIL("GED Obol root child-count sentinel group should be removable");
    if (owned_scene->setDatabaseSourceState("box.s",
	    TRUE,
	    77,
	    88,
	    FALSE,
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
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "box.s", 0, 0, 0, 0, 0.0,
	0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &box_record);
    if (!box_record.found ||
	    box_record.sourceRevision != 77 ||
	    box_record.inputsRevision != 88)
	FAIL("GED shape records should read source summary state from the owned Obol controller");
    if (!path_equal(ged_draw_registry_shape_ref_semantic_path(gedp,
	    box_record.ref), "box.s"))
	FAIL("GED shape refs should retain a semantic registry path for Obol lookup");
    ged_scene_node_ref box_ctx = ged_scene_shape_node(gedp, box_record.ref);
    struct ged_draw_scene_tree_summary box_context_tree;
    memset(&box_context_tree, 0, sizeof(box_context_tree));
    if (ged_scene_node_ref_is_null(box_ctx) ||
	    !ged_scene_node_has_state(gedp, box_ctx) ||
	    !ged_scene_node_tree_summary(gedp, box_ctx,
		&box_context_tree) ||
	    !box_context_tree.valid ||
	    !box_context_tree.is_group ||
	    box_context_tree.child_count <= 0)
	FAIL("GED shape-ref context should resolve to an owned Obol database-source context");
    ged_scene_node_ref box_registry_ctx =
	ged_scene_shape_node(gedp, box_record.ref);
    struct ged_draw_scene_tree_summary box_registry_context_tree;
    memset(&box_registry_context_tree, 0, sizeof(box_registry_context_tree));
    if (ged_scene_node_ref_is_null(box_registry_ctx) ||
	    !ged_scene_node_tree_summary(gedp, box_registry_ctx,
		&box_registry_context_tree) ||
	    !box_registry_context_tree.valid ||
	    !box_registry_context_tree.is_group ||
	    box_registry_context_tree.is_shape ||
	    !box_registry_context_tree.has_parent ||
	    !path_equal(box_registry_context_tree.name, "box.s") ||
	    !box_registry_context_tree.fullpath ||
	    !path_equal(DB_FULL_PATH_CUR_DIR(
		    box_registry_context_tree.fullpath)->d_namep, "box.s") ||
	    box_registry_context_tree.draw_tree_depth !=
		box_context_tree.draw_tree_depth ||
	    box_registry_context_tree.child_count !=
		box_context_tree.child_count)
	FAIL("GED registry scene-context tree summaries should prefer owned Obol source metadata");
    ged_scene_node_ref box_registry_parent_ctx =
	ged_scene_node_parent(gedp, box_registry_ctx);
    struct ged_draw_scene_tree_summary box_registry_parent_tree;
    memset(&box_registry_parent_tree, 0,
	    sizeof(box_registry_parent_tree));
    if (ged_scene_node_ref_is_null(box_registry_parent_ctx) ||
	    !ged_scene_node_tree_summary(gedp, box_registry_parent_ctx,
		&box_registry_parent_tree) ||
	    !box_registry_parent_tree.valid ||
	    !box_registry_parent_tree.is_group ||
	    box_registry_parent_tree.is_shape ||
	    box_registry_parent_tree.has_parent ||
	    !path_equal(ged_scene_node_name(gedp,
		    box_registry_parent_ctx), "/"))
	FAIL("GED registry semantic source parents should resolve to owned Obol parent contexts");
    if (!ged_scene_node_ref_equal(ged_scene_node_source(gedp, box_ctx),
	    box_ctx))
	FAIL("GED Obol shape-ref context source should be the owned Obol database-source context");
    struct ged_draw_database_source_summary box_context_source;
    memset(&box_context_source, 0, sizeof(box_context_source));
    if (!ged_draw_source_node_summary(gedp, box_ctx,
		&box_context_source) ||
	    !box_context_source.valid ||
	    box_context_source.source_revision != 77 ||
	    box_context_source.inputs_revision != 88)
	FAIL("GED scene-context source summaries should read owned Obol source state");
    struct ged_draw_scene_display_summary box_context_display;
    memset(&box_context_display, 0, sizeof(box_context_display));
    if (!ged_scene_node_display_summary(gedp, box_ctx,
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
    SbVec3f sentinel_points[2] = {
	SbVec3f(11.0f, 0.0f, 0.0f),
	SbVec3f(12.0f, 0.0f, 0.0f)
    };
    int32_t sentinel_commands[2] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW
    };
    BObolExternalLineSet external_line;
    external_line.points = sentinel_points;
    external_line.commands = sentinel_commands;
    external_line.count = 2;
    if (box_source->publishExternalLineSet(external_line) <= 0)
	FAIL("shape-context test should publish explicit external line geometry");
    SoBRLVListShape *box_shape = box_source->getRealizedShape();
    if (!box_shape)
	FAIL("external line publication should expose realized VLIST geometry");
    ged_draw_index_stats_reset(gedp);
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
    struct ged_draw_view_line_summary box_line;
    memset(&box_line, 0, sizeof(box_line));
    if (!ged_draw_shape_ref_line_summary(gedp, box_record.ref, &box_line) ||
	    !box_line.valid ||
	    box_line.point_count != 2)
	FAIL("GED shape line summary should read realized VLIST state from the owned Obol controller");
    struct ged_draw_view_line_summary box_context_line;
    memset(&box_context_line, 0, sizeof(box_context_line));
    if (!ged_scene_node_line_summary(gedp, box_ctx, &box_context_line) ||
	    !box_context_line.valid ||
	    box_context_line.point_count != 2)
	FAIL("GED shape-context line summary should read realized VLIST state from the owned Obol controller");
    ged_scene_node_ref box_source_ctx = ged_scene_node_source(gedp, box_ctx);
    if (ged_scene_node_ref_is_null(box_source_ctx))
	FAIL("GED shape contexts should expose a source context for Obol traversal");
    ged_scene_node_ref box_child_ctx = ged_scene_node_child_at(gedp,
	box_source_ctx, 0);
    struct ged_draw_scene_tree_summary box_child_tree;
    memset(&box_child_tree, 0, sizeof(box_child_tree));
    if (ged_scene_node_ref_is_null(box_child_ctx))
	FAIL("GED source scene-context traversal should create owned Obol realized child contexts");
    if (!ged_scene_node_tree_summary(gedp, box_child_ctx,
	    &box_child_tree) ||
	    !box_child_tree.valid)
	FAIL("GED source scene-context traversal should summarize owned Obol realized children");
    if (!box_child_tree.is_shape)
	FAIL("GED source scene-context traversal should classify owned Obol realized children as shapes");
    if (!ged_scene_node_ref_equal(ged_scene_node_parent(gedp, box_child_ctx),
	    box_source_ctx))
	FAIL("GED source scene-context traversal should return owned Obol realized children");
    point_t box_line_point;
    if (!ged_draw_shape_ref_line_point_at(gedp, box_record.ref, 1,
	    box_line_point) ||
	    fabs(box_line_point[0] - 12.0) > 0.001 ||
	    fabs(box_line_point[1]) > 0.001 ||
	    fabs(box_line_point[2]) > 0.001)
	FAIL("GED shape line point readback should read realized VLIST points from the owned Obol controller");
    point_t box_context_line_point;
    if (!ged_scene_node_line_point_at(gedp, box_ctx, 1,
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
    if (!ged_scene_node_line_command_at(gedp, box_ctx, 1,
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
    if (!ged_scene_node_geometry_summary(gedp, box_ctx,
	    &box_context_geometry) ||
	    !box_context_geometry.valid ||
	    !box_context_geometry.geometry_name ||
	    !BU_STR_EQUAL(box_context_geometry.geometry_name, "line-set") ||
	    box_context_geometry.point_count != 2 ||
	    box_context_geometry.index_count != 0)
	FAIL("GED shape-context geometry summary should read realized geometry from the owned Obol controller");
	    vect_t draw_bounds_min;
	    vect_t draw_bounds_max;
	    if (ged_scene_bounds(gedp, &draw_bounds_min, &draw_bounds_max,
		    GED_SCENE_BOUNDS_DATABASE) ||
		    draw_bounds_max[0] < 11.9)
		FAIL("GED draw bounds should read database-source bounds from the owned Obol controller");
	    vect_t context_bounds_min;
	    vect_t context_bounds_max;
	    if (ged_scene_node_subtree_bounds(gedp, box_ctx,
		    &context_bounds_min, &context_bounds_max, 0) ||
		    context_bounds_max[0] < 11.9)
		FAIL("GED source contexts should read owned Obol subtree bounds");
	    if (owned_scene->moveDatabaseSourceToGroup("box.s",
		    "group_only.s") != 1)
		FAIL("owned Obol source should move to the sentinel group for context bounds");
	    if (owned_scene->setGroupDrawIntent("group_only.s",
		    "ged-draw-group:group_only.s", BOBOL_LOD_DRAW_WIRE,
		    BOBOL_LOD_DRAW_WIRE, TRUE, 0) < 0 ||
		    !ged_scene_group_is_overlay(gedp, group_only_ctx))
		FAIL("owned Obol group overlay state should remain authoritative before bounds");
	    if (!ged_scene_node_subtree_bounds(gedp, group_only_ctx,
		    &context_bounds_min, &context_bounds_max, 0))
		FAIL("GED group context bounds should skip owned Obol overlay groups");
	    if (owned_scene->setGroupDrawIntent("group_only.s",
		    "ged-draw-group:group_only.s", BOBOL_LOD_DRAW_WIRE,
		    BOBOL_LOD_DRAW_WIRE, FALSE, 0) < 0 ||
		    ged_scene_group_is_overlay(gedp, group_only_ctx))
		FAIL("owned Obol group overlay state should clear before bounds");
	    if (ged_scene_node_subtree_bounds(gedp, group_only_ctx,
		    &context_bounds_min, &context_bounds_max, 0) ||
		    context_bounds_max[0] < 11.9)
		FAIL("GED group contexts should read owned Obol subtree bounds");
	    if (owned_scene->moveDatabaseSourceToGroup("box.s", "/") != 1)
		FAIL("owned Obol source should move back to the root group after context bounds");
	    vect_t box_xlate;
	    VSET(box_xlate, 5.0, 0.0, 0.0);
	    if (!ged_scene_internal_shape_translate_geometry(gedp, box_record.ref,
		    box_xlate))
		FAIL("GED shape geometry translation should succeed");
	    if (!ged_draw_shape_ref_line_point_at(gedp, box_record.ref, 1,
		    box_line_point) ||
		    fabs(box_line_point[0] - 17.0) > 0.001 ||
		    fabs(box_line_point[1]) > 0.001 ||
		    fabs(box_line_point[2]) > 0.001)
		FAIL("GED shape translation should mutate owned Obol VLIST points");
	    if (ged_scene_bounds(gedp, &draw_bounds_min, &draw_bounds_max,
		    GED_SCENE_BOUNDS_DATABASE) ||
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
	    if (!ged_draw_obol_database_source_publish_line_set_for_path(
		    gedp, "box.s", (const point_t *)published_points,
		    published_commands, 3))
		FAIL("GED Obol source line-set publish should succeed");
	    if (!ged_draw_shape_ref_line_summary(gedp, box_record.ref,
		    &box_line) ||
		    !box_line.valid ||
		    box_line.point_count != 3)
		FAIL("GED Obol source line-set publish should update owned Obol VLIST count");
	    if (!ged_draw_shape_ref_line_point_at(gedp, box_record.ref, 2,
		    box_line_point) ||
		    fabs(box_line_point[0] - 23.0) > 0.001 ||
		    fabs(box_line_point[1]) > 0.001 ||
		    fabs(box_line_point[2]) > 0.001)
		FAIL("GED Obol source line-set publish should update owned Obol VLIST points");
	    if (ged_scene_bounds(gedp, &draw_bounds_min, &draw_bounds_max,
		    GED_SCENE_BOUNDS_DATABASE) ||
		    draw_bounds_max[0] < 22.9)
		FAIL("GED draw bounds should reflect published owned Obol VLIST points");
	    box_shape = box_source->getRealizedShape();
	    if (!box_shape)
		FAIL("GED source line-set publish should retain realized VLIST geometry");
	    point_t explicit_center;
	    VSET(explicit_center, 30.0, 31.0, 32.0);
	    ged_draw_index_stats_reset(gedp);
	    if (!ged_scene_internal_shape_set_center(gedp, box_record.ref,
		    explicit_center))
		FAIL("GED shape center setter should succeed");
	    SbVec3f obol_center = box_shape->drawCenter.getValue();
	    if (!box_shape->drawCenterValid.getValue() ||
		    fabs(obol_center[0] - 30.0f) > 0.001f ||
		    fabs(obol_center[1] - 31.0f) > 0.001f ||
		    fabs(obol_center[2] - 32.0f) > 0.001f)
		FAIL("GED shape center setter should mutate owned Obol VLIST center");
	    BObolDatabaseSourceSummary box_placement_summary;
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
	    if (!ged_scene_internal_shape_bounds_update(gedp,
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
	    if (!ged_scene_internal_shape_geometry_clear(gedp, box_record.ref))
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
	    if (ged_scene_bounds(gedp, &draw_bounds_min, &draw_bounds_max,
		    GED_SCENE_BOUNDS_DATABASE) ||
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
	    struct ged_draw_scene_display_summary aux_display_state;
	    memset(&aux_display_state, 0, sizeof(aux_display_state));
	    aux_display_state.valid = 1;
	    aux_display_state.draw_mode = GED_DRAW_MODE_SHADED;
	    aux_display_state.visible = 0;
	    aux_display_state.highlighted = 1;
	    aux_display_state.line_style = 9;
	    aux_display_state.line_width = 11;
	    aux_display_state.transparency = 0.25;
	    aux_display_state.material_valid = 1;
	    aux_display_state.material_color[0] = 80;
	    aux_display_state.material_color[1] = 120;
	    aux_display_state.material_color[2] = 160;
	    if (!ged_draw_obol_database_source_publish_auxiliary_line_set_for_path(
		    gedp, "box.s", "obol_aux_clear_sentinel",
		    (const point_t *)aux_points, aux_commands, 2,
		    &aux_display_state))
		FAIL("GED Obol auxiliary line-set bridge should publish to owned source");
	    SoBRLVListShape *aux_shape =
		box_source->findAuxiliaryVListShape(
			"obol_aux_clear_sentinel");
	    SbColor aux_material_color;
	    if (aux_shape)
		aux_material_color = aux_shape->materialColor.getValue();
	    if (!aux_shape ||
		    aux_shape->point.getNum() != 2 ||
		    aux_shape->command.getNum() != 2 ||
		    bu_strcmp(aux_shape->recordRole.getValue().getString(),
			"auxiliary") != 0 ||
		    aux_shape->drawMode.getValue() !=
			BOBOL_LOD_DRAW_SHADED ||
		    aux_shape->visible.getValue() ||
		    !aux_shape->highlighted.getValue() ||
		    aux_shape->lineStyle.getValue() != 9 ||
		    aux_shape->lineWidth.getValue() != 11 ||
		    fabs(aux_shape->transparency.getValue() - 0.25) >
			0.001 ||
		    !aux_shape->materialColorValid.getValue() ||
		    fabs(aux_material_color[0] - (80.0f / 255.0f)) >
			0.001f ||
		    fabs(aux_material_color[1] - (120.0f / 255.0f)) >
			0.001f ||
		    fabs(aux_material_color[2] - (160.0f / 255.0f)) >
			0.001f)
		FAIL("GED Obol auxiliary line-set bridge should create an auxiliary VLIST with display state");
	    if (!aux_shape->drawMatrixValid.getValue() ||
		    !aux_shape->drawCenterValid.getValue() ||
		    !aux_shape->drawSizeValid.getValue())
		FAIL("GED Obol auxiliary line-set bridge should inherit owned source placement");
	    if (!ged_scene_internal_shape_geometry_clear(gedp, box_record.ref) ||
		    box_source->findAuxiliaryVListShape(
			"obol_aux_clear_sentinel"))
		FAIL("GED shape geometry clear should clear owned Obol auxiliary VLISTs");
	    point_t point_set_points[2] = {
		{24.0, 0.0, 0.0},
		{25.0, 1.0, 0.0}
	    };
	    if (!ged_draw_obol_database_source_publish_point_set_for_path(
		    gedp, "box.s", (const point_t *)point_set_points, 2))
		FAIL("GED Obol source point-set publish should succeed");
	    box_shape = box_source->getRealizedShape();
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
	    if (!ged_scene_internal_shape_geometry_clear(gedp, box_record.ref))
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
	    if (!ged_draw_obol_database_source_publish_indexed_face_set_for_path(
		    gedp, "box.s", (const point_t *)mesh_points, 4,
		    (const vect_t *)mesh_normals, 4, mesh_indices, 5))
		FAIL("GED Obol source indexed-face publish should succeed");
	    SoBRLMeshShape *box_mesh = box_source->getRealizedMesh();
	    if (!box_mesh ||
		    box_mesh->point.getNum() != 4 ||
		    box_mesh->coordIndex.getNum() != 6 ||
		    box_mesh->normal.getNum() != 6 ||
		    fabs(box_mesh->normal[0][2] - 1.0f) > 0.001f ||
		    fabs(box_mesh->normal[5][2] - 1.0f) > 0.001f ||
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
	    if (ged_scene_bounds(gedp, &draw_bounds_min, &draw_bounds_max,
		    GED_SCENE_BOUNDS_DATABASE) ||
		    draw_bounds_max[0] < 31.9)
		FAIL("GED draw bounds should reflect published owned Obol mesh points");
	    if (!ged_scene_internal_shape_geometry_clear(gedp, box_record.ref))
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
	    ged_draw_index_stats_reset(gedp);
	    if (ged_draw_obol_database_source_publish_primitive_wireframe_for_path(
		    gedp, "box.s", &box_intern, NULL, NULL) <= 0)
		FAIL("GED primitive wireframe publication should succeed");
	    rt_db_free_internal(&box_intern);
	    struct ged_draw_index_stats primitive_wire_stats;
	    memset(&primitive_wire_stats, 0, sizeof(primitive_wire_stats));
	    ged_draw_index_stats_get(gedp, &primitive_wire_stats);
	    if (primitive_wire_stats.path_queries ||
		    primitive_wire_stats.path_candidates)
		FAIL("GED Obol primitive wireframe publication should avoid registry path-index queries");
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
	    if (!ged_scene_internal_shape_geometry_clear(gedp, box_record.ref))
		FAIL("GED shape geometry clear should clear primitive wireframe publication");
	    const char *draw_box_shaded[3] = {"draw", "-m2", "box.s"};
	    if (ged_exec_draw(gedp, 3, draw_box_shaded) != BRLCAD_OK)
		FAIL("GED shaded draw should succeed for source realization test");
	    box_record.found = 0;
	    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
		    &box_record);
	    if (!box_record.found)
		FAIL("GED shape record should refresh after shaded draw source realization");
	    box_source = source_for_representation(owned_scene, "box.s",
		    SoBRLDatabaseSource::REPRESENTATION_SHADED);
	    if (!box_source) {
		box_source = owned_scene->findDatabaseSource("/box.s");
		if (!box_source)
		    box_source = source_for_path(owned_scene, "box.s");
	    }
	    if (!box_source)
		FAIL("owned Obol source should remain available after shaded draw source realization");
	    BObolRealizedShapeSummary box_mesh_summary;
	    /*
	     * With an attached progressive endpoint, draw success publishes the
	     * coarse source contract; detached mesh realization is intentionally
	     * asynchronous.  Wait on the actual carrier-free mesh predicate and
	     * service requested frames instead of assuming a particular worker
	     * scheduling speed (notably under TSan).
	     */
	    BObolViewController *box_view_controller = owned_controller;
	    for (int attempt = 0; attempt < 2000; attempt++) {
		if (box_source->getRealizedMeshCount() == 0 &&
		    box_source->getRealizedShapeSummary(0,
			box_mesh_summary) &&
		    box_mesh_summary.shapeKind ==
			BObolRealizedShapeSummary::SHAPE_MESH &&
		    box_mesh_summary.pointCount > 0 &&
		    box_mesh_summary.indexCount > 0)
		    break;
		if (box_view_controller) {
		    (void)box_view_controller->realizePending();
		    BObolProgressiveStatus shaded_status;
		    (void)box_view_controller->advanceProgressiveWork(
			NULL, &shaded_status);
		    SbString render_reason;
		    if (box_view_controller->consumeRenderRequest(
			    &render_reason)) {
			const uint64_t frame_started =
			    box_view_controller->beginRenderTiming();
			box_view_controller->completeRenderTiming(
			    frame_started, test_capacity_cad_timing());
		    }
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
		SoBRLDatabaseSource *updated_source =
		    source_for_representation(owned_scene, "box.s",
		    SoBRLDatabaseSource::REPRESENTATION_SHADED);
		if (!updated_source)
		    updated_source =
			owned_scene->findDatabaseSource("/box.s");
		if (!updated_source)
		    updated_source = source_for_path(owned_scene, "box.s");
		if (updated_source)
		    box_source = updated_source;
	    }
	    if (!box_source ||
		box_source->getRealizedMeshCount() != 0 ||
		!box_source->getRealizedShapeSummary(0, box_mesh_summary) ||
		box_mesh_summary.shapeKind !=
		    BObolRealizedShapeSummary::SHAPE_MESH ||
		box_mesh_summary.pointCount == 0 ||
		box_mesh_summary.indexCount == 0) {
		BObolDatabaseSourceSummary failed_summary;
		const int have_failed_summary =
		    box_source && box_source->getSummary(failed_summary);
		fprintf(stderr,
		    "shaded realization source=%p summary=%d status=%d "
		    "rep=%d compact=%d instances=%d shapes=%d meshes=%d "
		    "kind=%d points=%d indices=%d scene_sources=%d "
		    "diagnostic=%s\n",
		    (void *)box_source, have_failed_summary,
		    have_failed_summary ?
			failed_summary.realizationStatus : -1,
		    have_failed_summary ?
			failed_summary.representationMode : -1,
		    box_source ?
			box_source->isCompactOccurrenceRegistry() : -1,
		    box_source ? box_source->getCompactInstanceCount() : -1,
		    box_source ? box_source->getRealizedShapeCount() : -1,
		    box_source ? box_source->getRealizedMeshCount() : -1,
		    box_mesh_summary.shapeKind, box_mesh_summary.pointCount,
		    box_mesh_summary.indexCount,
		    owned_scene->getDatabaseSourceCount(),
		    have_failed_summary ?
			failed_summary.realizationDiagnostic.getString() : "");
		FAIL("GED shaded source realization should publish carrier-free Obol mesh geometry");
	    }
	    BObolDatabaseSourceSummary box_realized_summary;
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
	    if (box_source->prepareCompiledAssembly() != 1)
		FAIL("GED Obol source placement bridge should retain a compiled assembly");
	    SoGetBoundingBoxAction placement_bounds(SbViewportRegion(200, 200));
	    placement_bounds.apply(box_source);
	    const SbBox3f placed_box = placement_bounds.getBoundingBox();
	    if (placed_box.isEmpty() || placed_box.getMin()[0] > 23.1f ||
		placed_box.getMax()[0] < 24.9f)
		FAIL("GED Obol source placement bridge should transform compact geometry bounds");
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
	    ged_draw_index_stats_reset(gedp);
	    if (ged_draw_shape_ref_redraw_wireframe(gedp, box_record.ref,
		    gedp->dbip, NULL, NULL, NULL, 0) < 0)
		FAIL("GED wire redraw should succeed with owned Obol draw matrix");
	    struct ged_draw_index_stats obol_wire_redraw_stats;
	    memset(&obol_wire_redraw_stats, 0,
		    sizeof(obol_wire_redraw_stats));
	    ged_draw_index_stats_get(gedp, &obol_wire_redraw_stats);
	    if (obol_wire_redraw_stats.path_queries ||
		    obol_wire_redraw_stats.path_candidates)
		FAIL("GED Obol wire redraw should avoid registry path-index queries");
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
	    if (ged_scene_bounds(gedp, &draw_bounds_min, &draw_bounds_max,
		    GED_SCENE_BOUNDS_DATABASE) ||
		    draw_bounds_max[0] < 40.9)
		FAIL("GED redraw should use owned Obol draw matrix");
	    if (owned_scene->setDatabaseSourcePlacementState(
		    box_source_path.c_str(), FALSE, SbMatrix::identity(),
		    FALSE, SbVec3f(0.0f, 0.0f, 0.0f), FALSE, 0.0f) < 0)
		FAIL("owned Obol source draw matrix reset should succeed");
	    ged_draw_index_stats_reset(gedp);
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
	    struct ged_view_context *lod_view_ctx = ged_draw_shape_ref_view_context(gedp,
		    box_record.ref);
	    if (!lod_view_ctx)
		FAIL("GED LoD view context should be available");
	    if (!bv_scale_set(DRAW_TEST_BV(lod_view_ctx), 7.0))
		FAIL("GED LoD view context scale should be settable");
	    ged_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
	    lod_policy.csg_enabled = 1;
	    lod_policy.mesh_enabled = 0;
	    lod_policy.scale = 1.75;
	    lod_policy.bot_threshold = 77;
	    lod_policy.curve_scale = 2.25;
	    lod_policy.point_scale = 3.25;
	    if (!ged_view_lod_policy_apply(lod_view_ctx, &lod_policy))
		FAIL("GED LoD view policy should be settable");
	    struct ged_view_context *lod_view_ctxs[1] = {lod_view_ctx};
	    ged_draw_index_stats_reset(gedp);
	    if (!ged_draw_shape_ref_lod_ensure(gedp, box_record.ref,
		    lod_view_ctx, lod_view_ctxs, 1))
		FAIL("GED LoD ensure should succeed for Obol source policy test");
	    if (!box_source->getSummary(box_realized_summary) ||
		    !(box_realized_summary.realizationRoleFlags &
			SoBRLDatabaseSource::REALIZATION_ROLE_CSG) ||
		    !box_realized_summary.realizationViewDependent ||
		    fabs(box_realized_summary.realizationViewScale - 7.0f) >
			0.001 ||
		    fabs(box_realized_summary.realizationLodScale - 1.75f) >
			0.001 ||
		    box_realized_summary.realizationBotThreshold != 77 ||
		    fabs(box_realized_summary.realizationCurveScale - 2.25f) >
			0.001 ||
		    fabs(box_realized_summary.realizationPointScale - 3.25f) >
			0.001) {
		fprintf(stderr,
			"LoD Obol summary role=%d view=%d view_scale=%g lod_scale=%g bot=%u curve=%g point=%g\n",
			box_realized_summary.realizationRoleFlags,
			box_realized_summary.realizationViewDependent ? 1 : 0,
			(double)box_realized_summary.realizationViewScale,
			(double)box_realized_summary.realizationLodScale,
			(unsigned)box_realized_summary.realizationBotThreshold,
			(double)box_realized_summary.realizationCurveScale,
			(double)box_realized_summary.realizationPointScale);
		FAIL("GED LoD realization policy should update owned Obol source state");
	    }
	    struct ged_draw_shape_geometry_summary box_csg_geometry;
	    memset(&box_csg_geometry, 0, sizeof(box_csg_geometry));
	    if (!ged_draw_shape_ref_geometry_summary(gedp, box_record.ref,
		    &box_csg_geometry) ||
		    !box_csg_geometry.valid ||
		    !box_csg_geometry.geometry_name ||
		    !BU_STR_EQUAL(box_csg_geometry.geometry_name, "line-set") ||
		    box_csg_geometry.point_count == 0)
		FAIL("GED Obol adaptive CSG LoD draw should publish owned line geometry");
	    box_record.found = 0;
	    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
		    &box_record);
	    if (!box_record.found)
		FAIL("GED shape record should refresh after LoD policy test");
	    if (!ged_scene_internal_shape_geometry_clear(gedp, box_record.ref))
		FAIL("GED shape geometry clear should clear shaded source realization mesh");
	    ged_draw_index_stats_reset(gedp);
	    if (!ged_draw_shape_ref_set_visible(gedp, box_record.ref, 1))
		FAIL("GED visible setter should succeed");
	    struct ged_scene_path_request highlight_request;
	    ged_scene_path_request_init(&highlight_request);
	    highlight_request.path = "box.s";
	    ged_scene_occurrence_ref box_occurrence =
		ged_scene_occurrence_resolve(gedp, &highlight_request);
	    if (ged_scene_occurrence_ref_is_null(box_occurrence) ||
		ged_scene_occurrence_highlight_set(gedp, box_occurrence, 0,
		    NULL) != GED_SCENE_OK)
	FAIL("GED highlighted setter should succeed");
    struct ged_scene_occurrence_info cleared_occurrence;
    if (!ged_scene_occurrence_get(gedp, box_occurrence, &cleared_occurrence) ||
	!path_equal(cleared_occurrence.path, "box.s"))
	FAIL("cleared noncompact source should remain addressable after unrelated records");
    const unsigned char override_color[3] = {10, 20, 30};
    if (!ged_draw_shape_ref_set_color(gedp, box_record.ref, override_color))
	FAIL("GED color setter should succeed");
    box_source = source_for_path(owned_scene, "box.s");
    BObolDatabaseSourceSummary box_source_summary;
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
    if (!ged_draw_source_snapshot(gedp, box_record.ref,
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
    struct ged_view_context *stale_lod_view_ctx = ged_draw_shape_ref_view_context(gedp,
	    stale_box_ref);
    if (!stale_lod_view_ctx)
	FAIL("GED stale shape-ref view context should recover cached source state");
    if (!bv_scale_set(DRAW_TEST_BV(stale_lod_view_ctx), 9.0))
	FAIL("GED stale shape-ref LoD context scale should be settable");
    ged_view_lod_policy stale_lod_policy = BV_LOD_POLICY_INIT;
    stale_lod_policy.csg_enabled = 1;
    stale_lod_policy.mesh_enabled = 0;
    stale_lod_policy.scale = 2.75;
    stale_lod_policy.bot_threshold = 91;
    stale_lod_policy.curve_scale = 6.25;
    stale_lod_policy.point_scale = 7.25;
    if (!ged_view_lod_policy_apply(stale_lod_view_ctx,
	    &stale_lod_policy))
	FAIL("GED stale shape-ref LoD policy should be settable");
    struct ged_view_context *stale_lod_view_ctxs[1] = {stale_lod_view_ctx};
    if (!ged_draw_shape_ref_lod_ensure(gedp, stale_box_ref,
	    stale_lod_view_ctx, stale_lod_view_ctxs, 1))
	FAIL("GED stale shape-ref LoD ensure should recover cached Obol source runtime");
    if (!box_source->getSummary(box_realized_summary) ||
	    !(box_realized_summary.realizationRoleFlags &
		SoBRLDatabaseSource::REALIZATION_ROLE_CSG) ||
	    !box_realized_summary.realizationViewDependent ||
	    fabs(box_realized_summary.realizationViewScale - 9.0f) >
		0.001 ||
	    fabs(box_realized_summary.realizationLodScale - 2.75f) >
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
    ged_draw_index_stats_reset(gedp);
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
		BOBOL_LOD_DRAW_SHADED)
	FAIL("GED group mode setter should mutate the owned Obol group draw intent");
    if (!ged_draw_group_ref_set_mode(gedp, group_state.ref,
	    GED_DRAW_MODE_WIRE))
	FAIL("GED group mode restore should succeed");
    std::string obol_group_intent_path =
	std::string("ged-draw-group:") + group_path + "_intent";
    if (owned_scene->setGroupDrawIntent(group_path.c_str(),
	    obol_group_intent_path.c_str(), BOBOL_LOD_DRAW_SHADED,
	    BOBOL_LOD_DRAW_WIRE, TRUE, 501) <= 0)
	FAIL("owned Obol group draw-intent sentinel update should succeed");
    struct ged_draw_group_record intent_record;
    memset(&intent_record, 0, sizeof(intent_record));
    std::string obol_group_record_path = group_path + "_intent";
    if (!ged_draw_group_record_get(gedp, group_state.ref, &intent_record) ||
	    !intent_record.path ||
	    !BU_STR_EQUAL(intent_record.path,
		obol_group_record_path.c_str()) ||
	    intent_record.draw_mode != GED_DRAW_MODE_SHADED ||
	    !intent_record.is_overlay)
	FAIL("GED group records should read owned Obol draw-intent state");
    std::string original_obol_group_intent_path =
	std::string("ged-draw-group:") + group_path;
    if (owned_scene->setGroupDrawIntent(group_path.c_str(),
	    original_obol_group_intent_path.c_str(), BOBOL_LOD_DRAW_WIRE,
	    BOBOL_LOD_DRAW_WIRE, FALSE, 0) <= 0)
	FAIL("owned Obol group draw-intent sentinel restore should succeed");
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
	    fabs(scene_group->transparency.getValue() - 0.55) > 0.001 ||
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
	fabs(group_appearance_readback.transparency - 0.77) > 0.001 ||
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
	    0.55f,
	    TRUE,
	    SbColor(90.0f / 255.0f, 100.0f / 255.0f, 110.0f / 255.0f),
	    scene_group->materialColorValid.getValue(),
	    scene_group->materialColor.getValue(),
	    scene_group->materialRevision.getValue()) <= 0)
	FAIL("owned Obol group appearance sentinel restore should succeed");
    struct ged_draw_group_record group_record;
    memset(&group_record, 0, sizeof(group_record));
    if (!ged_draw_group_record_get(gedp, group_state.ref, &group_record) ||
	fabs(group_record.transparency - 0.55) > 0.001 ||
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
    /* Group mode/appearance mutations invalidate the retained source just as
     * they do in a GUI frame.  This headless test owns the frame pump, so
     * realize that pending state before asserting the current-record view. */
    (void)owned_controller->realizePending();
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
    source_record.realization_csg_lod_enabled = 0;
    source_record.realization_mesh_lod_enabled = 1;
    source_record.realization_view_scale = 11.0;
    source_record.realization_lod_scale = 3.5;
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
	    box_source_summary.realizationCsgLodEnabled ||
	    !box_source_summary.realizationMeshLodEnabled ||
	    fabs(box_source_summary.realizationViewScale - 11.0f) > 0.001f ||
	    fabs(box_source_summary.realizationLodScale - 3.5f) > 0.001f ||
	    box_source_summary.realizationBotThreshold != 88 ||
	    fabs(box_source_summary.realizationCurveScale - 4.5f) > 0.001f ||
	    fabs(box_source_summary.realizationPointScale - 5.5f) > 0.001f)
	FAIL("GED Obol source-record bridge should mutate owned source metadata");
    if (!ged_draw_obol_database_source_realize_for_path(gedp, "box.s"))
	FAIL("GED Obol source realization should realize the owned source");
    if (!box_source->getSummary(box_source_summary) ||
	    box_source_summary.stale)
	FAIL("GED shape realize-context should make the owned Obol source current before stale mutation check");
    if (ged_draw_source_mark_changed(gedp, "box.s",
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
    if (!source_for_path(owned_scene, "box.s") ||
	    source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror erase transactions");

    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller redraw command should succeed");
    if (!source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw transactions");
    if (owned_scene->replaceDatabaseSource("draft_move.s", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 6060) <= 0 ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("owned Obol transaction canary source should be created");
    ged_draw_index_stats_reset(gedp);
    struct ged_scene_reducer_request obol_index_txn =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_SOURCE_UPDATED,
		"draft_move.s");
    obol_index_txn.redraw = 0;
    struct ged_scene_reducer_result obol_index_result;
    ged_scene_reducer_result_init(&obol_index_result);
    if (ged_scene_reduce(gedp, &obol_index_txn,
	    &obol_index_result) <= 0)
	FAIL("GED Obol component-index stale transaction should succeed");
    ged_scene_reducer_result_free(&obol_index_result);
    struct ged_draw_index_stats obol_index_stats;
    memset(&obol_index_stats, 0, sizeof(obol_index_stats));
    ged_draw_index_stats_get(gedp, &obol_index_stats);
    if (obol_index_stats.slow_path_shape_scans ||
	    obol_index_stats.slow_path_group_scans)
	FAIL("GED Obol component indexes should avoid registry/index slow-path scans");
    SoBRLDatabaseSource *draft_move_source =
	source_for_path(owned_scene, "draft_move.s");
    BObolDatabaseSourceSummary draft_move_summary;
    if (!draft_move_source ||
	    !draft_move_source->getSummary(draft_move_summary) ||
	    !draft_move_summary.stale)
	FAIL("GED Obol component-index transaction should mark the owned source stale");
    struct ged_scene_reducer_request display_txn =
	ged_scene_reducer_request_make_value(GED_SCENE_REDUCER_VISIBILITY,
		"box.s", 0.0);
    struct ged_scene_reducer_result display_result;
    ged_scene_reducer_result_init(&display_result);
    if (ged_scene_reduce(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED visibility transaction should succeed");
    ged_scene_reducer_result_free(&display_result);
    box_source = source_for_path(owned_scene, "box.s");
    if (!box_source ||
	    !box_source->getSummary(box_source_summary) ||
	    box_source_summary.visible ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED visibility transaction should update owned Obol state without full-scene sync");
    display_txn = ged_scene_reducer_request_make_value(GED_SCENE_REDUCER_VISIBILITY,
	    "box.s", 1.0);
    ged_scene_reducer_result_init(&display_result);
    if (ged_scene_reduce(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED visibility restore transaction should succeed");
    ged_scene_reducer_result_free(&display_result);
    display_txn = ged_scene_reducer_request_make_value(GED_SCENE_REDUCER_TRANSPARENCY,
	    "box.s", 0.125);
    ged_scene_reducer_result_init(&display_result);
    if (ged_scene_reduce(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED transparency transaction should succeed");
    ged_scene_reducer_result_free(&display_result);
    if (!box_source->getSummary(box_source_summary) ||
	    fabs(box_source_summary.transparency - 0.125f) > 0.001f ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED transparency transaction should preserve Obol-only sources");
    display_txn = ged_scene_reducer_request_make_value(GED_SCENE_REDUCER_HIGHLIGHT,
	    "box.s", 1.0);
    ged_scene_reducer_result_init(&display_result);
    if (ged_scene_reduce(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED highlight transaction should succeed");
    ged_scene_reducer_result_free(&display_result);
    if (!box_source->getSummary(box_source_summary) ||
	    !box_source_summary.highlighted ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED highlight transaction should preserve Obol-only sources");
    display_txn = ged_scene_reducer_request_make_value(GED_SCENE_REDUCER_HIGHLIGHT,
	    "box.s", 0.0);
    ged_scene_reducer_result_init(&display_result);
    if (ged_scene_reduce(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED highlight restore transaction should succeed");
    ged_scene_reducer_result_free(&display_result);
    display_txn = ged_scene_reducer_request_make(GED_SCENE_REDUCER_STALE_SOURCE,
	    "box.s");
    display_txn.stale_reason = GED_DRAW_STALE_SETTINGS_CHANGED;
    ged_scene_reducer_result_init(&display_result);
    if (ged_scene_reduce(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED stale-source transaction should succeed");
    ged_scene_reducer_result_free(&display_result);
    if (!box_source->getSummary(box_source_summary) ||
	    !box_source_summary.stale ||
	    !(box_source_summary.staleReason &
		SoBRLDatabaseSource::STALE_DRAW) ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED stale-source transaction should target owned Obol state without full-scene sync");
    display_txn = ged_scene_reducer_request_make(GED_SCENE_REDUCER_MATERIAL_CHANGED,
	    NULL);
    ged_scene_reducer_result_init(&display_result);
    if (ged_scene_reduce(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED material-changed transaction should succeed");
    ged_scene_reducer_result_free(&display_result);
    if (!source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED material-changed transaction should preserve Obol-only sources");
    display_txn = ged_scene_reducer_request_make(GED_SCENE_REDUCER_REDRAW, NULL);
    ged_scene_reducer_result_init(&display_result);
    if (ged_scene_reduce(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED redraw transaction should succeed");
    ged_scene_reducer_result_free(&display_result);
    if (!source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s") ||
	    !source_for_path(owned_scene, "draft_move.s"))
	FAIL("GED redraw transaction should refresh draw sources without clearing Obol-only sources");
    display_txn = ged_scene_reducer_request_make(GED_SCENE_REDUCER_ERASE_PREFIX,
	    "box.s");
    ged_scene_reducer_result_init(&display_result);
    if (ged_scene_reduce(gedp, &display_txn,
	    &display_result) <= 0)
	FAIL("GED erase-prefix transaction should succeed");
    ged_scene_reducer_result_free(&display_result);
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
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "ball.s", 0, 0, 0, 0, 0.0,
	0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &ball_record);
    if (!ball_record.found)
	FAIL("GED ball shape record should be available before public erase");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, "ball.s", NULL, -1,
	    "public erase"))
	return 1;
    if (owned_scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(owned_scene, "box.s") ||
	    source_for_path(owned_scene, "ball.s"))
	FAIL("GED public erase should remove the owned Obol source");
    memset(&ball_record, 0, sizeof(ball_record));
    ball_record.matchPath = "ball.s";
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &ball_record);
    if (ball_record.found)
	FAIL("GED public erase should not expose stale shape draw records");
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller public-erase redraw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after public erase");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, "ball.s",
	    ged_draw_active_view_ctx(gedp), -1, "public active-scope erase"))
	return 1;
    if (owned_scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(owned_scene, "box.s") ||
	    source_for_path(owned_scene, "ball.s"))
	FAIL("GED public active-scope erase should remove the owned Obol group subtree");
    memset(&ball_record, 0, sizeof(ball_record));
    ball_record.matchPath = "ball.s";
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &ball_record);
    if (ball_record.found)
	FAIL("GED active-scope erase should not expose stale shape draw records");
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller active-scope erase redraw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
		FAIL("owned Obol scene controller should mirror redraw after public active-scope erase");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, "ball.s", NULL, -1,
	    "public root erase"))
	return 1;
    if (owned_scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(owned_scene, "box.s") ||
	    source_for_path(owned_scene, "ball.s"))
	FAIL("GED public root erase should remove the owned Obol source");
    memset(&ball_record, 0, sizeof(ball_record));
    ball_record.matchPath = "ball.s";
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &ball_record);
    if (ball_record.found)
	FAIL("GED public root erase should not expose stale shape draw records");
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller public-root-erase redraw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after public root erase");
    ged_draw_erase_name(gedp, "ball.s");
    if (owned_scene->getDatabaseSourceCount() != 1 ||
	    !source_for_path(owned_scene, "box.s") ||
	    source_for_path(owned_scene, "ball.s"))
	FAIL("GED public name erase should remove the owned Obol group");
    memset(&ball_record, 0, sizeof(ball_record));
    ball_record.matchPath = "ball.s";
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &ball_record);
    if (ball_record.found)
	FAIL("GED public name erase should not expose stale shape draw records");
    if (ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller public-name-erase redraw command should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after public name erase");
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
    (void)owned_controller->realizePending();
    SoBRLDatabaseSource *nested_leaf_source =
	source_for_path(owned_scene, nested_leaf_source_path);
    SoBRLDatabaseSource *nested_sibling_source =
	source_for_path(owned_scene, nested_sibling_source_path);
    if (!nested_leaf_source || !nested_sibling_source ||
	    owned_scene->getDatabaseSourceCount() != 4)
	FAIL("GED nested path draws should create owned Obol child sources");
    BObolDatabaseSourceSummary nested_sibling_initial_summary;
    BObolDatabaseSourceSummary nested_sibling_summary;
    if (!nested_sibling_source->getSummary(nested_sibling_initial_summary))
	FAIL("GED nested sibling source should expose a state summary");
    if (nested_sibling_initial_summary.realizationStatus !=
	    SoBRLDatabaseSource::REALIZED || nested_sibling_initial_summary.stale)
	FAIL("GED nested sibling source should initially be realized and current");
    if (owned_scene->setDatabaseSourceState(nested_sibling_source_path,
	    TRUE,
	    nested_sibling_initial_summary.sourceRevision,
	    nested_sibling_initial_summary.inputsRevision,
	    nested_sibling_initial_summary.visible,
	    nested_sibling_initial_summary.selected,
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
    const char *prefix_only_group_path = "prefix_owner.c/prefix_only.s";
    if (!ged_draw_obol_group_ensure_for_path(gedp,
	    prefix_only_group_path,
	    prefix_only_group_path,
	    GED_DRAW_MODE_WIRE,
	    0) ||
	    !owned_scene->findGroup(prefix_only_group_path))
	FAIL("GED root path-prefix group-only sentinel should be created");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE_PREFIX,
	    "prefix_owner.c", NULL, -1, "public root path-prefix group-only erase"))
	return 1;
    if (owned_scene->findGroup("prefix_owner.c") ||
	    owned_scene->findGroup(prefix_only_group_path))
	FAIL("GED root path-prefix group-only erase should remove owned Obol groups");
    group_source_state prefix_group_state = {0, GED_DRAW_GROUP_REF_NULL,
	NULL, prefix_only_group_path};
    ged_draw_foreach_group_record(gedp, group_source_state_cb,
	    &prefix_group_state);
    if (prefix_group_state.found)
	FAIL("GED root path-prefix group-only erase should not expose stale group draw records");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE_PREFIX,
	    "nested_parent.c/nested_child.c", ged_draw_active_view_ctx(gedp),
	    -1, "public active-scope path-prefix erase"))
	return 1;
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (source_for_path(owned_scene, nested_leaf_source_path) ||
	    !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED active-scope path-prefix erase should remove matching owned Obol sources only");
    if (owned_scene->findGroup("nested_parent.c/nested_child.c"))
	FAIL("GED active-scope path-prefix erase should remove the owned Obol group subtree");
    record_source_state nested_leaf_record = {0,
	GED_DRAW_SHAPE_REF_NULL_INIT, GED_DRAW_GROUP_REF_NULL_INIT, 0, 0,
	nested_leaf_source_path, 0, 0, 0, 0, 0.0, 0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &nested_leaf_record);
    if (nested_leaf_record.found)
	FAIL("GED active-scope path-prefix erase should not expose stale shape draw records");
    group_source_state nested_child_group_state = {0,
	GED_DRAW_GROUP_REF_NULL, NULL, "nested_parent.c/nested_child.c"};
    ged_draw_foreach_group_record(gedp, group_source_state_cb,
	    &nested_child_group_state);
    if (nested_child_group_state.found)
	FAIL("GED active-scope path-prefix erase should not expose stale group draw records");
    if (ged_exec_draw(gedp, 2, draw_nested_leaf) != BRLCAD_OK)
	FAIL("GED active-scope path-prefix erase redraw restore should succeed");
    nested_leaf_source = source_for_path(owned_scene,
	    nested_leaf_source_path);
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (!nested_leaf_source || !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED active-scope path-prefix erase redraw restore should preserve sibling owned Obol source state");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE_PREFIX,
	    "nested_parent.c/nested_child.c", NULL, -1,
	    "public root path-prefix erase"))
	return 1;
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (source_for_path(owned_scene, nested_leaf_source_path) ||
	    !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED root path-prefix erase should remove matching owned Obol sources only");
    memset(&nested_leaf_record, 0, sizeof(nested_leaf_record));
    nested_leaf_record.matchPath = nested_leaf_source_path;
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &nested_leaf_record);
    if (nested_leaf_record.found)
	FAIL("GED root path-prefix erase should not expose stale shape draw records");
    memset(&nested_child_group_state, 0, sizeof(nested_child_group_state));
    nested_child_group_state.matchPath = "nested_parent.c/nested_child.c";
    ged_draw_foreach_group_record(gedp, group_source_state_cb,
	    &nested_child_group_state);
    if (nested_child_group_state.found)
	FAIL("GED root path-prefix erase should not expose stale group draw records");
    if (ged_exec_draw(gedp, 2, draw_nested_leaf) != BRLCAD_OK)
	FAIL("GED root path-prefix erase redraw restore should succeed");
    nested_leaf_source = source_for_path(owned_scene,
	    nested_leaf_source_path);
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (!nested_leaf_source || !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED root path-prefix erase redraw restore should preserve sibling owned Obol source state");
    if (nested_sibling_summary.realizationStatus !=
	    SoBRLDatabaseSource::REALIZED || nested_sibling_summary.stale)
	FAIL("GED root path-prefix erase redraw restore should keep the sibling realization current");
    struct ged_scene_reducer_request reexpand_nested_child =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_SOURCE_UPDATED,
		"nested_child.c");
    reexpand_nested_child.redraw = 1;
    struct ged_scene_reducer_result nested_result;
    ged_scene_reducer_result_init(&nested_result);
    ged_draw_index_stats_reset(gedp);
    if (ged_scene_reduce(gedp, &reexpand_nested_child,
	    &nested_result) <= 0)
	FAIL("GED nested child reexpand transaction should succeed");
    ged_scene_reducer_result_free(&nested_result);
    struct ged_draw_index_stats nested_reexpand_stats;
    memset(&nested_reexpand_stats, 0, sizeof(nested_reexpand_stats));
    ged_draw_index_stats_get(gedp, &nested_reexpand_stats);
    if (nested_reexpand_stats.slow_path_shape_scans ||
	    nested_reexpand_stats.slow_path_group_scans)
	FAIL("GED nested child reexpand should avoid registry/index slow-path scans");
    struct ged_scene_reducer_request stale_nested_leaf =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_SOURCE_UPDATED,
		"nested_leaf.s");
    stale_nested_leaf.redraw = 0;
    ged_scene_reducer_result_init(&nested_result);
    if (ged_scene_reduce(gedp, &stale_nested_leaf,
	    &nested_result) <= 0)
	FAIL("GED nested leaf stale transaction should succeed");
    ged_scene_reducer_result_free(&nested_result);
    BObolDatabaseSourceSummary nested_leaf_summary;
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
    if (!nested_sibling_source->hasRealizedWireGeometry())
	FAIL("GED nested leaf stale transaction should retain sibling occurrence geometry");
    if (nested_sibling_summary.stale)
	FAIL("GED nested leaf stale transaction should not stale the sibling source");
    if (nested_sibling_summary.lineWidth != 23)
	FAIL("GED nested leaf stale transaction should preserve sibling source state");
    struct ged_scene_reducer_request redraw_nested_leaf =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_SOURCE_UPDATED,
		"nested_leaf.s");
    redraw_nested_leaf.redraw = 1;
    ged_scene_reducer_result_init(&nested_result);
    ged_draw_index_stats_reset(gedp);
    if (ged_scene_reduce(gedp, &redraw_nested_leaf,
	    &nested_result) <= 0)
	FAIL("GED nested leaf redraw transaction should succeed");
    ged_scene_reducer_result_free(&nested_result);
    struct ged_draw_index_stats nested_redraw_stats;
    memset(&nested_redraw_stats, 0, sizeof(nested_redraw_stats));
    ged_draw_index_stats_get(gedp, &nested_redraw_stats);
    if (nested_redraw_stats.slow_path_shape_scans ||
	    nested_redraw_stats.slow_path_group_scans)
	FAIL("GED nested leaf redraw should avoid registry/index slow-path scans");
    nested_leaf_source = source_for_path(owned_scene,
	    nested_leaf_source_path);
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (!nested_leaf_source ||
	    !nested_leaf_source->getSummary(nested_leaf_summary))
	FAIL("GED nested leaf redraw transaction should preserve its owned Obol source");
    if (nested_leaf_summary.stale) {
	fprintf(stderr, "nested leaf redraw remained stale: status=%d reason=0x%x source_rev=%u inputs_rev=%u shapes=%d meshes=%d\n",
	    nested_leaf_summary.realizationStatus,
	    nested_leaf_summary.staleReason,
	    nested_leaf_summary.sourceRevision,
	    nested_leaf_summary.inputsRevision,
	    nested_leaf_summary.realizedShapeCount,
	    nested_leaf_summary.realizedMeshCount);
	for (int source_index = 0;
	     source_index < owned_scene->getDatabaseSourceCount();
	     source_index++) {
	    BObolDatabaseSourceSummary source_summary;
	    if (!owned_scene->getDatabaseSourceSummary(source_index,
		    source_summary) || !source_summary.valid ||
		!path_equal(source_summary.path.getString(),
		    nested_leaf_source_path))
		continue;
	    fprintf(stderr, "nested leaf candidate[%d]: key=%s status=%d stale=%d reason=0x%x source_rev=%u realized_rev=%u shapes=%d meshes=%d\n",
		source_index, source_summary.instanceKey.getString(),
		source_summary.realizationStatus,
		source_summary.stale ? 1 : 0, source_summary.staleReason,
		source_summary.sourceRevision,
		source_summary.realizedSourceRevision,
		source_summary.realizedShapeCount,
		source_summary.realizedMeshCount);
	}
	FAIL("GED nested leaf redraw transaction should realize its owned Obol source");
    }
    if (!nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary))
	FAIL("GED nested leaf redraw transaction should retain the unrelated owned Obol source");
    if (nested_sibling_summary.lineWidth != 23)
	FAIL("GED nested leaf redraw transaction should preserve unrelated owned Obol source presentation");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_SOURCE_REFERENCES_REMOVED,
	    "nested_leaf.s", ged_draw_active_view_ctx(gedp), -1,
	    "public scoped component erase"))
	return 1;
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (source_for_path(owned_scene, nested_leaf_source_path) ||
	    !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED scoped component erase should remove matching owned Obol sources only");
    record_source_state component_leaf_record = {0,
	GED_DRAW_SHAPE_REF_NULL_INIT, GED_DRAW_GROUP_REF_NULL_INIT, 0, 0,
	nested_leaf_source_path, 0, 0, 0, 0, 0.0, 0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &component_leaf_record);
    if (component_leaf_record.found)
	FAIL("GED scoped component erase should not expose stale shape draw records");
    group_source_state component_child_group_state = {0,
	GED_DRAW_GROUP_REF_NULL, NULL, "nested_parent.c/nested_child.c"};
    ged_draw_foreach_group_record(gedp, group_source_state_cb,
	    &component_child_group_state);
    if (component_child_group_state.found)
	FAIL("GED scoped component erase should not expose stale group draw records");
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
    const char *nested_renamed_leaf_source_path =
	"nested_parent.c/nested_child_renamed.c/nested_leaf.s";
    const char *draw_nested_renamed_leaf[2] = {
	"draw", nested_renamed_leaf_source_path
    };
    if (ged_exec_draw(gedp, 2, draw_nested_renamed_leaf) != BRLCAD_OK ||
	    !source_for_path(owned_scene, nested_renamed_leaf_source_path))
	FAIL("GED scoped component mode-filter sentinel draw should succeed");
    if (owned_scene->setDatabaseSourceDrawMode(nested_leaf_source_path,
	    SoBRLDatabaseSource::SHADED) < 0)
	FAIL("GED scoped component mode-filter sentinel should set shaded owner mode");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_SOURCE_REFERENCES_REMOVED,
	    "nested_leaf.s", ged_draw_active_view_ctx(gedp),
	    GED_DRAW_MODE_WIRE, "public scoped component mode-filter erase"))
	return 1;
    nested_leaf_source = source_for_path(owned_scene,
	    nested_leaf_source_path);
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (!nested_leaf_source ||
	    source_for_path(owned_scene, nested_renamed_leaf_source_path) ||
	    !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED scoped component mode-filter erase should remove only matching-mode owned Obol sources");
    BObolDatabaseSourceSummary component_mode_summary;
    if (!nested_leaf_source->getSummary(component_mode_summary) ||
	    component_mode_summary.drawMode != SoBRLDatabaseSource::SHADED)
	FAIL("GED scoped component mode-filter erase should preserve nonmatching owner draw mode");
    memset(&component_leaf_record, 0, sizeof(component_leaf_record));
    component_leaf_record.matchPath = nested_renamed_leaf_source_path;
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &component_leaf_record);
    if (component_leaf_record.found)
	FAIL("GED scoped component mode-filter erase should not expose stale shape draw records");
    if (owned_scene->setDatabaseSourceDrawMode(nested_leaf_source_path,
	    SoBRLDatabaseSource::WIREFRAME) < 0)
	FAIL("GED scoped component mode-filter sentinel should restore wire owner mode");
    struct ged_scene_reducer_request remove_nested_leaf_ref =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_SOURCE_REFERENCES_REMOVED,
		"nested_leaf.s");
    ged_scene_reducer_result_init(&nested_result);
    if (ged_scene_reduce(gedp, &remove_nested_leaf_ref,
	    &nested_result) <= 0)
	FAIL("GED nested leaf reference removal transaction should succeed");
    ged_scene_reducer_result_free(&nested_result);
    nested_sibling_source = source_for_path(owned_scene,
	    nested_sibling_source_path);
    if (source_for_path(owned_scene, nested_leaf_source_path) ||
	    !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED nested leaf reference removal should remove only non-root owned Obol sources");
    struct ged_scene_reducer_request remove_nested_sibling =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_SOURCE_UPDATED,
		"nested_sibling.s");
    remove_nested_sibling.removed = 1;
    ged_scene_reducer_result_init(&nested_result);
    if (ged_scene_reduce(gedp, &remove_nested_sibling,
	    &nested_result) <= 0)
	FAIL("GED nested sibling source removal transaction should succeed");
    ged_scene_reducer_result_free(&nested_result);
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
    ged_scene_node_ref nested_child_ctx = ged_scene_group_node(gedp,
	    nested_child_group);
    struct db_full_path nested_child_ctx_path;
    db_full_path_init(&nested_child_ctx_path);
    if (ged_scene_node_ref_is_null(nested_child_ctx) ||
	    ged_scene_group_dbpath(gedp, nested_child_ctx,
		&nested_child_ctx_path) != 0)
	FAIL("GED nested group-ref context should expose owned Obol DB path");
    char *nested_child_ctx_path_str =
	db_path_to_string(&nested_child_ctx_path);
    if (!nested_child_ctx_path_str ||
	    !path_equal(nested_child_ctx_path_str,
		"nested_parent.c/nested_child.c"))
	FAIL("GED nested group-ref context should preserve owned Obol DB path");
    if (nested_child_ctx_path_str)
	bu_free(nested_child_ctx_path_str, "db_path_to_string");
    db_free_full_path(&nested_child_ctx_path);
    struct db_full_path nested_child_renamed_path;
    db_full_path_init(&nested_child_renamed_path);
    if (db_string_to_path(&nested_child_renamed_path, gedp->dbip,
	    "nested_parent.c/nested_child_renamed.c") != 0)
	FAIL("GED nested rename sentinel path should resolve");
    ged_draw_index_stats_reset(gedp);
    if (!ged_draw_group_ref_set_dbpath(gedp, nested_child_group,
	    &nested_child_renamed_path))
	FAIL("GED nested group set-dbpath should rename through owned Obol");
    ged_draw_group_ref nested_child_renamed_group =
	ged_draw_group_ref_lookup_or_create(gedp, &nested_child_renamed_path);
    db_free_full_path(&nested_child_renamed_path);
    if (ged_draw_group_ref_is_null(nested_child_renamed_group))
	FAIL("GED nested renamed group ref should resolve after owned Obol rename");
    if (owned_scene->findGroup("nested_parent.c/nested_child.c") ||
	    !owned_scene->findGroup(
		"nested_parent.c/nested_child_renamed.c"))
	FAIL("GED nested group set-dbpath should rename the owned Obol group in place");
    struct ged_draw_group_record nested_child_record;
    memset(&nested_child_record, 0, sizeof(nested_child_record));
    if (!ged_draw_group_record_get(gedp, nested_child_renamed_group,
	    &nested_child_record) ||
	    !nested_child_record.path ||
	    !path_equal(nested_child_record.path,
		"nested_parent.c/nested_child_renamed.c"))
	FAIL("GED nested group set-dbpath should expose the owned Obol renamed path");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE,
	    "nested_parent.c/nested_child_renamed.c",
	    ged_draw_active_view_ctx(gedp), -1,
	    "public active-scope nested group path erase"))
	return 1;
    if (!owned_scene->findGroup("nested_parent.c") ||
	    owned_scene->findGroup("nested_parent.c/nested_child_renamed.c"))
		FAIL("GED public active-scope nested group path erase should remove the owned Obol nested group");
    if (apply_path_transaction(gedp, GED_SCENE_REDUCER_ERASE, "nested_parent.c",
	    ged_draw_active_view_ctx(gedp), -1,
	    "public active-scope nested group path cleanup"))
	return 1;
    if (owned_scene->findGroup("nested_parent.c") ||
	    owned_scene->getDatabaseSourceCount() != 2)
	FAIL("GED public active-scope nested group path cleanup should restore owned Obol baseline");
    std::string clear_child_path = group_path + "/__obol_clear_child";
    if (!owned_scene->ensureGroup(clear_child_path.c_str()) ||
	    !owned_scene->findGroup(clear_child_path.c_str()))
	FAIL("owned Obol direct clear sentinel child group should be created");

    struct occurrence_visit_state occurrence_state = {0, 1};
    size_t occurrence_count = ged_scene_occurrence_count(gedp);
    if (!occurrence_count || !ged_scene_has_occurrences(gedp) ||
	ged_scene_occurrences_visit(gedp, occurrence_visit_cb,
	    &occurrence_state) != occurrence_count ||
	occurrence_state.count != occurrence_count || !occurrence_state.valid)
	FAIL("public occurrence enumeration should expose valid realized objects");
    ged_scene_occurrence_ref retired_occurrence =
	ged_scene_occurrence_first(gedp);
    struct ged_scene_occurrence_info retired_info;
    memset(&retired_info, 0, sizeof(retired_info));
    if (ged_scene_occurrence_ref_is_null(retired_occurrence) ||
	!ged_scene_occurrence_get(gedp, retired_occurrence, &retired_info) ||
	!retired_info.path || !retired_info.path[0])
	FAIL("public occurrence lookup should resolve the first realized object");
    if (ged_scene_occurrence_highlight_set(gedp, retired_occurrence, 1,
	    NULL) != GED_SCENE_OK ||
	!ged_scene_occurrence_get(gedp, retired_occurrence, &retired_info) ||
	!retired_info.highlighted)
	FAIL("public occurrence identity should survive a presentation mutation");
    if (ged_scene_occurrence_highlight_set(gedp, retired_occurrence, 0,
	    NULL) != GED_SCENE_OK)
	FAIL("public occurrence highlight cleanup should succeed");

    ged_draw_clear(gedp);
    if (owned_scene->getDatabaseSourceCount() != 0 ||
	    owned_scene->findGroup(clear_child_path.c_str()) ||
	    owned_scene->findGroup(group_path.c_str()))
	FAIL("GED direct draw clear should remove owned Obol group/source subtrees");
    record_source_state clear_box_record = {0, GED_DRAW_SHAPE_REF_NULL_INIT,
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "box.s", 0, 0, 0, 0, 0.0,
	0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &clear_box_record);
    if (clear_box_record.found)
	FAIL("GED direct draw clear should not expose stale shape draw records");
    group_source_state clear_box_group_state = {0,
	GED_DRAW_GROUP_REF_NULL, NULL, group_path.c_str()};
    ged_draw_foreach_group_record(gedp, group_source_state_cb,
	    &clear_box_group_state);
    if (clear_box_group_state.found)
	FAIL("GED direct draw clear should not expose stale group draw records");
    memset(&retired_info, 0, sizeof(retired_info));
    if (ged_scene_occurrence_get(gedp, retired_occurrence, &retired_info))
	FAIL("a retired public occurrence reference must not resolve after clear");
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller direct-clear redraw commands should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after direct clear");
    ged_scene_occurrence_ref replacement_occurrence =
	ged_scene_occurrence_first(gedp);
    if (ged_scene_occurrence_ref_is_null(replacement_occurrence) ||
	ged_scene_occurrence_ref_equal(retired_occurrence,
	    replacement_occurrence) ||
	!ged_scene_occurrence_get(gedp, replacement_occurrence, &retired_info))
	FAIL("redraw must allocate a valid occurrence identity that cannot alias a retired reference");
    std::string scoped_clear_child_path =
	group_path + "/__obol_scoped_clear_child";
    if (!owned_scene->ensureGroup(scoped_clear_child_path.c_str()) ||
	    !owned_scene->findGroup(scoped_clear_child_path.c_str()))
	FAIL("owned Obol scoped clear sentinel child group should be created");
    if (!ged_draw_clear_view(gedp, ged_draw_active_view_ctx(gedp)))
	FAIL("GED public scoped database-group clear should succeed");
    if (owned_scene->getDatabaseSourceCount() != 0 ||
	    owned_scene->findGroup(scoped_clear_child_path.c_str()) ||
	    owned_scene->findGroup(group_path.c_str()))
	FAIL("GED scoped database-group clear should remove owned Obol group/source subtrees");
    memset(&clear_box_record, 0, sizeof(clear_box_record));
    clear_box_record.matchPath = "box.s";
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &clear_box_record);
    if (clear_box_record.found)
	FAIL("GED scoped database-group clear should not expose stale shape draw records");
    memset(&clear_box_group_state, 0, sizeof(clear_box_group_state));
    clear_box_group_state.matchPath = group_path.c_str();
    ged_draw_foreach_group_record(gedp, group_source_state_cb,
	    &clear_box_group_state);
    if (clear_box_group_state.found)
	FAIL("GED scoped database-group clear should not expose stale group draw records");
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller scoped-clear redraw commands should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after scoped clear");
    if (exercise_multi_instance_transform_reuse(gedp, owned_scene))
	return 1;

    SoSeparator *root = new SoSeparator;
    root->ref();
    BObolViewController view_controller(root);
    bobol_display_endpoint_t *first_populated_endpoint =
	bobol_display_endpoint_create(&view_controller, 0);
    if (!first_populated_endpoint ||
	!ged_view_context_obol_endpoint_set(initial_view_ctx,
	    first_populated_endpoint, 1)) {
	if (first_populated_endpoint)
	    bobol_display_endpoint_destroy(first_populated_endpoint);
	FAIL("GED populated view should accept its first display endpoint");
    }
    SoNode *retained_render_root = view_controller.getRenderSceneRoot();
    bobol_display_endpoint_t *replacement_endpoint =
	bobol_display_endpoint_create(&view_controller, 0);
    if (!replacement_endpoint ||
	!ged_view_context_obol_endpoint_set(initial_view_ctx,
	    replacement_endpoint, 1)) {
	if (replacement_endpoint)
	    bobol_display_endpoint_destroy(replacement_endpoint);
	FAIL("GED populated view should accept a replacement display endpoint");
    }
    if (ged_view_context_obol_endpoint_get(initial_view_ctx) !=
	    replacement_endpoint ||
	ged_draw_obol_scene_controller(gedp) != owned_scene ||
	!view_controller.getRenderSceneRoot() ||
	view_controller.getRenderSceneRoot() != retained_render_root)
	FAIL("endpoint replacement should retain and rebind the shared scene");
    BObolSceneController *view_scene = owned_scene;
    if (view_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(view_scene, "box.s") ||
	    !source_for_path(view_scene, "ball.s"))
	FAIL("endpoint replacement should preserve populated shared draw state");

    BObolViewController progressive_controller(new SoBRLSceneGroup);
    struct ged_view_context *progressive_view_ctx = ged_view_context_create();
    bobol_display_endpoint_t *progressive_endpoint =
	bobol_display_endpoint_create(&progressive_controller, 0);
    if (!progressive_view_ctx || !progressive_endpoint ||
	!ged_view_context_host_attach(gedp, progressive_view_ctx) ||
	!ged_view_context_obol_endpoint_set(progressive_view_ctx,
	    progressive_endpoint, 1)) {
	if (progressive_endpoint && (!progressive_view_ctx ||
	    ged_view_context_obol_endpoint_get(progressive_view_ctx) !=
		progressive_endpoint))
	    bobol_display_endpoint_destroy(progressive_endpoint);
	if (progressive_view_ctx)
	    ged_view_context_free(progressive_view_ctx);
	FAIL("second endpoint should register per-view progressive services");
    }
    if (exercise_progressive_autoview_lifecycle(gedp,
	    &progressive_controller,
	    progressive_view_ctx, GED_DRAW_MODE_WIRE,
	    "progressive_root.c", "progressive_root_async.c", 4, false))
	return 1;
    (void)ged_view_context_obol_endpoint_set(progressive_view_ctx,
	NULL, 0);
    ged_view_context_free(progressive_view_ctx);

    if (!seed_view_lod_probe_payload(&view_controller, "box.s", "box.s"))
	FAIL("attached Obol view-controller LoD invalidation probe should seed draw payload");
    const char *attached_draw_draft[2] = {"draw", "draft_move.s"};
    const int attached_draw_ret = ged_exec_draw(gedp, 2, attached_draw_draft);
    if (attached_draw_ret != BRLCAD_OK ||
	    view_controller.getViewLodState()->payloadCount() != 0) {
	fprintf(stderr, "attached draw ret=%d lod_payloads=%zu\n",
	    attached_draw_ret,
	    view_controller.getViewLodState()->payloadCount());
	FAIL("attached Obol draw transaction should clear view-local LoD state");
    }
    if (view_scene->getDatabaseSourceCount() != 3 ||
	    !source_for_path(view_scene, "draft_move.s"))
	FAIL("attached Obol draw transaction should still sync the drawn source");

    if (!seed_view_lod_probe_payload(&view_controller, "box.s", "box.s"))
	FAIL("attached Obol view-controller LoD invalidation probe should seed erase payload");
    const char *attached_erase_draft[2] = {"erase", "draft_move.s"};
    if (ged_exec_erase(gedp, 2, attached_erase_draft) != BRLCAD_OK ||
	    view_controller.getViewLodState()->payloadCount() != 0)
	FAIL("attached Obol erase transaction should clear view-local LoD state");
    if (view_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(view_scene, "draft_move.s"))
	FAIL("attached Obol erase transaction should remove the erased source");

    if (!seed_view_lod_probe_payload(&view_controller, "box.s", "box.s"))
	FAIL("attached Obol view-controller LoD invalidation probe should seed source-update payload");
    struct ged_scene_reducer_request attached_source_update =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_SOURCE_UPDATED, "box.s");
    attached_source_update.redraw = 0;
    struct ged_scene_reducer_result attached_source_result;
    ged_scene_reducer_result_init(&attached_source_result);
    if (ged_scene_reduce(gedp, &attached_source_update,
	    &attached_source_result) <= 0 ||
	    view_controller.getViewLodState()->payloadCount() != 0) {
	ged_scene_reducer_result_free(&attached_source_result);
	FAIL("attached Obol source-update transaction should clear view-local LoD state");
    }
    ged_scene_reducer_result_free(&attached_source_result);
    const char *attached_redraw_box[2] = {"draw", "box.s"};
    if (ged_exec_draw(gedp, 2, attached_redraw_box) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(view_scene, "box.s"))
	FAIL("attached Obol source-update refresh should restore box source");

    if (!seed_view_lod_probe_payload(&view_controller, "box.s", "box.s"))
	FAIL("attached Obol view-controller LoD invalidation probe should seed clear payload");
    struct ged_scene_reducer_request attached_clear =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_CLEAR, NULL);
    struct ged_scene_reducer_result attached_clear_result;
    ged_scene_reducer_result_init(&attached_clear_result);
    if (ged_scene_reduce(gedp, &attached_clear,
	    &attached_clear_result) <= 0) {
	ged_scene_reducer_result_free(&attached_clear_result);
	FAIL("attached Obol clear transaction should succeed");
    }
    ged_scene_reducer_result_free(&attached_clear_result);
    if (view_controller.getViewLodState()->payloadCount() != 0 ||
	    view_scene->getDatabaseSourceCount() != 0)
	FAIL("attached Obol clear transaction should clear view-local LoD state and scene sources");
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(view_scene, "box.s") ||
	    !source_for_path(view_scene, "ball.s"))
	FAIL("attached Obol clear transaction redraw should restore current GED draw state");

    struct ged_scene_reducer_request attached_stale_source =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_STALE_SOURCE, "box.s");
    attached_stale_source.stale_reason = GED_DRAW_STALE_SETTINGS_CHANGED;
    if (apply_attached_view_lod_invalidation_probe(gedp, &view_controller,
	    &attached_stale_source, "stale-source"))
	return 1;
    if (view_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(view_scene, "box.s") ||
	    !source_for_path(view_scene, "ball.s"))
	FAIL("attached Obol stale-source transaction should preserve scene sources");

    if (ged_exec_draw(gedp, 2, attached_draw_draft) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 3 ||
	    !source_for_path(view_scene, "draft_move.s"))
	FAIL("attached Obol erase-prefix setup should draw draft source");
    struct ged_scene_reducer_request attached_erase_prefix =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_ERASE_PREFIX,
		"draft_move.s");
    if (apply_attached_view_lod_invalidation_probe(gedp, &view_controller,
	    &attached_erase_prefix, "erase-prefix"))
	return 1;
    if (view_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(view_scene, "draft_move.s"))
	FAIL("attached Obol erase-prefix transaction should remove matching source");

    const char *attached_draw_nested_leaf[2] = {
	"draw", "nested_parent.c/nested_child.c/nested_leaf.s"
    };
    if (ged_exec_draw(gedp, 2, attached_draw_nested_leaf) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 3 ||
	    !source_for_path(view_scene,
		"nested_parent.c/nested_child.c/nested_leaf.s"))
	FAIL("attached Obol reference-removal setup should draw nested source");
    struct ged_scene_reducer_request attached_reference_removal =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_SOURCE_REFERENCES_REMOVED,
		"nested_leaf.s");
    if (apply_attached_view_lod_invalidation_probe(gedp, &view_controller,
	    &attached_reference_removal, "source-reference removal"))
	return 1;
    if (view_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(view_scene,
		"nested_parent.c/nested_child.c/nested_leaf.s"))
	FAIL("attached Obol source-reference removal should remove matching nested source");

    const char *attached_draw_renamed_source[2] = {
	"draw", "renamed_source.s"
    };
    if (ged_exec_draw(gedp, 2, attached_draw_renamed_source) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 3 ||
	    !source_for_path(view_scene, "renamed_source.s"))
	FAIL("attached Obol rename setup should draw renamed source");
    if (!seed_view_lod_probe_payload(&view_controller, "box.s", "box.s"))
	FAIL("attached Obol view-controller LoD invalidation probe should seed source-rename payload");
    const char *attached_move_source[3] = {
	"move", "renamed_source.s", "__obol_attached_renamed_source.s"
    };
    if (ged_exec(gedp, 3, attached_move_source) != BRLCAD_OK ||
	    view_controller.getViewLodState()->payloadCount() != 0)
	FAIL("attached Obol source-rename command should clear view-local LoD state");
    if (view_scene->getDatabaseSourceCount() != 3 ||
	    source_for_path(view_scene, "renamed_source.s") ||
	    !source_for_path(view_scene, "__obol_attached_renamed_source.s"))
	FAIL("attached Obol source-rename command should rename source in place");
    const char *attached_move_source_back[3] = {
	"move", "__obol_attached_renamed_source.s", "renamed_source.s"
    };
    if (ged_exec(gedp, 3, attached_move_source_back) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 3 ||
	    !source_for_path(view_scene, "renamed_source.s") ||
	    source_for_path(view_scene, "__obol_attached_renamed_source.s"))
	FAIL("attached Obol source-rename restore should return to database source name");
    const char *attached_erase_renamed_source[2] = {
	"erase", "renamed_source.s"
    };
    if (ged_exec_erase(gedp, 2, attached_erase_renamed_source) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(view_scene, "renamed_source.s"))
	FAIL("attached Obol source-rename cleanup should restore baseline sources");

    struct ged_scene_reducer_request attached_clear_scope =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_CLEAR_SCOPE, NULL);
    attached_clear_scope.view = ged_draw_active_view_ctx(gedp);
    if (apply_attached_view_lod_invalidation_probe(gedp, &view_controller,
	    &attached_clear_scope, "clear-scope"))
	return 1;
    if (view_scene->getDatabaseSourceCount() != 0)
	FAIL("attached Obol clear-scope transaction should clear scoped scene sources");
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(view_scene, "box.s") ||
	    !source_for_path(view_scene, "ball.s"))
	FAIL("attached Obol clear-scope transaction redraw should restore baseline sources");

    struct ged_scene_reducer_request attached_teardown =
	ged_scene_reducer_request_make(GED_SCENE_REDUCER_TEARDOWN, NULL);
    if (apply_attached_view_lod_invalidation_probe(gedp, &view_controller,
	    &attached_teardown, "teardown"))
	return 1;
    if (view_scene->getDatabaseSourceCount() != 0)
	FAIL("attached Obol teardown transaction should clear scene sources");
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(view_scene, "box.s") ||
	    !source_for_path(view_scene, "ball.s"))
	FAIL("attached Obol teardown transaction redraw should restore baseline sources");

    struct ged_view_context *independent_view = ged_view_context_create();
    if (!independent_view)
	FAIL("Obol independent view source-owner test view should be created");
    if (!bv_name_set(DRAW_TEST_BV(independent_view), "V0"))
	FAIL("Obol independent view source-owner test view should be named");
    if (!ged_view_set_context_add(ged_view_set_ctx(gedp), independent_view))
	FAIL("Obol independent view source-owner test view should be registered");
    ged_view_context_owned_add(gedp, independent_view);

    BObolViewController independent_controller(new SoBRLSceneGroup);
    bobol_display_endpoint_t *independent_endpoint =
	bobol_display_endpoint_create(&independent_controller, 0);
    if (!independent_endpoint ||
	!ged_view_context_host_attach(gedp, independent_view) ||
	!ged_view_context_obol_endpoint_set(independent_view,
	    independent_endpoint, 1)) {
	if (independent_endpoint &&
	    ged_view_context_obol_endpoint_get(independent_view) !=
		independent_endpoint)
	    bobol_display_endpoint_destroy(independent_endpoint);
	FAIL("Obol independent view should accept a local display endpoint");
    }
    BObolSceneController *independent_scene =
	independent_controller.getSceneController();
    if (!independent_scene)
	FAIL("Obol independent view should own a local scene controller");

    const char *view_independent_on[5] = {"view", "independent", "V0",
	"1", NULL};
    if (ged_exec_view(gedp, 4, view_independent_on) != BRLCAD_OK)
	FAIL("real GED view independent command should succeed");

    const char *draw_v0_box[6] = {"draw", "-R", "-V", "V0", "box.s",
	NULL};
    const char *draw_v0_ball[6] = {"draw", "-R", "-V", "V0", "ball.s",
	NULL};
    if (ged_exec_draw(gedp, 5, draw_v0_box) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 5, draw_v0_ball) != BRLCAD_OK)
	FAIL("real GED independent-view draw commands should succeed");

    if (view_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_shared_path(view_scene, "box.s") ||
	    !source_for_shared_path(view_scene, "ball.s") ||
	    independent_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_view_path(independent_scene, "V0", "box.s") ||
	    !source_for_view_path(independent_scene, "V0", "ball.s"))
	FAIL("Obol independent view setup should isolate local source owners from shared owners");

    const char *erase_v0_box[5] = {"erase", "-V", "V0", "box.s", NULL};
    if (ged_exec_erase(gedp, 4, erase_v0_box) != BRLCAD_OK)
	FAIL("real GED independent-view erase command should succeed");

    if (view_scene->getDatabaseSourceCount() != 2 ||
	    independent_scene->getDatabaseSourceCount() != 1 ||
	    source_for_view_path(independent_scene, "V0", "box.s") ||
	    !source_for_view_path(independent_scene, "V0", "ball.s") ||
	    !source_for_shared_path(view_scene, "box.s") ||
	    !source_for_shared_path(view_scene, "ball.s"))
	FAIL("Obol independent-view erase should remove only the scoped source owner");

    if (ged_exec_draw(gedp, 5, draw_v0_box) != BRLCAD_OK)
	FAIL("real GED independent-view redraw command should succeed");
    if (view_scene->getDatabaseSourceCount() != 2 ||
	    independent_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_view_path(independent_scene, "V0", "box.s") ||
	    !source_for_view_path(independent_scene, "V0", "ball.s") ||
	    !source_for_shared_path(view_scene, "box.s") ||
	    !source_for_shared_path(view_scene, "ball.s"))
	FAIL("Obol independent-view draw should restore only the scoped source owner");

    const char *zap_v0[5] = {"zap", "-V", "V0", "-g", NULL};
    if (ged_exec_zap(gedp, 4, zap_v0) != BRLCAD_OK)
	FAIL("real GED independent-view zap command should succeed");
    if (view_scene->getDatabaseSourceCount() != 2 ||
	    independent_scene->getDatabaseSourceCount() != 0 ||
	    source_for_view_path(independent_scene, "V0", "box.s") ||
	    source_for_view_path(independent_scene, "V0", "ball.s") ||
	    !source_for_shared_path(view_scene, "box.s") ||
	    !source_for_shared_path(view_scene, "ball.s"))
	FAIL("Obol independent-view zap should clear only scoped source owners");

    const char *view_independent_off[5] = {"view", "independent", "V0",
	"0", NULL};
    if (ged_exec_view(gedp, 4, view_independent_off) != BRLCAD_OK ||
	    ged_view_context_is_independent(independent_view))
	FAIL("real GED view independent-off command should restore shared view semantics");

    const char *erase_box[2] = {"erase", "box.s"};
    if (ged_exec_erase(gedp, 2, erase_box) != BRLCAD_OK)
	FAIL("real GED erase command should succeed");
    if (view_scene->getDatabaseSourceCount() != 1 ||
	    source_for_path(view_scene, "box.s") ||
	    !source_for_path(view_scene, "ball.s"))
	FAIL("attached Obol scene controller should mirror GED erase command");

    point_t zap_polygon_origin = {0.0, 0.0, 0.0};
    struct ged_view_context *zap_view_ctx = ged_view_active_ctx(gedp);
    ged_view_polygon_ref zap_created =
	ged_view_polygon_create(zap_view_ctx,
		"zap::polygon", 0, GED_VIEW_POLYGON_SQUARE,
		zap_polygon_origin);
    ged_view_polygon_ref zap_found =
	ged_view_polygon_find(zap_view_ctx, "zap::polygon");
    if (ged_view_polygon_ref_is_null(zap_created) ||
	    ged_view_polygon_ref_is_null(zap_found))
	FAIL("Obol polygon store should publish a test polygon before zap");

    const char *zap_cmd[1] = {"zap"};
    if (ged_exec_zap(gedp, 1, zap_cmd) != BRLCAD_OK)
	FAIL("real GED zap command should succeed");
    if (view_scene->getDatabaseSourceCount() != 0)
	FAIL("attached Obol scene controller should mirror GED clear command");
    if (!ged_view_polygon_ref_is_null(
		ged_view_polygon_find(zap_view_ctx,
		    "zap::polygon")))
	FAIL("GED zap should clear Obol polygon store view features");

    if (exercise_ged_value_handle_lifetimes(gedp))
	return 1;

    (void)ged_view_context_obol_endpoint_set(initial_view_ctx, NULL, 0);
    ged_close(gedp);
    root->unref();
    bu_file_delete(dbpath);
    bu_dirclear(lcache);
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
