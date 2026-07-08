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
#include "brlobol/pick_detail.h"
#include "brlobol/scene_controller.h"
#include "brlobol/scene_group.h"
#include "brlobol/vlist_shape.h"
#include "brlobol/view_controller.h"
#include "brlobol/view_lod.h"
#include "brlobol/view_store.h"
#include "bg/line_layer.h"
#include "bg/plot3.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/str.h"
#include "ged.h"
#include "ged/commands.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "opennurbs_sphere.h"
#include "rt/db_internal.h"
#include "rt/view.h"
#include "view_test_util.h"
#include "wdb.h"

#include "../../ged_private.h"

#include <Inventor/SbMatrix.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/nodes/SoGroup.h>
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
    fastf_t mesh_owner_vertices[12] = {
	0.0, 0.0, 0.0,
	2.0, 0.0, 0.0,
	0.0, 2.0, 0.0,
	2.0, 2.0, 0.0
    };
    int mesh_owner_faces[6] = {
	0, 1, 2,
	1, 3, 2
    };
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
		4, 2, mesh_owner_vertices, mesh_owner_faces, NULL, NULL) == 0 &&
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

struct rt_preview_callback_state {
    uint64_t revision;
    int update_count;
    int pick_count;
};

static uint64_t
rt_preview_revision_cb(void *data)
{
    struct rt_preview_callback_state *ctx =
	(struct rt_preview_callback_state *)data;
    return ctx ? ctx->revision : 0;
}

static int
rt_preview_update_cb(void *data)
{
    struct rt_preview_callback_state *ctx =
	(struct rt_preview_callback_state *)data;
    if (!ctx)
	return 0;
    ctx->update_count++;
    return 1;
}

static int
rt_preview_pick_cb(void *data, int UNUSED(x), int UNUSED(y),
	void *UNUSED(pick_out))
{
    struct rt_preview_callback_state *ctx =
	(struct rt_preview_callback_state *)data;
    if (!ctx)
	return 0;
    ctx->pick_count++;
    return 1;
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
    int saw_custom_update;
    uint64_t line_layers_feature_id;
    uint64_t metadata_feature_id;
    uint64_t primitive_metadata_feature_id;
    uint64_t custom_feature_id;
};

struct custom_node_provider_state {
    int call_count;
    int saw_request;
    uint64_t generation;
    int local;
    SoNode *node;
};

static void
command_result_cb(const struct ged_draw_command_result *result, void *data)
{
    struct command_result_callback_state *ctx =
	(struct command_result_callback_state *)data;
    if (!ctx || !result)
	return;

    ctx->callback_count++;
    switch (result->status) {
	case GED_DRAW_COMMAND_RESULT_ACCEPTED:
	    ctx->accepted_count++;
	    break;
	case GED_DRAW_COMMAND_RESULT_UPDATED:
	    ctx->updated_count++;
	    break;
	case GED_DRAW_COMMAND_RESULT_REMOVED:
	    ctx->removed_count++;
	    break;
	case GED_DRAW_COMMAND_RESULT_FAILED:
	    ctx->failed_count++;
	    break;
	default:
	    break;
    }

    const char *name = result->feature_name ? result->feature_name : "";
    const char *command = result->command ? result->command : "";
    if (result->status == GED_DRAW_COMMAND_RESULT_UPDATED &&
	    BU_STR_EQUAL(name, "rtcheck::overlaps") &&
	    BU_STR_EQUAL(command, "lineLayersReplace")) {
	ctx->saw_line_layers_update = 1;
	ctx->line_layers_feature_id = result->feature_id;
    }
    if (result->status == GED_DRAW_COMMAND_RESULT_UPDATED &&
	    BU_STR_EQUAL(name, "rtcheck::overlaps") &&
	    BU_STR_EQUAL(command, "metadataReplace")) {
	ctx->saw_metadata_update = 1;
	ctx->metadata_feature_id = result->feature_id;
    }
    if (result->status == GED_DRAW_COMMAND_RESULT_UPDATED &&
	    BU_STR_EQUAL(name, "rtcheck::overlaps") &&
	    BU_STR_EQUAL(command, "primitiveMetadataReplace")) {
	ctx->saw_primitive_metadata_update = 1;
	ctx->primitive_metadata_feature_id = result->feature_id;
    }
    if (result->status == GED_DRAW_COMMAND_RESULT_REMOVED &&
	    BU_STR_EQUAL(command, "removePrefix") &&
	    strncmp(name, "rtcheck::", strlen("rtcheck::")) == 0)
	ctx->saw_remove_prefix = 1;
    if (result->status == GED_DRAW_COMMAND_RESULT_FAILED &&
	    BU_STR_EQUAL(name, "rtcheck::generation") &&
	    BU_STR_EQUAL(command, "lineLayersReplace"))
	ctx->saw_stale_failure = 1;
    if (result->status == GED_DRAW_COMMAND_RESULT_FAILED &&
	    BU_STR_EQUAL(command, "commit"))
	ctx->saw_commit_failure = 1;
    if (result->status == GED_DRAW_COMMAND_RESULT_UPDATED &&
	    BU_STR_EQUAL(name, "custom::node") &&
	    BU_STR_EQUAL(command, "customNodeReplace")) {
	ctx->saw_custom_update = 1;
	ctx->custom_feature_id = result->feature_id;
    }
}

static void *
custom_node_provider_cb(
	const struct ged_draw_command_scene_custom_node_request *request,
	void *data)
{
    struct custom_node_provider_state *ctx =
	(struct custom_node_provider_state *)data;
    if (!ctx || !request)
	return NULL;

    ctx->call_count++;
    ctx->saw_request =
	BU_STR_EQUAL(request->feature_name, "custom::node") &&
	BU_STR_EQUAL(request->owner_id, "custom") &&
	BU_STR_EQUAL(request->owner_role, "command-result");
    ctx->generation = request->generation;
    ctx->local = request->local;
    ctx->node = new SoSeparator;
    return ctx->node;
}

static int
feature_overlay_matches(BRLObolViewController *controller,
	const char *name,
	BRLObolOverlayClass overlay_class,
	BRLObolOverlayLifecycle lifecycle,
	BRLObolOverlayOrder order)
{
    if (!controller || !name)
	return 0;

    BRLObolFeatureHandle handle = controller->features().find(name);
    BRLObolOverlayInfo overlay;
    if (!handle.isValid() ||
	    !controller->features().overlayInfo(handle, overlay))
	return 0;

    return overlay.isOverlay &&
	overlay.role == BRLObolOverlayRole::Model &&
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

static int
seed_view_lod_probe_payload(BRLObolViewController *controller,
			    const char *path,
			    const char *name)
{
    if (!controller || !controller->getViewLodState() || !path || !name)
	return 0;

    SoBRLMeshShape *mesh = new SoBRLMeshShape;
    mesh->ref();
    mesh->sourcePath = path;
    mesh->sourceName = name;

    BRLObolLodRequest request;
    request.databaseId = "ged-obol-sync";
    request.sourceRevision = 1;
    request.sourceContentHash = 1;
    request.objectPath = path;
    request.objectName = name;
    request.viewRevision = 1;
    request.policyRevision = 1;
    request.drawMode = BRLOBOL_LOD_DRAW_WIRE;
    request.providerId = "ged-obol-sync-probe";
    request.providerVersion = "1";
    request.qualityTier = BRLOBOL_LOD_QUALITY_PROXY;
    request.bounds = SbBox3f(SbVec3f(-1.0f, -1.0f, -1.0f),
	    SbVec3f(1.0f, 1.0f, 1.0f));

    BRLObolLodCounts counts;
    counts.faceCount = 1;
    counts.pointCount = 2;
    BRLObolLodResult result =
	brlobol_lod_aabb_result(request, request.bounds, &counts);

    const int seeded = controller->getViewLodState()->applyProxyResult(
	    mesh, result) &&
	controller->getViewLodState()->payloadCount() > 0;
    mesh->unref();
    return seeded;
}

static int
apply_attached_view_lod_invalidation_probe(struct ged *gedp,
	BRLObolViewController *controller,
	struct ged_draw_transaction *txn,
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

    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int ret = ged_draw_apply_transaction(gedp, txn, &result);
    ged_draw_transaction_result_free(&result);
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
    return instance_key && strncmp(instance_key, prefix.c_str(),
	    prefix.length()) == 0;
}

static int
source_instance_is_any_view_scoped(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;
    const char *instance_key = source->instanceKey.getValue().getString();
    return instance_key && strncmp(instance_key, "ged-view:", 9) == 0;
}

static SoBRLDatabaseSource *
source_for_view_path(SoBRLSceneController *controller,
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
source_for_shared_path(SoBRLSceneController *controller, const char *path)
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
source_for_representation(SoBRLSceneController *controller,
	const char *path,
	int representation_mode)
{
    if (!controller || !path)
	return NULL;

    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (!source)
	    continue;
	BRLObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    continue;
	if (path_equal(summary.path.getString(), path) &&
		summary.representationMode == representation_mode)
	    return source;
    }

    return NULL;
}

static int
source_representation_count(SoBRLSceneController *controller,
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
	BRLObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    continue;
	if (path_equal(summary.path.getString(), path) &&
		summary.representationMode == representation_mode)
	    count++;
    }

    return count;
}

static int
source_path_count(SoBRLSceneController *controller,
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
verify_mode_source(SoBRLSceneController *controller,
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

    BRLObolDatabaseSourceSummary summary;
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
    if (expect_vlist && summary.realizedShapeCount <= 0)
	FAIL("mode-specific source should carry realized VLIST geometry");
    if (expect_mesh && summary.realizedMeshCount <= 0)
	FAIL("mode-specific source should carry realized mesh geometry");
    if (!summary.sourceBoundsValid || summary.sourceBounds.isEmpty())
	FAIL("mode-specific source should retain valid bounds");

    (void)label;
    return 0;
}

static int
apply_mode_value_transaction(struct ged *gedp,
	ged_draw_transaction_kind kind,
	const char *path,
	int mode,
	fastf_t value,
	const char *label)
{
    struct ged_draw_transaction txn =
	ged_draw_transaction_make_value(kind, path, value);
    txn.mode = mode;
    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int ret = ged_draw_apply_transaction(gedp, &txn, &result);
    ged_draw_transaction_result_free(&result);
    if (ret <= 0)
	FAIL("mode-specific value transaction should succeed");

    (void)label;
    return 0;
}

static int
apply_mode_path_transaction(struct ged *gedp,
	ged_draw_transaction_kind kind,
	const char *path,
	int mode,
	const char *label)
{
    struct ged_draw_transaction txn = ged_draw_transaction_make(kind, path);
    txn.mode = mode;
    if (kind == GED_DRAW_TXN_STALE_SOURCE)
	txn.stale_reason = GED_DRAW_STALE_SETTINGS_CHANGED;
    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int ret = ged_draw_apply_transaction(gedp, &txn, &result);
    ged_draw_transaction_result_free(&result);
    if (ret <= 0)
	FAIL("mode-specific path transaction should succeed");

    (void)label;
    return 0;
}

static int
apply_path_transaction(struct ged *gedp,
	ged_draw_transaction_kind kind,
	const char *path,
	void *view_ctx,
	int mode,
	const char *label)
{
    struct ged_draw_transaction txn = ged_draw_transaction_make(kind, path);
    txn.view = view_ctx;
    txn.mode = mode;
    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int ret = ged_draw_apply_transaction(gedp, &txn, &result);
    ged_draw_transaction_result_free(&result);
    if (ret <= 0)
	FAIL("public path transaction should succeed");

    (void)label;
    return 0;
}

static int
try_path_transaction(struct ged *gedp,
	ged_draw_transaction_kind kind,
	const char *path,
	void *view_ctx,
	int mode)
{
    struct ged_draw_transaction txn = ged_draw_transaction_make(kind, path);
    txn.view = view_ctx;
    txn.mode = mode;
    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int ret = ged_draw_apply_transaction(gedp, &txn, &result);
    ged_draw_transaction_result_free(&result);
    if (ret < 0)
	FAIL("public path transaction should not fail");

    return ret;
}

static int
exercise_mode_specific_source_lifecycle(struct ged *gedp,
	SoBRLSceneController *controller,
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

    (void)try_path_transaction(gedp, GED_DRAW_TXN_ERASE, path,
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

    if (apply_mode_value_transaction(gedp, GED_DRAW_TXN_VISIBILITY,
	    path, mode, 0.0, label))
	return 1;
    if (verify_mode_source(controller, path, representation_mode,
	    expect_vlist, expect_mesh, 0, 0, label))
	return 1;

    if (apply_mode_value_transaction(gedp, GED_DRAW_TXN_VISIBILITY,
	    path, mode, 1.0, label))
	return 1;
    if (apply_mode_path_transaction(gedp, GED_DRAW_TXN_STALE_SOURCE,
	    path, mode, label))
	return 1;
    if (verify_mode_source(controller, path, representation_mode,
	    expect_vlist, expect_mesh, 1, 1, label))
	return 1;

    if (apply_mode_path_transaction(gedp, GED_DRAW_TXN_REDRAW,
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

    if (apply_path_transaction(gedp, GED_DRAW_TXN_ERASE, path,
	    ged_draw_active_view_ctx(gedp), mode, "mode-specific erase"))
	return 1;
    if (source_representation_count(controller, path, representation_mode))
	FAIL("mode-specific erase should not leave target representation sources");
    if (!source_for_path(controller, path))
	FAIL("mode-specific erase should preserve the shared wire source");

    return 0;
}

static int
exercise_mesh_source_local_publication(struct ged *gedp,
	SoBRLSceneController *controller,
	const char *path)
{
    if (!gedp || !controller || !path)
	FAIL("mesh source-local test needs GED and Obol scene state");

    const int mode = GED_DRAW_MODE_SHADED_BOTS;
    const int representation = SoBRLDatabaseSource::REPRESENTATION_SHADED_BOTS;
    (void)try_path_transaction(gedp, GED_DRAW_TXN_ERASE, path,
	    ged_draw_active_view_ctx(gedp), mode);

    if (!ged_draw_obol_database_source_publication_begin(gedp,
	    ged_draw_active_view_ctx(gedp), mode))
	FAIL("mesh source-local publication scope should begin");
    if (!ged_draw_obol_database_source_ensure_for_path(gedp, path,
	    gedp->dbip, mode, 0)) {
	ged_draw_obol_database_source_publication_end(gedp);
	FAIL("mesh source-local test ensure should succeed");
    }

    SoBRLDatabaseSource *source =
	source_for_representation(controller, path, representation);
    BRLObolDatabaseSourceSummary summary;
    if (!source || !source->getSummary(summary) || !summary.valid ||
	    summary.instanceKey.getLength() == 0) {
	ged_draw_obol_database_source_publication_end(gedp);
	FAIL("mesh source-local test should find mode source");
    }

    SbMatrix placement = SbMatrix::identity();
    placement.setTranslate(SbVec3f(24.0f, 1.0f, 2.0f));
    if (controller->setDatabaseSourceInstancePlacementState(
		summary.instanceKey.getString(), TRUE, placement,
		FALSE, SbVec3f(0.0f, 0.0f, 0.0f), FALSE, 0.0f) < 0) {
	ged_draw_obol_database_source_publication_end(gedp);
	FAIL("mesh source-local test should set source placement");
    }

    point_t source_local_points[4] = {
	{1.0, 1.0, 1.0},
	{2.0, 1.0, 1.0},
	{1.0, 2.0, 1.0},
	{2.0, 2.0, 1.0}
    };
    vect_t normals[4] = {
	{0.0, 0.0, 1.0},
	{0.0, 0.0, 1.0},
	{0.0, 0.0, 1.0},
	{0.0, 0.0, 1.0}
    };
    int indices[5] = {0, 1, 3, 2, -1};

    if (!ged_draw_obol_database_source_publish_indexed_face_set_for_path(
	    gedp, path, (const point_t *)source_local_points, 4,
	    (const vect_t *)normals, 4, indices, 5)) {
	ged_draw_obol_database_source_publication_end(gedp);
	FAIL("mesh source-local publication should succeed");
    }
    ged_draw_obol_database_source_publication_end(gedp);

    source = source_for_representation(controller, path, representation);
    SoBRLMeshShape *mesh = source ? source->getRealizedMesh() : NULL;
    if (!mesh || mesh->point.getNum() != 4 ||
	    fabs(mesh->point[0][0] - 1.0f) > 0.001f ||
	    fabs(mesh->point[0][1] - 1.0f) > 0.001f ||
	    fabs(mesh->point[0][2] - 1.0f) > 0.001f ||
	    fabs(mesh->point[3][0] - 2.0f) > 0.001f ||
	    fabs(mesh->point[3][1] - 2.0f) > 0.001f ||
	    fabs(mesh->point[3][2] - 1.0f) > 0.001f)
	FAIL("mesh publication should store source-local geometry under source placement");

    if (!source->getSummary(summary) || !summary.valid ||
	    !summary.drawMatrixValid ||
	    !summary.drawMatrix.equals(placement, 0.0001f))
	FAIL("mesh source-local publication should preserve source placement");

    if (apply_path_transaction(gedp, GED_DRAW_TXN_ERASE, path,
	    ged_draw_active_view_ctx(gedp), mode,
	    "mesh source-local test cleanup"))
	return 1;

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
vlist_geometry_bounds(const SoBRLVListShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return 0;

    const SoBRLVListShape *geometry = shape->getGeometrySource();
    if (!geometry || geometry->point.getNum() <= 0)
	return 0;

    for (int i = 0; i < geometry->point.getNum(); i++)
	bounds.extendBy(geometry->point[i]);
    return bounds.isEmpty() ? 0 : 1;
}

static int
mesh_geometry_bounds(const SoBRLMeshShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return 0;

    const SoBRLMeshShape *geometry = shape->getGeometrySource();
    if (!geometry || geometry->point.getNum() <= 0)
	return 0;

    for (int i = 0; i < geometry->point.getNum(); i++)
	bounds.extendBy(geometry->point[i]);
    return bounds.isEmpty() ? 0 : 1;
}

static int
exercise_multi_instance_transform_reuse(struct ged *gedp,
	SoBRLSceneController *controller)
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

    SoBRLDatabaseSource *source_a = source_for_representation(controller,
	    path_a, SoBRLDatabaseSource::REPRESENTATION_WIRE);
    SoBRLDatabaseSource *source_b = source_for_representation(controller,
	    path_b, SoBRLDatabaseSource::REPRESENTATION_WIRE);
    if (!source_a || !source_b)
	FAIL("GED multi-instance root draw should create both transformed leaf sources");

    BRLObolDatabaseSourceSummary summary_a;
    BRLObolDatabaseSourceSummary summary_b;
    if (!source_a->getSummary(summary_a) || !summary_a.valid ||
	    !source_b->getSummary(summary_b) || !summary_b.valid)
	FAIL("GED multi-instance transformed sources should report summaries");

    SoBRLVListShape *shape_a = source_a->getRealizedShape();
    SoBRLVListShape *shape_b = source_b->getRealizedShape();
    const SoBRLVListShape *geometry_a =
	shape_a ? shape_a->getGeometrySource() : NULL;
    const SoBRLVListShape *geometry_b =
	shape_b ? shape_b->getGeometrySource() : NULL;
    if (!shape_a || !shape_b || !geometry_a || !geometry_b ||
	    geometry_a == shape_a || geometry_b == shape_b ||
	    geometry_a != geometry_b)
	FAIL("GED multi-instance transformed sources should reuse one local geometry node");

    SbBox3f local_bounds;
    if (!vlist_geometry_bounds(shape_a, local_bounds) ||
	    !box3f_near(local_bounds, -1.0f, -1.0f, -1.0f,
		1.0f, 1.0f, 1.0f))
	FAIL("GED multi-instance shared geometry should remain source-local");

    if (!summary_a.sourceBoundsValid ||
	    !box3f_near(summary_a.sourceBounds, -13.0f, 2.0f, -1.0f,
		-11.0f, 4.0f, 1.0f) ||
	    !summary_b.sourceBoundsValid ||
	    !box3f_near(summary_b.sourceBounds, 17.0f, 2.0f, -1.0f,
		19.0f, 4.0f, 1.0f))
	FAIL("GED multi-instance source bounds should include explicit tree-walk transforms");

    SbBox3f scene_bounds_a;
    SbBox3f scene_bounds_b;
    SbBox3f scene_bounds_root;
    if (!controller->getSceneSubtreeBounds(path_a, TRUE, scene_bounds_a) ||
	    !box3f_near(scene_bounds_a, -13.0f, 2.0f, -1.0f,
		-11.0f, 4.0f, 1.0f) ||
	    !controller->getSceneSubtreeBounds(path_b, TRUE, scene_bounds_b) ||
	    !box3f_near(scene_bounds_b, 17.0f, 2.0f, -1.0f,
		19.0f, 4.0f, 1.0f) ||
	    !controller->getSceneSubtreeBounds("reuse_root.c", TRUE,
		scene_bounds_root) ||
	    !box3f_near(scene_bounds_root, -13.0f, 2.0f, -1.0f,
		19.0f, 4.0f, 1.0f))
	FAIL("GED multi-instance scene bounds should apply explicit transforms");

    if (apply_mode_value_transaction(gedp, GED_DRAW_TXN_VISIBILITY,
	    path_a, GED_DRAW_MODE_WIRE, 0.0, "multi-instance visibility"))
	return 1;
    if (!source_a->getSummary(summary_a) || !summary_a.valid ||
	    !source_b->getSummary(summary_b) || !summary_b.valid ||
	    summary_a.visible ||
	    summary_b.visible ||
	    shape_a->visible.getValue() ||
	    shape_b->visible.getValue())
	FAIL("GED multi-instance visibility should update repeated logical instances consistently");
    if (apply_mode_value_transaction(gedp, GED_DRAW_TXN_VISIBILITY,
	    path_a, GED_DRAW_MODE_WIRE, 1.0, "multi-instance visibility restore"))
	return 1;
    if (!source_a->getSummary(summary_a) || !summary_a.valid ||
	    !source_b->getSummary(summary_b) || !summary_b.valid ||
	    !summary_a.visible ||
	    !summary_b.visible ||
	    !shape_a->visible.getValue() ||
	    !shape_b->visible.getValue())
	FAIL("GED multi-instance visibility restore should update repeated logical instances consistently");

    if (apply_mode_value_transaction(gedp, GED_DRAW_TXN_HIGHLIGHT,
	    "reuse_root.c", GED_DRAW_MODE_WIRE, 1.0,
	    "multi-instance highlight"))
	return 1;
    if (!source_a->getSummary(summary_a) || !summary_a.valid ||
	    !source_b->getSummary(summary_b) || !summary_b.valid ||
	    !summary_a.highlighted ||
	    !summary_b.highlighted ||
	    !shape_a->highlighted.getValue() ||
	    !shape_b->highlighted.getValue())
	FAIL("GED multi-instance highlight should update repeated logical instances consistently");
    if (apply_mode_value_transaction(gedp, GED_DRAW_TXN_HIGHLIGHT,
	    "reuse_root.c", GED_DRAW_MODE_WIRE, 0.0,
	    "multi-instance highlight restore"))
	return 1;
    if (!source_a->getSummary(summary_a) || !summary_a.valid ||
	    !source_b->getSummary(summary_b) || !summary_b.valid ||
	    summary_a.highlighted ||
	    summary_b.highlighted ||
	    shape_a->highlighted.getValue() ||
	    shape_b->highlighted.getValue())
	FAIL("GED multi-instance highlight restore should update repeated logical instances consistently");

    const SoBRLVListShape *seed_geometry = geometry_b;
    if (summary_a.instanceKey.getLength() <= 0 ||
	    controller->markDatabaseSourceInstanceStale(
		summary_a.instanceKey.getString(),
		SoBRLDatabaseSource::STALE_DRAW) <= 0)
	FAIL("GED multi-instance direct source stale should target one transformed source");
    if (!source_a->getSummary(summary_a) || !summary_a.valid ||
	    !source_b->getSummary(summary_b) || !summary_b.valid ||
	    !summary_a.stale ||
	    !(summary_a.staleReason & SoBRLDatabaseSource::STALE_DRAW) ||
	    summary_b.stale)
	FAIL("GED multi-instance direct source stale should preserve sibling current state");
    if (!controller->realizePending())
	FAIL("GED multi-instance direct source redraw should realize pending source");
    source_a = source_for_representation(controller, path_a,
	    SoBRLDatabaseSource::REPRESENTATION_WIRE);
    source_b = source_for_representation(controller, path_b,
	    SoBRLDatabaseSource::REPRESENTATION_WIRE);
    shape_a = source_a ? source_a->getRealizedShape() : NULL;
    shape_b = source_b ? source_b->getRealizedShape() : NULL;
    geometry_a = shape_a ? shape_a->getGeometrySource() : NULL;
    geometry_b = shape_b ? shape_b->getGeometrySource() : NULL;
    if (!source_a || !source_b ||
	    !source_a->getSummary(summary_a) || !summary_a.valid ||
	    !source_b->getSummary(summary_b) || !summary_b.valid ||
	    summary_a.stale ||
	    summary_b.stale ||
	    !shape_a || !shape_b || !geometry_a || !geometry_b ||
	    geometry_a != seed_geometry ||
	    geometry_b != seed_geometry)
	FAIL("GED multi-instance direct source redraw should reuse seeded sibling geometry");

    ged_draw_index_stats_reset(gedp);
    if (apply_mode_path_transaction(gedp, GED_DRAW_TXN_REDRAW,
	    path_a, GED_DRAW_MODE_WIRE, "multi-instance redraw"))
	return 1;
    struct ged_draw_index_stats redraw_stats;
    memset(&redraw_stats, 0, sizeof(redraw_stats));
    ged_draw_index_stats_get(gedp, &redraw_stats);
    if (redraw_stats.slow_path_shape_scans ||
	    redraw_stats.slow_path_group_scans)
	FAIL("GED multi-instance logical redraw should avoid registry/index slow-path scans");
    source_a = source_for_representation(controller, path_a,
	    SoBRLDatabaseSource::REPRESENTATION_WIRE);
    source_b = source_for_representation(controller, path_b,
	    SoBRLDatabaseSource::REPRESENTATION_WIRE);
    shape_a = source_a ? source_a->getRealizedShape() : NULL;
    shape_b = source_b ? source_b->getRealizedShape() : NULL;
    geometry_a = shape_a ? shape_a->getGeometrySource() : NULL;
    geometry_b = shape_b ? shape_b->getGeometrySource() : NULL;
    if (!source_a || !source_b ||
	    !source_a->getSummary(summary_a) || !summary_a.valid ||
	    !source_b->getSummary(summary_b) || !summary_b.valid ||
	    summary_a.stale ||
	    summary_b.stale ||
	    !shape_a || !shape_b || !geometry_a || !geometry_b ||
	    geometry_a == shape_a || geometry_b == shape_b ||
	    geometry_a != geometry_b)
	FAIL("GED multi-instance logical redraw should preserve shared local geometry");

    const char *draw_reuse_root_shaded[3] = {
	"draw", "-m2", "reuse_root.c"
    };
    if (ged_exec_draw(gedp, 3, draw_reuse_root_shaded) != BRLCAD_OK)
	FAIL("GED multi-instance shaded root draw should succeed");
    SoBRLDatabaseSource *mesh_source_a =
	source_for_representation(controller, path_a,
		SoBRLDatabaseSource::REPRESENTATION_SHADED);
    SoBRLDatabaseSource *mesh_source_b =
	source_for_representation(controller, path_b,
		SoBRLDatabaseSource::REPRESENTATION_SHADED);
    BRLObolDatabaseSourceSummary mesh_summary_a;
    BRLObolDatabaseSourceSummary mesh_summary_b;
    if (!mesh_source_a || !mesh_source_b ||
	    !mesh_source_a->getSummary(mesh_summary_a) ||
	    !mesh_summary_a.valid ||
	    !mesh_source_b->getSummary(mesh_summary_b) ||
	    !mesh_summary_b.valid ||
	    mesh_summary_a.realizedMeshCount <= 0 ||
	    mesh_summary_b.realizedMeshCount <= 0)
	FAIL("GED multi-instance shaded draw should create transformed mesh sources");
    SoBRLMeshShape *mesh_a = mesh_source_a->getRealizedMesh();
    SoBRLMeshShape *mesh_b = mesh_source_b->getRealizedMesh();
    const SoBRLMeshShape *mesh_geometry_a =
	mesh_a ? mesh_a->getGeometrySource() : NULL;
    const SoBRLMeshShape *mesh_geometry_b =
	mesh_b ? mesh_b->getGeometrySource() : NULL;
    if (!mesh_a || !mesh_b || !mesh_geometry_a || !mesh_geometry_b ||
	    mesh_geometry_a == mesh_a || mesh_geometry_b == mesh_b ||
	    mesh_geometry_a != mesh_geometry_b)
	FAIL("GED multi-instance shaded sources should reuse one local mesh node");
    SbBox3f mesh_local_bounds;
    if (!mesh_geometry_bounds(mesh_a, mesh_local_bounds) ||
	    !box3f_near(mesh_local_bounds, -1.0f, -1.0f, -1.0f,
		1.0f, 1.0f, 1.0f))
	FAIL("GED multi-instance shaded shared mesh should remain source-local");
    if (!mesh_summary_a.sourceBoundsValid ||
	    !box3f_near(mesh_summary_a.sourceBounds, -13.0f, 2.0f, -1.0f,
		-11.0f, 4.0f, 1.0f) ||
	    !mesh_summary_b.sourceBoundsValid ||
	    !box3f_near(mesh_summary_b.sourceBounds, 17.0f, 2.0f, -1.0f,
		19.0f, 4.0f, 1.0f))
	FAIL("GED multi-instance shaded source bounds should apply transforms");

    const char *autoview_cmd[2] = {"autoview", NULL};
    if (ged_exec_autoview(gedp, 1, autoview_cmd) != BRLCAD_OK)
	FAIL("GED multi-instance autoview should succeed");
    mat_t view_center_mat;
    point_t view_center;
    bv_center_mat_get(view_center_mat, DRAW_TEST_BV_CONST(ged_draw_active_view_ctx(gedp)));
    MAT_DELTAS_GET_NEG(view_center, view_center_mat);
    if (bv_size_get(DRAW_TEST_BV_CONST(ged_draw_active_view_ctx(gedp))) < 31.9 ||
	    fabs(view_center[X] - 3.0) > 0.1)
	FAIL("GED multi-instance autoview should use transformed scene bounds");

    const char *erase_reuse_root[2] = {"erase", "reuse_root.c"};
    if (ged_exec_erase(gedp, 2, erase_reuse_root) != BRLCAD_OK)
	FAIL("GED multi-instance transform root erase should succeed");
    if (source_for_path(controller, path_a) ||
	    source_for_path(controller, path_b) ||
	    source_for_representation(controller, path_a,
		SoBRLDatabaseSource::REPRESENTATION_SHADED) ||
	    source_for_representation(controller, path_b,
		SoBRLDatabaseSource::REPRESENTATION_SHADED) ||
	    controller->getDatabaseSourceCount() != initial_source_count)
	FAIL("GED multi-instance transform cleanup should restore prior source state");

    return 0;
}

static int
exercise_duplicate_occurrence_pick_identity(struct ged *gedp,
	SoBRLSceneController *controller)
{
    if (!gedp || !controller)
	FAIL("duplicate occurrence pick test needs GED and Obol scene state");

    const int initial_source_count = controller->getDatabaseSourceCount();
    const char *path_a = "dup_twice.c/dup_leaf.s";
    const char *path_b = "dup_twice.c/dup_leaf.s@1";
    const char *draw_duplicate[2] = {"draw", "dup_twice.c"};
    if (ged_exec_draw(gedp, 2, draw_duplicate) != BRLCAD_OK)
	FAIL("GED duplicate occurrence draw should succeed");

    SoBRLDatabaseSource *source_a = source_for_representation(controller,
	    path_a, SoBRLDatabaseSource::REPRESENTATION_WIRE);
    SoBRLDatabaseSource *source_b = source_for_representation(controller,
	    path_b, SoBRLDatabaseSource::REPRESENTATION_WIRE);
    if (!source_a || !source_b)
	FAIL("GED duplicate occurrence draw should preserve both source instances");

    BRLObolDatabaseSourceSummary summary_a;
    BRLObolDatabaseSourceSummary summary_b;
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
    if (strcmp(shape_a->ownerSourceInstanceKey.getValue().getString(),
	    summary_a.instanceKey.getString()) != 0 ||
	    strcmp(shape_b->ownerSourceInstanceKey.getValue().getString(),
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
    if (!raw_detail ||
	    !raw_detail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("GED duplicate occurrence ray pick should return BRL-CAD pick detail");
    const SoBRLPickDetail *pick_detail =
	static_cast<const SoBRLPickDetail *>(raw_detail);
    if (strcmp(pick_detail->getPath().getString(),
	    summary_b.path.getString()) != 0 ||
	    strcmp(pick_detail->getSourceInstanceKey().getString(),
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
	if (strcmp(candidate->recordRole.getValue().getString(),
		"auxiliary") == 0)
	    return candidate;
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
	SoBRLSceneController *controller,
	const char *path)
{
    if (!gedp || !controller || !path)
	FAIL("evaluated-wire realize-context test needs GED and Obol scene state");

    const int mode = GED_DRAW_MODE_EVAL_WIRE;
    const int representation = SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE;
    (void)try_path_transaction(gedp, GED_DRAW_TXN_ERASE, path,
	    ged_draw_active_view_ctx(gedp), mode);

    char mode_arg[16] = {0};
    snprintf(mode_arg, sizeof(mode_arg), "-m%d", mode);
    const char *draw_mode_cmd[4] = {"draw", mode_arg, path, NULL};
    if (ged_exec_draw(gedp, 3, draw_mode_cmd) != BRLCAD_OK)
	FAIL("evaluated-wire realize-context setup draw should succeed");

    SoBRLDatabaseSource *source =
	source_for_representation(controller, path, representation);
    BRLObolDatabaseSourceSummary summary;
    if (!source || !source->getSummary(summary) || !summary.valid ||
	    summary.realizationStatus != SoBRLDatabaseSource::REALIZED ||
	    summary.realizedShapeCount <= 0 ||
	    summary.instanceKey.getLength() == 0)
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
	    summary.realizedShapeCount <= 0)
	FAIL("evaluated-wire shape-ref realize-context should refresh the mode source");

    (void)try_path_transaction(gedp, GED_DRAW_TXN_ERASE, path,
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
    BRLObolViewController *owned_controller = ged_draw_obol_controller(gedp);
    if (!owned_scene || ged_draw_obol_scene_controller(gedp) != owned_scene ||
	    !ged_draw_obol_scene_controller_owned(gedp) ||
	    !owned_controller ||
	    owned_controller->getSceneController() != owned_scene)
	FAIL("GED should create and report an owned Obol view scene");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should full-sync current GED draw state");
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

    void *feature_view_ctx = ged_view_active_ctx(gedp);
    if (!feature_view_ctx ||
	    !ged_draw_obol_view_context_feature_store_active(feature_view_ctx))
	FAIL("GED active view should expose the owned Obol feature store");
    point_t feature_points[2] = {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}};
    int feature_cmds[2] = {
	GED_DRAW_VIEW_LINE_MOVE,
	GED_DRAW_VIEW_LINE_DRAW
    };
    struct ged_draw_view_feature_style feature_style =
	GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    feature_style.visible = 1;
    feature_style.color_valid = 1;
    feature_style.color[0] = 20;
    feature_style.color[1] = 40;
    feature_style.color[2] = 60;
    feature_style.line_width = 3;
    if (!ged_draw_view_context_lines_replace(feature_view_ctx, "cap2::line",
	    0, feature_points, feature_cmds, 2, &feature_style) ||
	    !owned_controller->features().exists("cap2::line"))
	FAIL("GED feature line replacement should publish into the owned Obol feature store");
    struct ged_draw_view_feature_summary feature_summary =
	GED_DRAW_VIEW_FEATURE_SUMMARY_INIT;
    if (!ged_draw_view_context_feature_summary(feature_view_ctx,
	    "cap2::line", &feature_summary) ||
	    !feature_summary.exists ||
	    feature_summary.geometry_command_count != 2 ||
	    feature_summary.color[0] != 20 ||
	    feature_summary.color[1] != 40 ||
	    feature_summary.color[2] != 60)
	FAIL("GED feature summary should read the owned Obol feature store");
    point_t *copied_points = NULL;
    size_t copied_point_count = 0;
    if (!ged_draw_view_context_feature_points_copy(feature_view_ctx,
	    "cap2::line", &copied_points, &copied_point_count) ||
	    copied_point_count != 2 ||
	    copied_points[1][X] != 2.0)
	FAIL("GED feature point readback should copy owned Obol line geometry");
    bu_free(copied_points, "GED Obol feature test copied points");
    int copied_cmd = 0;
    if (!ged_draw_view_context_feature_line_command_at(feature_view_ctx,
	    "cap2::line", 1, &copied_cmd) ||
	    copied_cmd != GED_DRAW_VIEW_LINE_DRAW)
	FAIL("GED feature command readback should copy owned Obol line commands");
    point_t append_point = {3.0, 0.0, 0.0};
    if (!ged_draw_view_context_lines_append_point(feature_view_ctx,
	    "cap2::line", append_point) ||
	    !ged_draw_view_context_feature_points_copy(feature_view_ctx,
		"cap2::line", &copied_points, &copied_point_count) ||
	    copied_point_count != 3)
	FAIL("GED line append should mutate owned Obol line geometry");
    bu_free(copied_points, "GED Obol feature test appended points");
    if (!ged_draw_view_context_feature_visible_set(feature_view_ctx,
	    "cap2::line", 0) ||
	    ged_draw_view_context_feature_visible(feature_view_ctx,
		"cap2::line") != 0)
	FAIL("GED feature visibility should mutate owned Obol feature style");
    if (!ged_draw_view_context_line_color_set(feature_view_ctx,
	    "cap2::line", 90, 80, 70) ||
	    !ged_draw_view_context_line_width_set(feature_view_ctx,
		"cap2::line", 5))
	FAIL("GED line style setters should mutate owned Obol feature style");
    struct ged_draw_view_line_style line_style;
    memset(&line_style, 0, sizeof(line_style));
    if (!ged_draw_view_context_line_style_get(feature_view_ctx,
	    "cap2::line", &line_style) ||
	    line_style.color[0] != 90 ||
	    line_style.color[1] != 80 ||
	    line_style.color[2] != 70 ||
	    line_style.line_width != 5)
	FAIL("GED line style readback should read owned Obol feature style");

    point_t tcl_line_points[2] = {{0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}};
    struct ged_draw_view_line_style tcl_line_style;
    memset(&tcl_line_style, 0, sizeof(tcl_line_style));
    tcl_line_style.color[0] = 101;
    tcl_line_style.color[1] = 102;
    tcl_line_style.color[2] = 103;
    tcl_line_style.line_width = 4;
    if (!ged_draw_view_context_tcl_lines_replace(feature_view_ctx,
	    "cap2::tcl-line", tcl_line_points, 2, &tcl_line_style) ||
	    !feature_overlay_matches(owned_controller, "cap2::tcl-line",
		BRLObolOverlayClass::TclOverlay,
		BRLObolOverlayLifecycle::PerCommand,
		BRLObolOverlayOrder::PostTransparent))
	FAIL("GED Tcl line replacement should publish typed Obol overlay metadata");
    struct ged_draw_view_feature_summary tcl_line_summary =
	GED_DRAW_VIEW_FEATURE_SUMMARY_INIT;
    if (!ged_draw_view_context_feature_summary(feature_view_ctx,
	    "cap2::tcl-line", &tcl_line_summary) ||
	    !tcl_line_summary.exists ||
	    !tcl_line_summary.is_overlay ||
	    tcl_line_summary.geometry_command_count != 2)
	FAIL("GED Tcl line summary should read typed Obol overlay feature state");

    point_t annotation_point = {5.0, 5.0, 0.0};
    if (!ged_draw_view_context_lines_create_model_annotation(
	    feature_view_ctx, "cap2::annotation", 1, annotation_point) ||
	    !feature_overlay_matches(owned_controller, "cap2::annotation",
		BRLObolOverlayClass::UserAnnotation,
		BRLObolOverlayLifecycle::Persistent,
		BRLObolOverlayOrder::Model))
	FAIL("GED model annotation creation should publish typed Obol overlay metadata");

    point_t polygon_points[2] = {{1.0, 0.0, 0.0}, {1.0, 2.0, 0.0}};
    int polygon_cmds[2] = {
	GED_DRAW_VIEW_LINE_MOVE,
	GED_DRAW_VIEW_LINE_DRAW
    };
    if (!ged_draw_view_context_tcl_polygons_replace(feature_view_ctx,
	    "cap2::polygon-overlay", polygon_points, polygon_cmds, 2,
	    &feature_style) ||
	    !feature_overlay_matches(owned_controller,
		"cap2::polygon-overlay",
		BRLObolOverlayClass::PolygonEdit,
		BRLObolOverlayLifecycle::PerTool,
		BRLObolOverlayOrder::PostTransparent))
	FAIL("GED Tcl polygon replacement should publish typed Obol overlay metadata");

    struct ged_draw_view_label_data label = GED_DRAW_VIEW_LABEL_DATA_INIT;
    label.text = "cap2 label";
    VSET(label.point, 1.0, 2.0, 3.0);
    label.color_valid = 1;
    label.color[0] = 7;
    label.color[1] = 8;
    label.color[2] = 9;
    label.font_size = 18.0;
    if (!ged_draw_view_context_labels_replace(feature_view_ctx,
	    "cap2::label", 0, &label, 1) ||
	    !owned_controller->features().exists("cap2::label") ||
	    ged_draw_view_context_label_count(feature_view_ctx,
		"cap2::label") != 1)
	FAIL("GED label replacement should publish into the owned Obol feature store");
    BRLObolFeatureHandle label_handle =
	owned_controller->features().find("cap2::label");
    BRLObolFeatureRecord label_record;
    if (!label_handle.isValid() ||
	    !owned_controller->features().record(label_handle,
		label_record) ||
	    label_record.labels.size() != 1 ||
	    fabs(label_record.labels[0].fontSize - 18.0f) > 0.001f)
	FAIL("GED label replacement should preserve explicit Obol font size");
    struct bu_vls label_text = BU_VLS_INIT_ZERO;
    point_t label_point = VINIT_ZERO;
    unsigned char label_rgb[3] = {0, 0, 0};
    if (!ged_draw_view_context_label_copy(feature_view_ctx, "cap2::label",
	    0, &label_text, label_point, label_rgb) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&label_text), "cap2 label") ||
	    label_point[X] != 1.0 ||
	    label_rgb[0] != 7 ||
	    label_rgb[1] != 8 ||
	    label_rgb[2] != 9) {
	bu_vls_free(&label_text);
	FAIL("GED label copy should read owned Obol label data");
    }
    bu_vls_free(&label_text);
    point_t moved_label_point = {4.0, 5.0, 6.0};
    if (!ged_draw_view_context_label_point_set(feature_view_ctx,
	    "cap2::label", 0, moved_label_point) ||
	    !ged_draw_view_context_label_copy(feature_view_ctx, "cap2::label",
		0, NULL, label_point, NULL) ||
	    label_point[X] != 4.0 ||
	    label_point[Y] != 5.0 ||
	    label_point[Z] != 6.0)
	FAIL("GED label point mutation should update owned Obol label data");

    point_t created_label_point = {1.0, 1.0, 1.0};
    point_t created_label_target = {0.0, 0.0, 0.0};
    if (!ged_draw_view_context_label_create(feature_view_ctx,
	    "cap2::created-label", 0, "created", created_label_point,
	    created_label_target, 1) ||
	    !ged_draw_view_context_label_copy(feature_view_ctx,
		"cap2::created-label", 0, NULL, NULL, label_rgb) ||
	    label_rgb[0] != 255 ||
	    label_rgb[1] != 255 ||
	    label_rgb[2] != 0)
	FAIL("GED label create should preserve the legacy yellow feature color in Obol");

    point_t arrow_points[2] = {{0.0, 0.0, 0.0}, {0.0, 3.0, 0.0}};
    if (!ged_draw_view_context_tcl_arrows_replace(feature_view_ctx,
	    "cap2::arrow", arrow_points, 2, &feature_style) ||
	    !owned_controller->features().exists("cap2::arrow") ||
	    !feature_overlay_matches(owned_controller, "cap2::arrow",
		BRLObolOverlayClass::TclOverlay,
		BRLObolOverlayLifecycle::PerCommand,
		BRLObolOverlayOrder::PostTransparent) ||
	    !ged_draw_view_context_arrow_tip_set(feature_view_ctx,
		"cap2::arrow", 0.25, 0.5))
	FAIL("GED arrow replacement should publish into the owned Obol feature store");
    fastf_t tip_length = 0.0;
    fastf_t tip_width = 0.0;
    if (!ged_draw_view_context_arrow_tip_get(feature_view_ctx,
	    "cap2::arrow", &tip_length, &tip_width) ||
	    fabs(tip_length - 0.25) > 0.001 ||
	    fabs(tip_width - 0.5) > 0.001)
	FAIL("GED arrow tip readback should read owned Obol feature style");

    struct ged_draw_view_axes_state axes_state;
    memset(&axes_state, 0, sizeof(axes_state));
    VSET(axes_state.position, 1.0, 1.0, 1.0);
    axes_state.size = 4.0;
    axes_state.line_width = 2;
    axes_state.color[0] = 11;
    axes_state.color[1] = 22;
    axes_state.color[2] = 33;
    if (!ged_draw_view_context_axes_create(feature_view_ctx, "cap2::axes",
	    0, &axes_state) ||
	    !owned_controller->features().exists("cap2::axes"))
	FAIL("GED axes creation should publish into the owned Obol feature store");
    struct ged_draw_view_axes_state axes_readback;
    memset(&axes_readback, 0, sizeof(axes_readback));
    if (!ged_draw_view_context_axes_state_get(feature_view_ctx,
	    "cap2::axes", &axes_readback) ||
	    axes_readback.position[X] != 1.0 ||
	    axes_readback.size != 4.0 ||
	    axes_readback.line_width != 2 ||
	    axes_readback.color[0] != 11 ||
	    axes_readback.color[1] != 22 ||
	    axes_readback.color[2] != 33)
	FAIL("GED axes readback should read owned Obol axes state");
    point_t *axis_centers = NULL;
    size_t axis_center_count = 0;
    if (!ged_draw_view_context_feature_axes_centers_copy(feature_view_ctx,
	    "cap2::axes", &axis_centers, &axis_center_count) ||
	    axis_center_count != 1 ||
	    axis_centers[0][X] != 1.0)
	FAIL("GED axes center readback should copy owned Obol axes centers");
    bu_free(axis_centers, "GED Obol feature test axes centers");

    point_t tcl_axes_centers[1] = {{2.0, 2.0, 0.0}};
    if (!ged_draw_view_context_tcl_axes_replace(feature_view_ctx,
	    "cap2::tcl-axes", tcl_axes_centers, 1, 2.5, &feature_style) ||
	    !feature_overlay_matches(owned_controller, "cap2::tcl-axes",
		BRLObolOverlayClass::TclOverlay,
		BRLObolOverlayLifecycle::PerCommand,
		BRLObolOverlayOrder::PostTransparent))
	FAIL("GED Tcl axes replacement should publish typed Obol overlay metadata");

    point_t face_points[4] = {
	{0.0, 0.0, 0.0},
	{1.0, 0.0, 0.0},
	{1.0, 1.0, 0.0},
	{0.0, 1.0, 0.0}
    };
    int face_indices[6] = {0, 1, 2, 0, 2, 3};
    if (!ged_draw_view_context_indexed_face_set_replace(feature_view_ctx,
	    "cap2::mesh", 0, face_points, 4, NULL, 0, face_indices, 6,
	    &feature_style) ||
	    !owned_controller->features().exists("cap2::mesh"))
	FAIL("GED indexed-face replacement should publish into the owned Obol feature store");

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
    if (!ged_draw_view_context_diagnostic_line_layer_builder_replace(
		feature_view_ctx, "cap2::diagnostic", diagnostic_builder)) {
	bg_line_layer_builder_free(diagnostic_builder);
	FAIL("GED diagnostic line-layer replacement should publish into the owned Obol feature store");
    }
    bg_line_layer_builder_free(diagnostic_builder);
    BRLObolFeatureRecord diagnostic_record;
    BRLObolFeatureHandle diagnostic_handle =
	owned_controller->features().find("cap2::diagnostic");
    if (!diagnostic_handle.isValid() ||
	    !owned_controller->features().record(diagnostic_handle,
		diagnostic_record) ||
	    diagnostic_record.kind != BRLObolFeatureKind::LineLayer ||
	    !feature_overlay_matches(owned_controller, "cap2::diagnostic",
		BRLObolOverlayClass::Diagnostic,
		BRLObolOverlayLifecycle::PerCommand,
		BRLObolOverlayOrder::PostTransparent))
	FAIL("GED diagnostic line-layer replacement should stamp typed Obol diagnostic metadata");

    struct command_result_callback_state command_callback_state = {};
    struct ged_draw_command_scene_desc command_scene_desc =
	GED_DRAW_COMMAND_SCENE_DESC_INIT;
    command_scene_desc.owner_id = "rtcheck";
    command_scene_desc.owner_role = "command-result";
    command_scene_desc.result_cb = command_result_cb;
    command_scene_desc.result_cb_data = &command_callback_state;
    struct ged_draw_command_scene *command_scene =
	ged_draw_command_scene_begin(feature_view_ctx, &command_scene_desc);
    if (!command_scene)
	FAIL("GED command-scene begin should create an Obol-backed publication context");
    struct ged_draw_view_feature_style command_style =
	GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    command_style.visible = 1;
    command_style.selectable = 0;
    struct ged_draw_view_line_layer_data command_layer =
	GED_DRAW_VIEW_LINE_LAYER_DATA_INIT;
    command_layer.name = "rtcheck::overlaps/yellow";
    command_layer.points = feature_points;
    command_layer.commands = feature_cmds;
    command_layer.point_count = 2;
    struct ged_draw_command_scene_metadata command_metadata[2] = {
	{"result.kind", "overlap"},
	{"result.count", "1"}
    };
    struct ged_draw_command_scene_metadata command_primitive_metadata[2] = {
	{"overlap.objects", "box.s cone.s"},
	{"overlap.depth", "0.25"}
    };
    if (!ged_draw_command_scene_line_layers_replace(command_scene,
	    "rtcheck::overlaps", &command_layer, 1, &command_style) ||
	    !ged_draw_command_scene_feature_metadata_replace(command_scene,
		"rtcheck::overlaps", command_metadata, 2) ||
	    !ged_draw_command_scene_feature_primitive_metadata_replace(
		command_scene, "rtcheck::overlaps", 0,
		command_primitive_metadata, 2) ||
	    !ged_draw_command_scene_commit(command_scene))
	FAIL("GED command-scene line-layer replacement should commit");
    BRLObolFeatureHandle command_handle =
	owned_controller->features().find("rtcheck::overlaps",
		BRLOBOL_FEATURE_SCOPE_SHARED);
    BRLObolFeatureRecord command_record;
    if (!command_handle.isValid() ||
	    !owned_controller->features().record(command_handle,
		command_record) ||
	    command_record.scope != BRLObolFeatureScope::Shared ||
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
		BRLObolOverlayClass::CommandResult,
		BRLObolOverlayLifecycle::PerCommand,
		BRLObolOverlayOrder::PostTransparent))
	FAIL("GED command-scene result should be shared, owned, selectable-aware command content");
    if (command_callback_state.accepted_count < 2 ||
	    command_callback_state.updated_count < 2 ||
	    !command_callback_state.saw_line_layers_update ||
	    !command_callback_state.saw_metadata_update ||
	    !command_callback_state.saw_primitive_metadata_update ||
	    command_callback_state.line_layers_feature_id != command_handle.id ||
	    command_callback_state.metadata_feature_id != command_handle.id ||
	    command_callback_state.primitive_metadata_feature_id !=
		command_handle.id)
	FAIL("GED command-scene result callback should report line-layer and metadata updates with feature handles");
    struct bu_vls primitive_key = BU_VLS_INIT_ZERO;
    struct bu_vls primitive_value = BU_VLS_INIT_ZERO;
    if (ged_draw_view_context_feature_primitive_metadata_count(
		feature_view_ctx, "rtcheck::overlaps", 0) != 2 ||
	    !ged_draw_view_context_feature_primitive_metadata_copy(
		feature_view_ctx, "rtcheck::overlaps", 0, 0,
		&primitive_key, &primitive_value) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&primitive_key), "overlap.objects") ||
	    !BU_STR_EQUAL(bu_vls_cstr(&primitive_value), "box.s cone.s")) {
	bu_vls_free(&primitive_key);
	bu_vls_free(&primitive_value);
	FAIL("GED command-scene primitive metadata should read back through neutral view APIs");
    }
    bu_vls_free(&primitive_key);
    bu_vls_free(&primitive_value);

    struct bu_vls resolved_feature = BU_VLS_INIT_ZERO;
    int resolved_primitive = -1;
    int primitive_index = -1;
    if (!ged_draw_view_context_feature_pick_primitive_resolve(
		feature_view_ctx, "rtcheck::overlaps/yellow", 0, 1, 1,
		&resolved_feature, &resolved_primitive) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&resolved_feature),
		"rtcheck::overlaps") ||
	    resolved_primitive != 0 ||
	    ged_draw_view_context_feature_selected_primitive_count(
		feature_view_ctx, "rtcheck::overlaps") != 1 ||
	    ged_draw_view_context_feature_highlighted_primitive_count(
		feature_view_ctx, "rtcheck::overlaps") != 1 ||
	    !ged_draw_view_context_feature_selected_primitive_at(
		feature_view_ctx, "rtcheck::overlaps", 0, &primitive_index) ||
	    primitive_index != 0 ||
	    !ged_draw_view_context_feature_highlighted_primitive_at(
		feature_view_ctx, "rtcheck::overlaps", 0, &primitive_index) ||
	    primitive_index != 0) {
	bu_vls_free(&resolved_feature);
	FAIL("GED command-scene child primitive picks should resolve and set parent primitive state");
    }
    bu_vls_free(&resolved_feature);
    if (ged_draw_view_context_feature_pick_primitive_resolve(
	    feature_view_ctx, "rtcheck::overlaps/yellow", 1, 0, 0, NULL,
	    &resolved_primitive))
	FAIL("GED command-scene child primitive resolver should reject out-of-range picks");
    struct ged_draw_view_feature_summary command_summary =
	GED_DRAW_VIEW_FEATURE_SUMMARY_INIT;
    if (!ged_draw_view_context_feature_summary(feature_view_ctx,
	    "rtcheck::overlaps", &command_summary) ||
	    command_summary.primitive_metadata_count != 1 ||
	    command_summary.selected_primitive_count != 1 ||
	    command_summary.highlighted_primitive_count != 1)
	FAIL("GED command-scene result summary should report primitive state");
    if (!owned_controller->features().record(command_handle,
		command_record) ||
	    command_record.selectedPrimitives.size() != 1 ||
	    command_record.highlightedPrimitives.size() != 1)
	FAIL("GED command-scene result record should preserve primitive state");
    SoBRLVListShape *command_vlist = first_feature_vlist(
	    owned_controller->features().node(command_handle));
    if (!command_vlist ||
	    command_vlist->selectedPrimitive.getNum() != 1 ||
	    command_vlist->selectedPrimitive[0] != 0 ||
	    command_vlist->highlightedPrimitive.getNum() != 1 ||
	    command_vlist->highlightedPrimitive[0] != 0)
	FAIL("GED command-scene result primitive state should reach realized Coin VLIST");

    command_scene = ged_draw_command_scene_begin(feature_view_ctx,
	    &command_scene_desc);
    if (!command_scene ||
	    ged_draw_command_scene_features_remove_prefix(command_scene,
		"rtcheck::") != 1 ||
	    !ged_draw_command_scene_commit(command_scene) ||
	    owned_controller->features().exists("rtcheck::overlaps"))
	FAIL("GED command-scene remove-prefix should remove owned shared command results");
    if (command_callback_state.removed_count < 1 ||
	    !command_callback_state.saw_remove_prefix)
	FAIL("GED command-scene result callback should report owner-scoped feature removal");

    command_scene_desc.generation = 10;
    struct ged_draw_command_scene *stale_scene =
	ged_draw_command_scene_begin(feature_view_ctx, &command_scene_desc);
    if (!stale_scene)
	FAIL("GED command-scene stale generation test should create old scene");

    point_t latest_points[2] = {
	{0.0, 0.0, 0.0},
	{0.0, 2.0, 0.0}
    };
    struct ged_draw_view_line_layer_data latest_layer =
	GED_DRAW_VIEW_LINE_LAYER_DATA_INIT;
    latest_layer.name = "rtcheck::generation/latest";
    latest_layer.points = latest_points;
    latest_layer.commands = feature_cmds;
    latest_layer.point_count = 2;
    command_scene_desc.generation = 11;
    command_scene = ged_draw_command_scene_begin(feature_view_ctx,
	    &command_scene_desc);
    if (!command_scene ||
	    !ged_draw_command_scene_line_layers_replace(command_scene,
		"rtcheck::generation", &latest_layer, 1, &command_style) ||
	    !ged_draw_command_scene_commit(command_scene))
	FAIL("GED command-scene latest generation should publish");

    point_t stale_points[2] = {
	{0.0, 0.0, 0.0},
	{0.0, 3.0, 0.0}
    };
    struct ged_draw_view_line_layer_data stale_layer =
	GED_DRAW_VIEW_LINE_LAYER_DATA_INIT;
    stale_layer.name = "rtcheck::generation/stale";
    stale_layer.points = stale_points;
    stale_layer.commands = feature_cmds;
    stale_layer.point_count = 2;
    if (ged_draw_command_scene_line_layers_replace(stale_scene,
	    "rtcheck::generation", &stale_layer, 1, &command_style) ||
	    ged_draw_command_scene_commit(stale_scene))
	FAIL("GED command-scene stale generation should be rejected");
    if (command_callback_state.failed_count < 2 ||
	    !command_callback_state.saw_stale_failure ||
	    !command_callback_state.saw_commit_failure)
	FAIL("GED command-scene result callback should report stale generation rejection");

    command_handle = owned_controller->features().find("rtcheck::generation",
	    BRLOBOL_FEATURE_SCOPE_SHARED);
    if (!command_handle.isValid() ||
	    !owned_controller->features().record(command_handle,
		command_record) ||
	    command_record.owner.generation != 11 ||
	    command_record.points.size() != 2 ||
	    command_record.points[1][1] != 2.0f)
	FAIL("GED command-scene stale generation should not replace latest result");
    command_scene = ged_draw_command_scene_begin(feature_view_ctx,
	    &command_scene_desc);
    if (!command_scene ||
	    ged_draw_command_scene_features_remove_prefix(command_scene,
		"rtcheck::generation") != 1 ||
	    !ged_draw_command_scene_commit(command_scene) ||
	    owned_controller->features().exists("rtcheck::generation"))
	FAIL("GED command-scene generation cleanup should remove latest result");

    struct custom_node_provider_state custom_provider_state = {};
    struct ged_draw_view_feature_style custom_style =
	GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    custom_style.visible = 1;
    custom_style.selectable = 1;
    struct ged_draw_command_scene_metadata custom_metadata[1] = {
	{"result.kind", "custom-node"}
    };
    command_scene_desc.owner_id = "custom";
    command_scene_desc.owner_role = "command-result";
    command_scene_desc.generation = 33;
    command_scene = ged_draw_command_scene_begin(feature_view_ctx,
	    &command_scene_desc);
    if (!command_scene ||
	    !ged_draw_command_scene_custom_node_replace(command_scene,
		"custom::node", custom_node_provider_cb,
		&custom_provider_state, &custom_style) ||
	    !ged_draw_command_scene_feature_metadata_replace(command_scene,
		"custom::node", custom_metadata, 1) ||
	    !ged_draw_command_scene_commit(command_scene))
	FAIL("GED command-scene custom Coin node provider should commit");
    if (custom_provider_state.call_count != 1 ||
	    !custom_provider_state.saw_request ||
	    custom_provider_state.generation != 33 ||
	    custom_provider_state.local)
	FAIL("GED command-scene custom Coin provider should receive owner/scope request metadata");
    BRLObolFeatureHandle custom_handle =
	owned_controller->features().find("custom::node",
		BRLOBOL_FEATURE_SCOPE_SHARED);
    BRLObolFeatureRecord custom_record;
    if (!custom_handle.isValid() ||
	    !owned_controller->features().record(custom_handle,
		custom_record) ||
	    custom_record.kind != BRLObolFeatureKind::CustomNode ||
	    custom_record.scope != BRLObolFeatureScope::Shared ||
	    custom_record.owner.generation != 33 ||
	    custom_record.metadata.size() != 1 ||
	    !BU_STR_EQUAL(custom_record.metadata[0].value.getString(),
		"custom-node") ||
	    owned_controller->features().node(custom_handle) !=
		custom_provider_state.node ||
	    !feature_overlay_matches(owned_controller, "custom::node",
		BRLObolOverlayClass::CommandResult,
		BRLObolOverlayLifecycle::PerCommand,
		BRLObolOverlayOrder::PostTransparent))
	FAIL("GED command-scene custom Coin node should be an owned shared command-result feature");
    int custom_primitive = 7;
    if (!ged_draw_view_context_feature_selected_primitives_replace(
		feature_view_ctx, "custom::node", &custom_primitive, 1) ||
	    !ged_draw_view_context_feature_highlighted_primitives_replace(
		feature_view_ctx, "custom::node", &custom_primitive, 1) ||
	    !owned_controller->features().record(custom_handle,
		custom_record) ||
	    custom_record.selectedPrimitives.size() != 1 ||
	    custom_record.selectedPrimitives[0] != 7 ||
	    custom_record.highlightedPrimitives.size() != 1 ||
	    custom_record.highlightedPrimitives[0] != 7 ||
	    owned_controller->features().node(custom_handle) !=
		custom_provider_state.node)
	FAIL("GED command-scene custom Coin node should preserve primitive state without replacing the provider node");
    if (!command_callback_state.saw_custom_update ||
	    command_callback_state.custom_feature_id != custom_handle.id)
	FAIL("GED command-scene custom Coin node should report result callbacks");
    command_scene = ged_draw_command_scene_begin(feature_view_ctx,
	    &command_scene_desc);
    if (!command_scene ||
	    ged_draw_command_scene_features_remove_prefix(command_scene,
		"custom::") != 1 ||
	    !ged_draw_command_scene_commit(command_scene) ||
	    owned_controller->features().exists("custom::node"))
	FAIL("GED command-scene custom Coin node cleanup should use owner-scoped removal");

    command_scene_desc.generation = 0;
    command_scene_desc.result_cb = NULL;
    command_scene_desc.result_cb_data = NULL;

    FILE *nirt_plot = tmpfile();
    if (!nirt_plot)
	FAIL("NIRT/qray command-scene uplot test should create a temporary plot stream");
    int old_plot_mode = pl_getOutputMode();
    pl_setOutputMode(PL_OUTPUT_MODE_BINARY);
    point_t nirt_a = VINIT_ZERO;
    point_t nirt_b = {1.0, 0.0, 0.0};
    pl_color(nirt_plot, 0, 255, 0);
    pdv_3line(nirt_plot, nirt_a, nirt_b);
    pl_setOutputMode(old_plot_mode);
    rewind(nirt_plot);
    if (_ged_draw_uplot_to_command_scene_feature(gedp, nirt_plot,
	    "query_ray", 1.0, PL_OUTPUT_MODE_BINARY, "nirt",
	    "command-result", "query_ray", "query-ray", 0) != BRLCAD_OK) {
	fclose(nirt_plot);
	FAIL("NIRT/qray uplot import should publish through command-scene ownership");
    }
    fclose(nirt_plot);
    BRLObolFeatureHandle nirt_handle =
	owned_controller->features().find("query_ray",
		BRLOBOL_FEATURE_SCOPE_SHARED);
    BRLObolFeatureRecord nirt_record;
    if (!nirt_handle.isValid() ||
	    !owned_controller->features().record(nirt_handle, nirt_record) ||
	    nirt_record.kind != BRLObolFeatureKind::LineLayer ||
	    nirt_record.scope != BRLObolFeatureScope::Shared ||
	    !BU_STR_EQUAL(nirt_record.owner.ownerId.getString(), "nirt") ||
	    !BU_STR_EQUAL(nirt_record.owner.ownerRole.getString(),
		"command-result") ||
	    nirt_record.layers.size() != 1 ||
	    nirt_record.points.size() != 2 ||
	    nirt_record.metadata.size() < 6 ||
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
	    !feature_overlay_matches(owned_controller, "query_ray",
		BRLObolOverlayClass::CommandResult,
		BRLObolOverlayLifecycle::PerCommand,
		BRLObolOverlayOrder::PostTransparent))
	FAIL("NIRT/qray command-scene uplot result should be shared owned command content");
    command_scene_desc.owner_id = "nirt";
    command_scene = ged_draw_command_scene_begin(feature_view_ctx,
	    &command_scene_desc);
    if (!command_scene ||
	    ged_draw_command_scene_features_remove_prefix(command_scene,
		"query_ray") != 1 ||
	    !ged_draw_command_scene_commit(command_scene) ||
	    owned_controller->features().exists("query_ray"))
	FAIL("NIRT/qray command-scene cleanup should remove owned shared command results");

    struct bg_line_layer_builder *builder_publish =
	bg_line_layer_builder_create();
    if (!builder_publish)
	FAIL("command-scene builder publish test should allocate a builder");
    point_t builder_a = {0.0, 0.0, 0.0};
    point_t builder_b = {0.0, 1.0, 0.0};
    if (!bg_line_layer_builder_add(builder_publish, 255, 255, 0,
		builder_a, BG_GEOMETRY_LINE_MOVE) ||
	    !bg_line_layer_builder_add(builder_publish, 255, 255, 0,
		builder_b, BG_GEOMETRY_LINE_DRAW)) {
	bg_line_layer_builder_free(builder_publish);
	FAIL("command-scene builder publish test should accept line geometry");
    }
    if (_ged_line_layer_builder_publish_command_scene_feature(gedp,
	    "nmg::_helper_test", builder_publish, "nmg",
	    "command-result", "nmg::_helper_test", "nmg-test", 0) != BRLCAD_OK) {
	bg_line_layer_builder_free(builder_publish);
	FAIL("line-layer builder helper should publish through command-scene ownership");
    }
    bg_line_layer_builder_free(builder_publish);
    BRLObolFeatureHandle builder_handle =
	owned_controller->features().find("nmg::_helper_test",
		BRLOBOL_FEATURE_SCOPE_SHARED);
    BRLObolFeatureRecord builder_record;
    if (!builder_handle.isValid() ||
	    !owned_controller->features().record(builder_handle,
		builder_record) ||
	    builder_record.kind != BRLObolFeatureKind::LineLayer ||
	    builder_record.scope != BRLObolFeatureScope::Shared ||
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
		BRLObolOverlayClass::CommandResult,
		BRLObolOverlayLifecycle::PerCommand,
		BRLObolOverlayOrder::PostTransparent))
	FAIL("line-layer builder helper should preserve shared owned command-result metadata");
    command_scene_desc.owner_id = "nmg";
    command_scene = ged_draw_command_scene_begin(feature_view_ctx,
	    &command_scene_desc);
    if (!command_scene ||
	    ged_draw_command_scene_features_remove_prefix(command_scene,
		"nmg::_helper_test") != 1 ||
	    !ged_draw_command_scene_commit(command_scene) ||
	    owned_controller->features().exists("nmg::_helper_test"))
	FAIL("line-layer builder helper result should clean up by owner-scoped prefix");

    struct bg_line_layer_builder *builder_generation =
	bg_line_layer_builder_create();
    if (!builder_generation)
	FAIL("command-scene builder generation test should create builder");
    point_t builder_latest_b = {0.0, 2.0, 0.0};
    if (!bg_line_layer_builder_add(builder_generation, 255, 255, 0,
		builder_a, BG_GEOMETRY_LINE_MOVE) ||
	    !bg_line_layer_builder_add(builder_generation, 255, 255, 0,
		builder_latest_b, BG_GEOMETRY_LINE_DRAW)) {
	bg_line_layer_builder_free(builder_generation);
	FAIL("command-scene builder generation test should accept latest geometry");
    }
    if (_ged_line_layer_builder_publish_command_scene_feature(gedp,
	    "nmg::_helper_generation", builder_generation, "nmg",
	    "command-result", "nmg::_helper_generation", "nmg-generation",
	    22) != BRLCAD_OK) {
	bg_line_layer_builder_free(builder_generation);
	FAIL("line-layer builder helper should publish latest generation");
    }
    bg_line_layer_builder_free(builder_generation);

    builder_generation = bg_line_layer_builder_create();
    if (!builder_generation)
	FAIL("command-scene builder stale generation test should create builder");
    point_t builder_stale_b = {0.0, 3.0, 0.0};
    if (!bg_line_layer_builder_add(builder_generation, 255, 255, 0,
		builder_a, BG_GEOMETRY_LINE_MOVE) ||
	    !bg_line_layer_builder_add(builder_generation, 255, 255, 0,
		builder_stale_b, BG_GEOMETRY_LINE_DRAW)) {
	bg_line_layer_builder_free(builder_generation);
	FAIL("command-scene builder stale generation test should accept stale geometry");
    }
    if (_ged_line_layer_builder_publish_command_scene_feature(gedp,
	    "nmg::_helper_generation", builder_generation, "nmg",
	    "command-result", "nmg::_helper_generation", "nmg-generation",
	    21) == BRLCAD_OK) {
	bg_line_layer_builder_free(builder_generation);
	FAIL("line-layer builder helper should reject stale generation without fallback");
    }
    bg_line_layer_builder_free(builder_generation);

    builder_handle = owned_controller->features().find(
	    "nmg::_helper_generation", BRLOBOL_FEATURE_SCOPE_SHARED);
    if (!builder_handle.isValid() ||
	    !owned_controller->features().record(builder_handle,
		builder_record) ||
	    builder_record.owner.generation != 22 ||
	    builder_record.points.size() != 2 ||
	    builder_record.points[1][1] != 2.0f ||
	    builder_record.metadata.size() < 7 ||
	    !BU_STR_EQUAL(builder_record.metadata[5].value.getString(),
		"nmg-generation") ||
	    !BU_STR_EQUAL(builder_record.metadata[6].key.getString(),
		"result.generation") ||
	    !BU_STR_EQUAL(builder_record.metadata[6].value.getString(), "22"))
	FAIL("line-layer builder helper stale generation should not replace latest result");
    command_scene_desc.owner_id = "nmg";
    command_scene_desc.generation = 22;
    command_scene = ged_draw_command_scene_begin(feature_view_ctx,
	    &command_scene_desc);
    if (!command_scene ||
	    ged_draw_command_scene_features_remove_prefix(command_scene,
		"nmg::_helper_generation") != 1 ||
	    !ged_draw_command_scene_commit(command_scene) ||
	    owned_controller->features().exists("nmg::_helper_generation"))
	FAIL("line-layer builder helper generation result should clean up by owner-scoped prefix");
    command_scene_desc.generation = 0;

    struct rt_preview_callback_state ged_preview_state = {77, 0, 0};
    struct ged_draw_view_edit_preview_callbacks ged_preview_callbacks =
	GED_DRAW_VIEW_EDIT_PREVIEW_CALLBACKS_INIT;
    ged_preview_callbacks.revision_cb = rt_preview_revision_cb;
    ged_preview_callbacks.update_cb = rt_preview_update_cb;
    ged_preview_callbacks.pick_cb = rt_preview_pick_cb;
    int ged_preview_owner = 0;
    ged_draw_view_feature_ref ged_preview_ref =
	ged_draw_view_context_feature_overlay_ensure(feature_view_ctx,
		"cap2::ged-preview", &ged_preview_owner, &ged_preview_state,
		&ged_preview_callbacks, "cap2::ged-source.s");
    if (ged_draw_view_feature_ref_is_null(ged_preview_ref))
	FAIL("GED feature overlay ensure should return an Obol feature ref");
    BRLObolFeatureHandle ged_preview_handle =
	owned_controller->features().find("cap2::ged-preview");
    BRLObolFeatureSummary ged_preview_summary;
    if (!ged_preview_handle.isValid() ||
	    !owned_controller->features().summary("cap2::ged-preview",
		ged_preview_summary) ||
	    !ged_preview_summary.exists ||
	    ged_preview_summary.kind != BRLObolFeatureKind::EditPreview ||
	    !ged_preview_summary.overlay.isOverlay ||
	    ged_preview_summary.overlay.overlayClass !=
		BRLObolOverlayClass::EditHandle ||
	    ged_preview_summary.overlay.lifecycle !=
		BRLObolOverlayLifecycle::PerTool ||
	    ged_preview_summary.overlay.order !=
		BRLObolOverlayOrder::PostTransparent ||
	    ged_preview_summary.overlay.ownerToken != &ged_preview_owner ||
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
    if (!ged_draw_view_feature_points_replace(ged_preview_ref,
		GED_DRAW_VIEW_FEATURE_TRANSIENT_PREVIEW, ged_preview_points,
		ged_preview_cmds, 3))
	FAIL("GED feature points replacement should update Obol edit-preview geometry");
    BRLObolFeatureRecord ged_preview_record;
    if (!owned_controller->features().record(ged_preview_handle,
		ged_preview_record) ||
	    ged_preview_record.kind != BRLObolFeatureKind::EditPreview ||
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
    if (!ged_draw_view_feature_primitive_wireframe_replace(ged_preview_ref,
		gedp->dbip, &preview_box_intern, preview_xform, NULL, NULL)) {
	rt_db_free_internal(&preview_box_intern);
	FAIL("GED feature primitive wireframe helper should publish transformed primitive geometry");
    }
    rt_db_free_internal(&preview_box_intern);
    if (!owned_controller->features().record(ged_preview_handle,
		ged_preview_record) ||
	    ged_preview_record.kind != BRLObolFeatureKind::EditPreview ||
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
    ged_draw_view_feature_set_visible(ged_preview_ref, 0);
    ged_draw_view_feature_set_color(ged_preview_ref, 12, 34, 56);
    BRLObolFeatureStyle ged_preview_style;
    if (!owned_controller->features().style(ged_preview_handle,
		ged_preview_style) ||
	    !ged_preview_style.hasVisible ||
	    ged_preview_style.visible ||
	    !ged_preview_style.hasColor)
	FAIL("GED feature style mutations should update Obol feature style");
    if (!ged_draw_view_feature_touch(ged_preview_ref) ||
	    ged_preview_state.update_count != 1)
	FAIL("GED feature touch should dispatch through Obol edit-preview callbacks");
    if (!ged_draw_view_context_edit_preview_publish_event(feature_view_ctx,
		ged_preview_ref, GED_DRAW_VIEW_EDIT_PREVIEW_UPDATE,
		"cap2::ged-source.s"))
	FAIL("GED feature preview events should route through the Obol feature API");
    if (!ged_draw_view_feature_clear_geometry(ged_preview_ref) ||
	    !owned_controller->features().record(ged_preview_handle,
		ged_preview_record) ||
	    !ged_preview_record.points.empty())
	FAIL("GED feature clear geometry should clear the Obol feature record");

    int ged_label_owner = 0;
    ged_draw_view_feature_ref ged_label_ref =
	ged_draw_view_context_feature_label_ensure(feature_view_ctx,
		"cap2::ged-label", &ged_label_owner);
    struct ged_draw_view_feature_label ged_label;
    memset(&ged_label, 0, sizeof(ged_label));
    ged_label.text = "ged label";
    VSET(ged_label.point, 3.0, 4.0, 5.0);
    ged_label.color_valid = 1;
    ged_label.color[0] = 210;
    ged_label.color[1] = 211;
    ged_label.color[2] = 212;
    ged_label.font_size = 16.0;
    if (ged_draw_view_feature_ref_is_null(ged_label_ref) ||
	    !ged_draw_view_feature_labels_replace(ged_label_ref, &ged_label, 1))
	FAIL("GED feature label replacement should route into Obol labels");
    BRLObolFeatureHandle ged_label_handle =
	owned_controller->features().find("cap2::ged-label");
    BRLObolFeatureRecord ged_label_record;
    if (!ged_label_handle.isValid() ||
	    !owned_controller->features().record(ged_label_handle,
		ged_label_record) ||
	    ged_label_record.kind != BRLObolFeatureKind::Labels ||
	    ged_label_record.labels.size() != 1 ||
	    !BU_STR_EQUAL(ged_label_record.labels[0].text.getString(),
		"ged label") ||
	    fabs(ged_label_record.labels[0].fontSize - 16.0f) > 0.001f ||
	    !ged_label_record.overlay.isOverlay ||
	    ged_label_record.overlay.ownerToken != &ged_label_owner)
	FAIL("GED feature label API should publish typed Obol label records");
    if (!ged_draw_view_context_feature_remove(feature_view_ctx, "cap2::ged-label") ||
	    owned_controller->features().exists("cap2::ged-label"))
	FAIL("GED feature remove should delete Obol-backed feature records");

    if (ged_draw_view_context_features_remove_prefix(feature_view_ctx,
	    "cap2::") < 11 ||
	    owned_controller->features().exists("cap2::line") ||
	    owned_controller->features().exists("cap2::tcl-line") ||
	    owned_controller->features().exists("cap2::annotation") ||
	    owned_controller->features().exists("cap2::polygon-overlay") ||
	    owned_controller->features().exists("cap2::label") ||
	    owned_controller->features().exists("cap2::arrow") ||
	    owned_controller->features().exists("cap2::axes") ||
	    owned_controller->features().exists("cap2::tcl-axes") ||
	    owned_controller->features().exists("cap2::mesh") ||
	    owned_controller->features().exists("cap2::diagnostic") ||
	    owned_controller->features().exists("cap2::rt-preview"))
	FAIL("GED feature prefix removal should clear owned Obol feature store entries");

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
		drawMode.getValue() != BRLOBOL_LOD_DRAW_SHADED)
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
		drawMode.getValue() != BRLOBOL_LOD_DRAW_SHADED)
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
	    SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE, 1, 0, 1,
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
	    SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS, 0, 1, 1,
	    "evaluated-points"))
	return 1;
    const char *draw_eval_points_option[4] = {"draw",
	"--evaluated-points", "--add-mode", "box.s"};
    if (ged_exec_draw(gedp, 4, draw_eval_points_option) != BRLCAD_OK)
	FAIL("GED evaluated-points option draw should succeed");
    if (source_representation_count(owned_scene, "box.s",
	    SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS) != 1)
	FAIL("GED evaluated-points option should create one mode source");
    (void)try_path_transaction(gedp, GED_DRAW_TXN_ERASE, "box.s",
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
    SoBRLVListShape *annot_shape = annot_source ?
	annot_source->getRealizedShape() : NULL;
    const SoBRLVListShape *annot_geom_source = annot_shape ?
	annot_shape->getGeometrySource() : NULL;
    if (owned_scene->getDatabaseSourceCount() != 3 ||
	    !annot_shape || !annot_geom_source ||
	    annot_geom_source->point.getNum() != 2 ||
	    annot_geom_source->command.getNum() != 2 ||
	    annot_geom_source->command[0] != SoBRLVListShape::MOVE ||
	    annot_geom_source->command[1] != SoBRLVListShape::DRAW ||
	    strcmp(annot_shape->sourceType.getValue().getString(),
		"annotation") != 0 ||
	    strcmp(annot_shape->geometryKind.getValue().getString(),
		"annotation") != 0 ||
	    fabs(annot_geom_source->point[1][0] - 50.25f) > 0.001f ||
	    fabs(annot_geom_source->point[1][1] - 0.5f) > 0.001f ||
	    fabs(annot_geom_source->point[1][2]) > 0.001f)
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
    SoBRLVListShape *submodel_shape = submodel_source->getRealizedShape();
    const SoBRLVListShape *submodel_geom = submodel_shape ?
	submodel_shape->getGeometrySource() : NULL;
    if (!submodel_shape || !submodel_geom ||
	    submodel_geom->point.getNum() == 0 ||
	    submodel_geom->command.getNum() != submodel_geom->point.getNum() ||
	    strcmp(submodel_shape->recordRole.getValue().getString(),
		"database") != 0 ||
	    strcmp(submodel_shape->sourceType.getValue().getString(),
		"submodel") != 0 ||
	    !path_equal(submodel_shape->ownerSourcePath.getValue().getString(),
		"submodel_owner.s") ||
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
    SoBRLVListShape *submodel_temp_shape =
	submodel_temp_source->getRealizedShape();
    const SoBRLVListShape *submodel_temp_geom = submodel_temp_shape ?
	submodel_temp_shape->getGeometrySource() : NULL;
    if (!submodel_temp_shape || !submodel_temp_geom ||
	    submodel_temp_geom->point.getNum() == 0 ||
	    strcmp(submodel_temp_shape->recordRole.getValue().getString(),
		"database") != 0 ||
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
    void *mesh_lod_view_ctx = ged_draw_shape_ref_view_context(gedp,
	    mesh_record.ref);
    if (!mesh_lod_view_ctx)
	FAIL("GED BoT mesh LoD view context should be available");
    void *mesh_lod_view_ctxs[1] = {mesh_lod_view_ctx};
    if (!ged_draw_shape_ref_lod_ensure(gedp, mesh_record.ref,
	    mesh_lod_view_ctx, mesh_lod_view_ctxs, 1))
	FAIL("GED BoT mesh LoD ensure should succeed for owned Obol mesh update");
    SoBRLDatabaseSource *mesh_source =
	source_for_path(owned_scene, "mesh_owner.bot");
    SoBRLMeshShape *mesh_shape = mesh_source ?
	mesh_source->getRealizedMesh() : NULL;
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
    const char *erase_mesh_owner[2] = {"erase", "mesh_owner.bot"};
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
	    strcmp(brep_wire_shape->sourceType.getValue().getString(),
		"line-set") != 0)
	FAIL("GED BREP wireframe draw should publish owned Obol line geometry");
    const char *erase_brep_wire[2] = {"erase", "brep_owner.brep"};
    if (ged_exec_erase(gedp, 2, erase_brep_wire) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "brep_owner.brep"))
	FAIL("GED BREP wireframe erase should restore the owned Obol source baseline");
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
    void *brep_lod_view_ctx = ged_draw_shape_ref_view_context(gedp,
	    brep_record.ref);
    if (!brep_lod_view_ctx)
	FAIL("GED BREP mesh LoD view context should be available");
    void *brep_lod_view_ctxs[1] = {brep_lod_view_ctx};
    if (!ged_draw_shape_ref_lod_ensure(gedp, brep_record.ref,
	    brep_lod_view_ctx, brep_lod_view_ctxs, 1))
	FAIL("GED BREP mesh LoD ensure should succeed for owned Obol mesh update");
    SoBRLDatabaseSource *brep_source =
	source_for_path(owned_scene, "brep_owner.brep");
    SoBRLMeshShape *brep_shape = brep_source ?
	brep_source->getRealizedMesh() : NULL;
    SbVec3f brep_lod_bmin;
    SbVec3f brep_lod_bmax;
    if (!brep_source || !brep_shape ||
	    brep_shape->point.getNum() == 0 ||
	    brep_shape->coordIndex.getNum() == 0)
	FAIL("GED BREP mesh LoD update should publish owned Obol mesh fields");
    if (!brep_source->getMeshLod())
	FAIL("GED Obol BREP mesh LoD draw should store the runtime handle on the owned source");
    if (!brep_source->getMeshLodBounds(brep_lod_bmin, brep_lod_bmax))
	FAIL("GED Obol BREP mesh LoD draw should store runtime bounds on the owned source");
    const char *erase_brep_owner[2] = {"erase", "brep_owner.brep"};
    if (ged_exec_erase(gedp, 2, erase_brep_owner) != BRLCAD_OK ||
	    owned_scene->getDatabaseSourceCount() != 2 ||
	    source_for_path(owned_scene, "brep_owner.brep"))
	FAIL("GED BREP mesh LoD erase should restore the owned Obol source baseline");
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
    void *group_only_ctx = ged_draw_group_ref_context(gedp, group_only_ref);
    struct ged_draw_scene_tree_summary group_only_context_tree;
    memset(&group_only_context_tree, 0, sizeof(group_only_context_tree));
    if (!group_only_ctx ||
	    !ged_draw_scene_context_tree_summary(group_only_ctx,
		&group_only_context_tree) ||
	    group_only_context_tree.child_count != 1)
	FAIL("GED group-ref context should resolve to an owned Obol group context");
    void *group_only_parent_ctx =
	ged_draw_scene_context_parent(group_only_ctx);
    struct ged_draw_scene_tree_summary group_only_parent_tree;
    memset(&group_only_parent_tree, 0, sizeof(group_only_parent_tree));
    if (!group_only_parent_ctx ||
	    !ged_draw_scene_context_tree_summary(group_only_parent_ctx,
		&group_only_parent_tree) ||
	    !group_only_parent_tree.valid ||
	    !group_only_parent_tree.is_group ||
	    group_only_parent_tree.has_parent ||
	    !path_equal(ged_draw_scene_context_name(group_only_parent_ctx),
		"/"))
	FAIL("GED group-ref context should expose an owned Obol parent context");
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
    ged_draw_erase_name(gedp, "group_only.s");
    memset(&group_only_tree, 0, sizeof(group_only_tree));
    if (!ged_draw_group_ref_tree_summary(gedp, group_only_ref,
	    &group_only_tree) ||
	    !group_only_tree.valid ||
	    !group_only_tree.is_group)
	FAIL("GED public name erase should preserve overlay groups from owned Obol state");
    if (owned_scene->setGroupDrawIntent("group_only.s",
	    "ged-draw-group:group_only.s", BRLOBOL_LOD_DRAW_WIRE,
	    BRLOBOL_LOD_DRAW_WIRE, FALSE, 0) < 0)
	FAIL("owned Obol group overlay erase sentinel restore should succeed");
    struct ged_draw_scene_tree_summary original_root_count_tree;
    memset(&original_root_count_tree, 0, sizeof(original_root_count_tree));
    if (!ged_draw_scene_context_tree_summary(group_only_parent_ctx,
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
	    !ged_draw_scene_context_tree_summary(group_only_parent_ctx,
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
    void *box_registry_ctx = ged_draw_scene_handle_context(
	    ged_draw_registry_shape_ref_scene_handle(gedp, box_record.ref));
    struct ged_draw_scene_tree_summary box_registry_context_tree;
    memset(&box_registry_context_tree, 0, sizeof(box_registry_context_tree));
    if (!box_registry_ctx ||
	    !ged_draw_scene_context_tree_summary(box_registry_ctx,
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
    void *box_registry_parent_ctx =
	ged_draw_scene_context_parent(box_registry_ctx);
    struct ged_draw_scene_tree_summary box_registry_parent_tree;
    memset(&box_registry_parent_tree, 0,
	    sizeof(box_registry_parent_tree));
    if (!box_registry_parent_ctx ||
	    !ged_draw_scene_context_tree_summary(box_registry_parent_ctx,
		&box_registry_parent_tree) ||
	    !box_registry_parent_tree.valid ||
	    !box_registry_parent_tree.is_group ||
	    box_registry_parent_tree.is_shape ||
	    box_registry_parent_tree.has_parent ||
	    !path_equal(ged_draw_scene_context_name(
		    box_registry_parent_ctx), "/"))
	FAIL("GED registry semantic source parents should resolve to owned Obol parent contexts");
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
	    if (ged_draw_bounds(gedp, &draw_bounds_min, &draw_bounds_max, 0) ||
		    draw_bounds_max[0] < 22.9)
		FAIL("GED draw bounds should reflect published owned Obol VLIST points");
	    point_t explicit_center;
	    VSET(explicit_center, 30.0, 31.0, 32.0);
	    ged_draw_index_stats_reset(gedp);
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
		    strcmp(aux_shape->recordRole.getValue().getString(),
			"auxiliary") != 0 ||
		    aux_shape->drawMode.getValue() !=
			BRLOBOL_LOD_DRAW_SHADED ||
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
	    if (!ged_draw_shape_ref_geometry_clear(gedp, box_record.ref) ||
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
	    if (!ged_draw_obol_database_source_publish_indexed_face_set_for_path(
		    gedp, "box.s", (const point_t *)mesh_points, 4,
		    (const vect_t *)mesh_normals, 4, mesh_indices, 5))
		FAIL("GED Obol source indexed-face publish should succeed");
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
	    box_source = source_for_representation(owned_scene, "box.s",
		    SoBRLDatabaseSource::REPRESENTATION_SHADED);
	    if (!box_source) {
		box_source = owned_scene->findDatabaseSource("/box.s");
		if (!box_source)
		    box_source = source_for_path(owned_scene, "box.s");
	    }
	    if (!box_source)
		FAIL("owned Obol source should remain available after shaded draw source realization");
	    box_mesh = box_source->getRealizedMesh();
	    {
		const SoBRLMeshShape *box_mesh_geom = box_mesh ?
		    box_mesh->getGeometrySource() : NULL;
		if (!box_mesh_geom ||
			box_mesh_geom->point.getNum() == 0 ||
			box_mesh_geom->coordIndex.getNum() == 0)
		FAIL("GED shaded source realization should publish owned Obol mesh fields");
	    }
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
	    if (ged_draw_bounds(gedp, &draw_bounds_min, &draw_bounds_max, 0) ||
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
	    void *lod_view_ctx = ged_draw_shape_ref_view_context(gedp,
		    box_record.ref);
	    if (!lod_view_ctx)
		FAIL("GED LoD view context should be available");
	    if (!bv_scale_set(DRAW_TEST_BV(lod_view_ctx), 7.0))
		FAIL("GED LoD view context scale should be settable");
	    ged_draw_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
	    lod_policy.csg_enabled = 1;
	    lod_policy.mesh_enabled = 0;
	    lod_policy.scale = 1.75;
	    lod_policy.bot_threshold = 77;
	    lod_policy.curve_scale = 2.25;
	    lod_policy.point_scale = 3.25;
	    if (!ged_draw_view_context_lod_policy_apply(lod_view_ctx, &lod_policy))
		FAIL("GED LoD view policy should be settable");
	    void *lod_view_ctxs[1] = {lod_view_ctx};
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
		    if (!ged_draw_shape_ref_geometry_clear(gedp, box_record.ref))
			FAIL("GED shape geometry clear should clear shaded source realization mesh");
		    ged_draw_index_stats_reset(gedp);
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
    if (!bv_scale_set(DRAW_TEST_BV(stale_lod_view_ctx), 9.0))
	FAIL("GED stale shape-ref LoD context scale should be settable");
    ged_draw_view_lod_policy stale_lod_policy = BV_LOD_POLICY_INIT;
    stale_lod_policy.csg_enabled = 1;
    stale_lod_policy.mesh_enabled = 0;
    stale_lod_policy.scale = 2.75;
    stale_lod_policy.bot_threshold = 91;
    stale_lod_policy.curve_scale = 6.25;
    stale_lod_policy.point_scale = 7.25;
    if (!ged_draw_view_context_lod_policy_apply(stale_lod_view_ctx,
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
		BRLOBOL_LOD_DRAW_SHADED)
	FAIL("GED group mode setter should mutate the owned Obol group draw intent");
    if (!ged_draw_group_ref_set_mode(gedp, group_state.ref,
	    GED_DRAW_MODE_WIRE))
	FAIL("GED group mode restore should succeed");
    std::string obol_group_intent_path =
	std::string("ged-draw-group:") + group_path + "_intent";
    if (owned_scene->setGroupDrawIntent(group_path.c_str(),
	    obol_group_intent_path.c_str(), BRLOBOL_LOD_DRAW_SHADED,
	    BRLOBOL_LOD_DRAW_WIRE, TRUE, 501) <= 0)
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
	    original_obol_group_intent_path.c_str(), BRLOBOL_LOD_DRAW_WIRE,
	    BRLOBOL_LOD_DRAW_WIRE, FALSE, 0) <= 0)
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
    struct ged_draw_transaction obol_index_txn =
	ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED,
		"draft_move.s");
    obol_index_txn.redraw = 0;
    struct ged_draw_transaction_result obol_index_result;
    ged_draw_transaction_result_init(&obol_index_result);
    if (ged_draw_apply_transaction(gedp, &obol_index_txn,
	    &obol_index_result) <= 0)
	FAIL("GED Obol component-index stale transaction should succeed");
    ged_draw_transaction_result_free(&obol_index_result);
    struct ged_draw_index_stats obol_index_stats;
    memset(&obol_index_stats, 0, sizeof(obol_index_stats));
    ged_draw_index_stats_get(gedp, &obol_index_stats);
    if (obol_index_stats.slow_path_shape_scans ||
	    obol_index_stats.slow_path_group_scans)
	FAIL("GED Obol component indexes should avoid registry/index slow-path scans");
    SoBRLDatabaseSource *draft_move_source =
	source_for_path(owned_scene, "draft_move.s");
    BRLObolDatabaseSourceSummary draft_move_summary;
    if (!draft_move_source ||
	    !draft_move_source->getSummary(draft_move_summary) ||
	    !draft_move_summary.stale)
	FAIL("GED Obol component-index transaction should mark the owned source stale");
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
	FAIL("GED redraw transaction should refresh draw sources without clearing Obol-only sources");
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
	GED_DRAW_GROUP_REF_NULL_INIT, 0, 0, "ball.s", 0, 0, 0, 0, 0.0,
	0};
    ged_draw_foreach_shape_record(gedp, record_source_state_cb,
	    &ball_record);
    if (!ball_record.found)
	FAIL("GED ball shape record should be available before public erase");
    if (apply_path_transaction(gedp, GED_DRAW_TXN_ERASE, "ball.s", NULL, -1,
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
    if (apply_path_transaction(gedp, GED_DRAW_TXN_ERASE, "ball.s",
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
    if (apply_path_transaction(gedp, GED_DRAW_TXN_ERASE, "ball.s", NULL, -1,
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
    SoBRLDatabaseSource *nested_leaf_source =
	source_for_path(owned_scene, nested_leaf_source_path);
    SoBRLDatabaseSource *nested_sibling_source =
	source_for_path(owned_scene, nested_sibling_source_path);
    if (!nested_leaf_source || !nested_sibling_source ||
	    owned_scene->getDatabaseSourceCount() != 4)
	FAIL("GED nested path draws should create owned Obol child sources");
    BRLObolDatabaseSourceSummary nested_sibling_initial_summary;
    BRLObolDatabaseSourceSummary nested_sibling_summary;
    if (!nested_sibling_source->getSummary(nested_sibling_initial_summary))
	FAIL("GED nested sibling source should expose a state summary");
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
    if (apply_path_transaction(gedp, GED_DRAW_TXN_ERASE_PREFIX,
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
    if (apply_path_transaction(gedp, GED_DRAW_TXN_ERASE_PREFIX,
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
    if (apply_path_transaction(gedp, GED_DRAW_TXN_ERASE_PREFIX,
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
    struct ged_draw_transaction reexpand_nested_child =
	ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED,
		"nested_child.c");
    reexpand_nested_child.redraw = 1;
    struct ged_draw_transaction_result nested_result;
    ged_draw_transaction_result_init(&nested_result);
    ged_draw_index_stats_reset(gedp);
    if (ged_draw_apply_transaction(gedp, &reexpand_nested_child,
	    &nested_result) <= 0)
	FAIL("GED nested child reexpand transaction should succeed");
    ged_draw_transaction_result_free(&nested_result);
    struct ged_draw_index_stats nested_reexpand_stats;
    memset(&nested_reexpand_stats, 0, sizeof(nested_reexpand_stats));
    ged_draw_index_stats_get(gedp, &nested_reexpand_stats);
    if (nested_reexpand_stats.slow_path_shape_scans ||
	    nested_reexpand_stats.slow_path_group_scans)
	FAIL("GED nested child reexpand should avoid registry/index slow-path scans");
    struct ged_draw_transaction stale_nested_leaf =
	ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED,
		"nested_leaf.s");
    stale_nested_leaf.redraw = 0;
    ged_draw_transaction_result_init(&nested_result);
    if (ged_draw_apply_transaction(gedp, &stale_nested_leaf,
	    &nested_result) <= 0)
	FAIL("GED nested leaf stale transaction should succeed");
    ged_draw_transaction_result_free(&nested_result);
    BRLObolDatabaseSourceSummary nested_leaf_summary;
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
    ged_draw_index_stats_reset(gedp);
    if (ged_draw_apply_transaction(gedp, &redraw_nested_leaf,
	    &nested_result) <= 0)
	FAIL("GED nested leaf redraw transaction should succeed");
    ged_draw_transaction_result_free(&nested_result);
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
	    !nested_leaf_source->getSummary(nested_leaf_summary) ||
	    nested_leaf_summary.stale ||
	    !nested_sibling_source ||
	    !nested_sibling_source->getSummary(nested_sibling_summary) ||
	    nested_sibling_summary.lineWidth != 23)
	FAIL("GED nested leaf redraw transaction should preserve unrelated owned Obol source state");
    if (apply_path_transaction(gedp, GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED,
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
    if (apply_path_transaction(gedp, GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED,
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
    BRLObolDatabaseSourceSummary component_mode_summary;
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
    void *nested_child_ctx = ged_draw_group_ref_context(gedp,
	    nested_child_group);
    struct db_full_path nested_child_ctx_path;
    db_full_path_init(&nested_child_ctx_path);
    if (!nested_child_ctx ||
	    ged_draw_group_context_dbpath(gedp, nested_child_ctx,
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
    if (apply_path_transaction(gedp, GED_DRAW_TXN_ERASE,
	    "nested_parent.c/nested_child_renamed.c",
	    ged_draw_active_view_ctx(gedp), -1,
	    "public active-scope nested group path erase"))
	return 1;
    if (!owned_scene->findGroup("nested_parent.c") ||
	    owned_scene->findGroup("nested_parent.c/nested_child_renamed.c"))
		FAIL("GED public active-scope nested group path erase should remove the owned Obol nested group");
    if (apply_path_transaction(gedp, GED_DRAW_TXN_ERASE, "nested_parent.c",
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
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK)
	FAIL("owned-controller direct-clear redraw commands should succeed");
    if (owned_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(owned_scene, "box.s") ||
	    !source_for_path(owned_scene, "ball.s"))
	FAIL("owned Obol scene controller should mirror redraw after direct clear");
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
    if (exercise_mesh_source_local_publication(gedp, owned_scene,
	    "box.s"))
	return 1;
    if (exercise_multi_instance_transform_reuse(gedp, owned_scene))
	return 1;

    ged_draw_obol_scene_controller_detach(gedp);
    if (ged_draw_obol_scene_controller(gedp) ||
	    ged_draw_obol_controller(gedp) ||
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

    if (!seed_view_lod_probe_payload(&view_controller, "box.s", "box.s"))
	FAIL("attached Obol view-controller LoD invalidation probe should seed draw payload");
    const char *attached_draw_draft[2] = {"draw", "draft_move.s"};
    if (ged_exec_draw(gedp, 2, attached_draw_draft) != BRLCAD_OK ||
	    view_controller.getViewLodState()->payloadCount() != 0)
	FAIL("attached Obol draw transaction should clear view-local LoD state");
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
    struct ged_draw_transaction attached_source_update =
	ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED, "box.s");
    attached_source_update.redraw = 0;
    struct ged_draw_transaction_result attached_source_result;
    ged_draw_transaction_result_init(&attached_source_result);
    if (ged_draw_apply_transaction(gedp, &attached_source_update,
	    &attached_source_result) <= 0 ||
	    view_controller.getViewLodState()->payloadCount() != 0) {
	ged_draw_transaction_result_free(&attached_source_result);
	FAIL("attached Obol source-update transaction should clear view-local LoD state");
    }
    ged_draw_transaction_result_free(&attached_source_result);
    const char *attached_redraw_box[2] = {"draw", "box.s"};
    if (ged_exec_draw(gedp, 2, attached_redraw_box) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(view_scene, "box.s"))
	FAIL("attached Obol source-update refresh should restore box source");

    if (!seed_view_lod_probe_payload(&view_controller, "box.s", "box.s"))
	FAIL("attached Obol view-controller LoD invalidation probe should seed clear payload");
    struct ged_draw_transaction attached_clear =
	ged_draw_transaction_make(GED_DRAW_TXN_CLEAR, NULL);
    struct ged_draw_transaction_result attached_clear_result;
    ged_draw_transaction_result_init(&attached_clear_result);
    if (ged_draw_apply_transaction(gedp, &attached_clear,
	    &attached_clear_result) <= 0) {
	ged_draw_transaction_result_free(&attached_clear_result);
	FAIL("attached Obol clear transaction should succeed");
    }
    ged_draw_transaction_result_free(&attached_clear_result);
    if (view_controller.getViewLodState()->payloadCount() != 0 ||
	    view_scene->getDatabaseSourceCount() != 0)
	FAIL("attached Obol clear transaction should clear view-local LoD state and scene sources");
    if (ged_exec_draw(gedp, 2, draw_box) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 2, draw_ball) != BRLCAD_OK ||
	    view_scene->getDatabaseSourceCount() != 2 ||
	    !source_for_path(view_scene, "box.s") ||
	    !source_for_path(view_scene, "ball.s"))
	FAIL("attached Obol clear transaction redraw should restore current GED draw state");

    struct ged_draw_transaction attached_stale_source =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "box.s");
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
    struct ged_draw_transaction attached_erase_prefix =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE_PREFIX,
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
    struct ged_draw_transaction attached_reference_removal =
	ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED,
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

    struct ged_draw_transaction attached_clear_scope =
	ged_draw_transaction_make(GED_DRAW_TXN_CLEAR_SCOPE, NULL);
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

    struct ged_draw_transaction attached_teardown =
	ged_draw_transaction_make(GED_DRAW_TXN_TEARDOWN, NULL);
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

    void *independent_view = ged_view_context_create();
    if (!independent_view)
	FAIL("Obol independent view source-owner test view should be created");
    if (!bv_name_set(DRAW_TEST_BV(independent_view), "V0"))
	FAIL("Obol independent view source-owner test view should be named");
    if (!ged_view_set_context_add(ged_view_set_ctx(gedp), independent_view))
	FAIL("Obol independent view source-owner test view should be registered");
    ged_view_context_owned_add(gedp, independent_view);

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

    if (view_scene->getDatabaseSourceCount() != 4 ||
	    !source_for_shared_path(view_scene, "box.s") ||
	    !source_for_shared_path(view_scene, "ball.s") ||
	    !source_for_view_path(view_scene, "V0", "box.s") ||
	    !source_for_view_path(view_scene, "V0", "ball.s"))
	FAIL("Obol independent view setup should create scoped source owners without replacing shared owners");

    const char *erase_v0_box[5] = {"erase", "-V", "V0", "box.s", NULL};
    if (ged_exec_erase(gedp, 4, erase_v0_box) != BRLCAD_OK)
	FAIL("real GED independent-view erase command should succeed");
    if (view_scene->getDatabaseSourceCount() != 3 ||
	    source_for_view_path(view_scene, "V0", "box.s") ||
	    !source_for_view_path(view_scene, "V0", "ball.s") ||
	    !source_for_shared_path(view_scene, "box.s") ||
	    !source_for_shared_path(view_scene, "ball.s"))
	FAIL("Obol independent-view erase should remove only the scoped source owner");

    if (ged_exec_draw(gedp, 5, draw_v0_box) != BRLCAD_OK)
	FAIL("real GED independent-view redraw command should succeed");
    if (view_scene->getDatabaseSourceCount() != 4 ||
	    !source_for_view_path(view_scene, "V0", "box.s") ||
	    !source_for_view_path(view_scene, "V0", "ball.s") ||
	    !source_for_shared_path(view_scene, "box.s") ||
	    !source_for_shared_path(view_scene, "ball.s"))
	FAIL("Obol independent-view draw should restore only the scoped source owner");

    const char *zap_v0[5] = {"zap", "-V", "V0", "-g", NULL};
    if (ged_exec_zap(gedp, 4, zap_v0) != BRLCAD_OK)
	FAIL("real GED independent-view zap command should succeed");
    if (view_scene->getDatabaseSourceCount() != 2 ||
	    source_for_view_path(view_scene, "V0", "box.s") ||
	    source_for_view_path(view_scene, "V0", "ball.s") ||
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
    void *zap_view_ctx = ged_view_active_ctx(gedp);
    ged_draw_view_polygon_ref zap_created =
	ged_draw_view_context_polygon_create(zap_view_ctx,
		"zap::polygon", 0, GED_DRAW_VIEW_POLYGON_SQUARE,
		zap_polygon_origin);
    ged_draw_view_polygon_ref zap_found =
	ged_draw_view_context_polygon_find(zap_view_ctx, "zap::polygon");
    if (ged_draw_view_polygon_ref_is_null(zap_created) ||
	    ged_draw_view_polygon_ref_is_null(zap_found))
	FAIL("Obol polygon store should publish a test polygon before zap");

    const char *zap_cmd[1] = {"zap"};
    if (ged_exec_zap(gedp, 1, zap_cmd) != BRLCAD_OK)
	FAIL("real GED zap command should succeed");
    if (view_scene->getDatabaseSourceCount() != 0)
	FAIL("attached Obol scene controller should mirror GED clear command");
    if (!ged_draw_view_polygon_ref_is_null(
		ged_draw_view_context_polygon_find(zap_view_ctx,
		    "zap::polygon")))
	FAIL("GED zap should clear Obol polygon store view features");

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
