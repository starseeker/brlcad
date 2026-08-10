/*        T E S T _ Q T C A D _ O B O L _ E D I T _ P R E V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "bu/app.h"

#include "bu/str.h"

#include "bv.h"

#include "BObol/BEditPreview.h"
#include "BObol/BExportAction.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BSnapAction.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "BObol/BVListShape.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "ged.h"
#include "ged/draw.h"
#include "ged/view.h"
#include "ged/view_feature_internal.h"
#include "qtcad/QgGedEventBatch.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"
#include "wdb.h"

#include <Inventor/nodes/SoGroup.h>

#include <QApplication>

#include <stdio.h>
#include <string.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static SoBRLEditPreview *
find_preview(BObolViewController *controller, const char *previewId)
{
    if (!controller || !previewId)
	return NULL;

    SoNode *root = controller->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    SoGroup *group = static_cast<SoGroup *>(root);
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(SoBRLEditPreview::getClassTypeId()))
	    continue;
	SoBRLEditPreview *preview = static_cast<SoBRLEditPreview *>(node);
	if (bu_strcmp(preview->previewId.getValue().getString(), previewId) == 0)
	    return preview;
    }

    return NULL;
}

static SoBRLVListShape *
preview_shape(SoBRLEditPreview *preview)
{
    if (!preview || preview->getNumChildren() != 1)
	return NULL;

    SoNode *node = preview->getChild(0);
    if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	return NULL;

    return static_cast<SoBRLVListShape *>(node);
}

static int
check_shape_metadata(SoBRLVListShape *shape,
	const char *path,
	const char *previewId,
	const char *editIntentId,
	const char *editIntentRole,
	uint32_t sourceRevision,
	int segmentCount)
{
    if (!shape)
	return 0;
    if (bu_strcmp(shape->sourcePath.getValue().getString(), path) != 0 ||
	    bu_strcmp(shape->sourceName.getValue().getString(), previewId) != 0 ||
	    bu_strcmp(shape->sourceType.getValue().getString(), "edit-preview") != 0 ||
	    shape->sourceId.getValue() != sourceRevision ||
	    !shape->editEmphasis.getValue() ||
	    bu_strcmp(shape->editIntentId.getValue().getString(), editIntentId) != 0 ||
	    bu_strcmp(shape->editIntentRole.getValue().getString(), editIntentRole) != 0 ||
	    shape->getSegmentCount() != segmentCount)
	return 0;
    return 1;
}

static int
make_edit_preview_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0;
    wdb_close(wdbp);
    return ret;
}

static int
edit_preview_update(QgView *view, const char *name, const char *path,
    const char *intent_id, const char *intent_role,
    const SbVec3f *points, const int *commands, int count,
    uint32_t source_revision, uint32_t inputs_revision)
{
    if (!view || !name || !points || !commands || count <= 0)
	return 0;
    point_t *ged_points = (point_t *)bu_calloc((size_t)count,
	sizeof(point_t), "edit transaction test points");
    for (int i = 0; i < count; i++) {
	ged_points[i][X] = points[i][0];
	ged_points[i][Y] = points[i][1];
	ged_points[i][Z] = points[i][2];
    }
    struct ged_view_edit_transaction transaction =
	GED_VIEW_EDIT_TRANSACTION_INIT;
    transaction.event = GED_VIEW_EDIT_PREVIEW_UPDATE;
    transaction.feature_name = name;
    transaction.owner = view;
    transaction.source_path = path;
    transaction.edit_intent_id = intent_id ? intent_id : name;
    transaction.edit_intent_role = intent_role ? intent_role : "preview";
    transaction.points = (const point_t *)ged_points;
    transaction.commands = commands;
    transaction.point_count = (size_t)count;
    transaction.source_revision = source_revision;
    transaction.inputs_revision = inputs_revision;
    const int ret = ged_view_feature_edit_transaction_apply(
	ged_view_context_from_bv(view->viewContext()), &transaction, NULL);
    bu_free(ged_points, "edit transaction test points");
    if (ret)
	view->need_update(QG_VIEW_DRAWN);
    return ret;
}

static int
edit_preview_clear(QgView *view, const char *name)
{
    struct ged_view_edit_transaction transaction =
	GED_VIEW_EDIT_TRANSACTION_INIT;
    transaction.event = GED_VIEW_EDIT_PREVIEW_CANCEL;
    transaction.feature_name = name;
    const int ret = ged_view_feature_edit_transaction_apply(
	ged_view_context_from_bv(view->viewContext()), &transaction, NULL);
    if (ret)
	view->need_update(QG_VIEW_DRAWN);
    return ret;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    QApplication app(argc, argv);

    QgView view(NULL, QgViewType::SW);
    view.resize(160, 120);

    QgPluginContext context;
    context.viewWidgetAccessor = [&view]() -> QgView * { return &view; };
    if (context.getViewWidget() != &view)
	FAIL("plugin context should expose the active QgView widget");
    if (context.activeViewContext() != view.viewContext())
	FAIL("plugin context should expose the active view context");

    BObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");

    const char *dbpath = "qtcad_obol_edit_preview_tmp.g";
    if (!make_edit_preview_db(dbpath))
	FAIL("failed to create qtcad Obol edit-preview test database");
    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol edit-preview test database");
    struct ged_view_context *view_ctx =
	ged_view_context_from_bv(view.viewContext());
    ged_view_active_ctx_set(gedp, view_ctx);
    (void)ged_view_context_host_attach(gedp, view_ctx);
    if (!ged_view_context_obol_endpoint_set(view_ctx,
	    view.displayEndpoint(), 0))
	FAIL("qtcad test should attach the QgView display endpoint to GED");
    SoNode *attachedRenderRoot = controller->getRenderSceneRoot();
    if (!ged_view_context_obol_endpoint_set(view_ctx,
	    view.displayEndpoint(), 0) ||
	controller->getRenderSceneRoot() != attachedRenderRoot)
	FAIL("reaffirming an Obol endpoint should preserve its retained render root");

    QgView secondView(NULL, QgViewType::SW);
    secondView.resize(160, 120);
    struct ged_view_context *second_view_ctx =
	ged_view_context_from_bv(secondView.viewContext());
    BObolViewController *second_controller =
	secondView.obolViewController();
    if (!second_controller ||
	!ged_view_context_host_attach(gedp, second_view_ctx) ||
	!ged_view_set_context_add(ged_view_set_ctx(gedp), second_view_ctx) ||
	!ged_view_context_obol_endpoint_set(second_view_ctx,
	    secondView.displayEndpoint(), 0))
	FAIL("qtcad test should attach a second hosted Obol edit view");
    ged_view_active_ctx_set(gedp, view_ctx);

    const char *previewId = "_test_edit_preview";
    const char *identity = "/box.s::edit-preview";
    SbVec3f points[4] = {
	SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 1.0f, 0.0f),
	SbVec3f(0.0f, 1.0f, 0.0f)
    };
    int commands[4] = {
	GED_DRAW_VIEW_LINE_MOVE,
	GED_DRAW_VIEW_LINE_DRAW,
	GED_DRAW_VIEW_LINE_DRAW,
	GED_DRAW_VIEW_LINE_DRAW
    };

    point_t multi_points[2] = {
	{0.0, 0.0, 0.0},
	{1.0, 1.0, 0.0}
    };
    int multi_commands[2] = {
	GED_DRAW_VIEW_LINE_MOVE,
	GED_DRAW_VIEW_LINE_DRAW
    };
    struct ged_view_edit_transaction multi_transaction =
	GED_VIEW_EDIT_TRANSACTION_INIT;
    multi_transaction.event = GED_VIEW_EDIT_PREVIEW_BEGIN;
    multi_transaction.feature_name = "_test_multi_edit_preview";
    multi_transaction.owner = &view;
    multi_transaction.source_path = "/box.s::multi-edit";
    multi_transaction.edit_intent_id = "edit::box.s/multi";
    multi_transaction.edit_intent_role = "translate";
    multi_transaction.points = (const point_t *)multi_points;
    multi_transaction.commands = multi_commands;
    multi_transaction.point_count = 2;
    if (ged_draw_edit_transaction_apply(gedp, &multi_transaction) < 2 ||
	!controller->features().exists("_test_multi_edit_preview") ||
	!second_controller->features().exists("_test_multi_edit_preview"))
	FAIL("neutral edit transactions should publish to every hosted view");
    multi_transaction.event = GED_VIEW_EDIT_PREVIEW_CANCEL;
    multi_transaction.points = NULL;
    multi_transaction.commands = NULL;
    multi_transaction.point_count = 0;
    if (ged_draw_edit_transaction_apply(gedp, &multi_transaction) < 2 ||
	controller->features().exists("_test_multi_edit_preview") ||
	second_controller->features().exists("_test_multi_edit_preview"))
	FAIL("neutral edit cancellation should retire previews from every hosted view");

    struct directory *multi_box_dp = db_lookup(gedp->dbip, "box.s",
	LOOKUP_QUIET);
    struct rt_db_internal multi_box_intern = RT_DB_INTERNAL_INIT_ZERO;
    if (!multi_box_dp || rt_db_get_internal(&multi_box_intern, multi_box_dp,
	    gedp->dbip, NULL) < 0)
	FAIL("multi-view primitive edit test should load its database object");
    mat_t multi_matrix;
    MAT_IDN(multi_matrix);
    MAT_DELTAS(multi_matrix, 3.0, 0.0, 0.0);
    struct ged_view_edit_transaction primitive_transaction =
	GED_VIEW_EDIT_TRANSACTION_INIT;
    primitive_transaction.event = GED_VIEW_EDIT_PREVIEW_REPLACE_SOURCE;
    primitive_transaction.feature_name = "_test_multi_primitive_edit";
    primitive_transaction.owner = &view;
    primitive_transaction.source_path = "box.s";
    primitive_transaction.edit_intent_id = "edit::box.s/primitive";
    primitive_transaction.edit_intent_role = "primitive-edit";
    primitive_transaction.dbip = gedp->dbip;
    primitive_transaction.internal = &multi_box_intern;
    primitive_transaction.matrix = multi_matrix;
    if (ged_draw_edit_transaction_apply(gedp, &primitive_transaction) < 2 ||
	!controller->features().exists("_test_multi_primitive_edit") ||
	!second_controller->features().exists("_test_multi_primitive_edit")) {
	rt_db_free_internal(&multi_box_intern);
	FAIL("neutral primitive edits should publish transformed feedback to every hosted view");
    }
    primitive_transaction.event = GED_VIEW_EDIT_PREVIEW_COMMIT;
    primitive_transaction.internal = NULL;
    primitive_transaction.matrix = NULL;
    if (ged_draw_edit_transaction_apply(gedp, &primitive_transaction) < 2 ||
	controller->features().exists("_test_multi_primitive_edit") ||
	second_controller->features().exists("_test_multi_primitive_edit")) {
	rt_db_free_internal(&multi_box_intern);
	FAIL("neutral primitive commit should retire feedback from every hosted view");
    }

    const char *draw_box_av[2] = {"draw", "box.s"};
    if (ged_exec_draw(gedp, 2, draw_box_av) != BRLCAD_OK) {
	rt_db_free_internal(&multi_box_intern);
	FAIL("primitive revision test should draw its edited object");
    }
    const uint64_t revision_before_edit = ged_scene_revision(gedp);
    {
	QgGedEventBatch event_batch(gedp);
	if (rt_db_put_internal(multi_box_dp, gedp->dbip,
		&multi_box_intern) < 0) {
	    rt_db_free_internal(&multi_box_intern);
	    FAIL("accepted primitive edit should write through a database event batch");
	}
	(void)ged_event_notify_object_modified(gedp, "box.s", 1, NULL);
    }
    if (ged_scene_revision(gedp) <= revision_before_edit)
	FAIL("accepted primitive edit should advance the GED scene revision");

    controller->clearRenderRequest();
    if (edit_preview_update(&view, previewId, identity, NULL, NULL, points,
	    commands, 4, 7, 11) != 1)
	FAIL("qtcad helper should publish edit preview geometry into Obol");
    if (!controller->isRenderRequested())
	FAIL("edit preview updates should request an Obol render");

    SoBRLEditPreview *preview = find_preview(controller, previewId);
    if (!preview)
	FAIL("Obol scene should contain the published edit preview node");
    if (preview->previewStatus.getValue() != SoBRLEditPreview::CURRENT ||
	    preview->sourceRevision.getValue() != 7 ||
	    preview->inputsRevision.getValue() != 11 ||
	    preview->realizedSourceRevision.getValue() != 7 ||
	    preview->realizedInputsRevision.getValue() != 11)
	FAIL("edit preview node should record explicit source/input revisions");

    SoBRLVListShape *shape = preview_shape(preview);
    if (!check_shape_metadata(shape, identity, previewId, previewId,
	    "preview", 7, 3))
	FAIL("edit preview shape should carry selection/export/edit-intent metadata");

    /*
     * The primitive-revision setup above deliberately leaves box.s drawn.
     * Exercise the preview adapter itself here: applying these query actions
     * to the complete scene makes the assertion depend on asynchronous draw
     * publication order and may count or select the unrelated source.
     */
    SoBRLExportAction previewExport;
    previewExport.apply(preview);
    if (previewExport.getLineCount() != 3 ||
	    bu_strcmp(previewExport.getLine(0).path.getString(), identity) != 0 ||
	    !previewExport.getLine(0).editEmphasis ||
	    bu_strcmp(previewExport.getLine(0).editIntentId.getString(), previewId) != 0 ||
	    bu_strcmp(previewExport.getLine(0).editIntentRole.getString(), "preview") != 0)
	FAIL("qtcad edit preview export should preserve edit-intent metadata");

    SoBRLSnapAction previewSnap;
    previewSnap.setEnabledKinds(SoBRLSnapAction::LINE_NEAREST);
    previewSnap.setQueryPoint(SbVec3f(1.0f, 0.25f, 0.0f));
    previewSnap.setTolerance(0.01f);
    previewSnap.apply(preview);
    if (!previewSnap.hasCandidate() ||
	    bu_strcmp(previewSnap.getPath().getString(), identity) != 0 ||
	    bu_strcmp(previewSnap.getEditIntentId().getString(), previewId) != 0 ||
	    bu_strcmp(previewSnap.getEditIntentRole().getString(), "preview") != 0)
	FAIL("qtcad edit preview snap should preserve edit-intent metadata");

    SoBRLMeasureAction previewMeasure;
    previewMeasure.setQueryPoint(SbVec3f(1.0f, 0.25f, 0.0f));
    previewMeasure.apply(preview);
    if (!previewMeasure.hasNearestSegment() ||
	    bu_strcmp(previewMeasure.getNearestPath().getString(), identity) != 0 ||
	    bu_strcmp(previewMeasure.getNearestEditIntentId().getString(), previewId) != 0 ||
	    bu_strcmp(previewMeasure.getNearestEditIntentRole().getString(), "preview") != 0)
	FAIL("qtcad edit preview measure should preserve edit-intent metadata");

    SbVec3f replacementPoints[2] = {
	SbVec3f(2.0f, 0.0f, 0.0f),
	SbVec3f(2.0f, 1.0f, 0.0f)
    };
    int replacementCommands[2] = {
	GED_DRAW_VIEW_LINE_MOVE,
	GED_DRAW_VIEW_LINE_DRAW
    };

    controller->clearRenderRequest();
    if (edit_preview_update(&view, previewId, identity, NULL, NULL,
	    replacementPoints, replacementCommands, 2, 0, 0) != 1)
	FAIL("qtcad helper should replace existing edit preview geometry");

    preview = find_preview(controller, previewId);
    if (!preview)
	FAIL("replacement should keep the edit preview node present");
    if (preview->getNumChildren() != 1 ||
	    preview->sourceRevision.getValue() != 8 ||
	    preview->inputsRevision.getValue() != 12 ||
	    preview->realizedSourceRevision.getValue() != 8 ||
	    preview->realizedInputsRevision.getValue() != 12)
	FAIL("implicit edit preview revisions should advance from the previous values");

    shape = preview_shape(preview);
    if (!check_shape_metadata(shape, identity, previewId, previewId,
	    "preview", 8, 1))
	FAIL("replacement should drop stale geometry and update shape metadata");

    SoBRLExportAction replacementExport;
    replacementExport.apply(preview);
    if (replacementExport.getLineCount() != 1 ||
	    bu_strcmp(replacementExport.getLine(0).editIntentId.getString(), previewId) != 0 ||
	    bu_strcmp(replacementExport.getLine(0).editIntentRole.getString(), "preview") != 0)
	FAIL("replacement edit preview export should keep edit-intent metadata current");

    const char *customIntentId = "edit::box.s/translate";
    const char *customIntentRole = "move-handle";
    controller->clearRenderRequest();
    if (edit_preview_update(&view, previewId, identity,
	    customIntentId, customIntentRole, replacementPoints,
	    replacementCommands, 2, 21, 22) != 1)
	FAIL("qtcad helper should publish explicit edit intent metadata");
    if (!controller->isRenderRequested())
	FAIL("custom edit preview updates should request an Obol render");

    preview = find_preview(controller, previewId);
    if (!preview)
	FAIL("custom intent replacement should keep the edit preview node present");
    if (bu_strcmp(preview->editIntentId.getValue().getString(), customIntentId) != 0 ||
	    bu_strcmp(preview->editIntentRole.getValue().getString(), customIntentRole) != 0 ||
	    preview->sourceRevision.getValue() != 21 ||
	    preview->inputsRevision.getValue() != 22)
	FAIL("edit preview node should record explicit live edit intent fields");

    shape = preview_shape(preview);
    if (!check_shape_metadata(shape, identity, previewId, customIntentId,
	    customIntentRole, 21, 1))
	FAIL("custom intent replacement should update shape edit-intent metadata");

    SoBRLExportAction customIntentExport;
    customIntentExport.apply(preview);
    if (customIntentExport.getLineCount() != 1 ||
	    bu_strcmp(customIntentExport.getLine(0).editIntentId.getString(),
		customIntentId) != 0 ||
	    bu_strcmp(customIntentExport.getLine(0).editIntentRole.getString(),
		customIntentRole) != 0)
	FAIL("custom edit preview export should preserve explicit intent metadata");

    SoBRLSnapAction customIntentSnap;
    customIntentSnap.setEnabledKinds(SoBRLSnapAction::LINE_NEAREST);
    customIntentSnap.setQueryPoint(SbVec3f(2.0f, 0.25f, 0.0f));
    customIntentSnap.setTolerance(0.01f);
    customIntentSnap.apply(preview);
    if (!customIntentSnap.hasCandidate() ||
	    bu_strcmp(customIntentSnap.getEditIntentId().getString(),
		customIntentId) != 0 ||
	    bu_strcmp(customIntentSnap.getEditIntentRole().getString(),
		customIntentRole) != 0)
	FAIL("custom edit preview snap should preserve explicit intent metadata");

    SoBRLMeasureAction customIntentMeasure;
    customIntentMeasure.setQueryPoint(SbVec3f(2.0f, 0.25f, 0.0f));
    customIntentMeasure.apply(preview);
    if (!customIntentMeasure.hasNearestSegment() ||
	    bu_strcmp(customIntentMeasure.getNearestEditIntentId().getString(),
		customIntentId) != 0 ||
	    bu_strcmp(customIntentMeasure.getNearestEditIntentRole().getString(),
		customIntentRole) != 0)
	FAIL("custom edit preview measure should preserve explicit intent metadata");

    controller->clearRenderRequest();
    if (edit_preview_clear(&view, previewId) != 1)
	FAIL("qtcad helper should clear edit preview geometry");
    if (find_preview(controller, previewId))
	FAIL("Obol scene should remove cleared edit preview nodes");
    if (!controller->isRenderRequested())
	FAIL("edit preview clears should request an Obol render");

    const char *gedPreviewId = "_test_ged_edit_preview";
    const char *gedIdentity = "/box.s::ged-edit-preview";
    const char *gedIntentId = "edit::box.s/ged-preview";
    const char *gedIntentRole = "scale-handle";
    controller->clearRenderRequest();
    (void)bv_context_refresh_complete(static_cast<struct bv_context *>(view.viewContext()));
    if (edit_preview_update(&view, gedPreviewId,
	    gedIdentity, gedIntentId, gedIntentRole, points, commands, 4,
	    41, 42) != 1)
	FAIL("qtcad helper should publish GED-owned edit preview geometry");
    if (!bv_refresh_dirty_get(bv_context_view_const(static_cast<const struct bv_context *>(view.viewContext()))))
	FAIL("GED-routed edit preview updates should request a qtcad view refresh");

    struct ged_view_feature_summary gedSummary =
	GED_VIEW_FEATURE_SUMMARY_INIT;
    if (!ged_view_feature_get_summary(view_ctx, gedPreviewId,
	    &gedSummary) ||
	    !gedSummary.exists ||
	    gedSummary.kind != GED_VIEW_FEATURE_KIND_EDIT_PREVIEW ||
	    gedSummary.scope != GED_VIEW_FEATURE_SCOPE_LOCAL ||
	    gedSummary.overlay_class !=
	    GED_VIEW_FEATURE_OVERLAY_CLASS_EDIT_HANDLE ||
	    gedSummary.lifecycle != GED_VIEW_FEATURE_LIFECYCLE_PER_TOOL ||
	    !gedSummary.is_transient_preview ||
	    !gedSummary.owner_id[0])
	FAIL("GED-routed qtcad edit preview should be a local transient feature with owner metadata");

    BObolFeatureHandle gedHandle =
	controller->features().find(gedPreviewId);
    BObolFeatureRecord gedRecord;
    if (!gedHandle.isValid() ||
	    !controller->features().record(gedHandle, gedRecord) ||
	    gedRecord.kind != BObolFeatureKind::EditPreview ||
	    gedRecord.identity != gedIdentity ||
	    gedRecord.editIntentId != gedIntentId ||
	    gedRecord.editIntentRole != gedIntentRole ||
	    gedRecord.sourceRevision != 41 ||
	    gedRecord.inputsRevision != 42 ||
	    gedRecord.scope != BObolFeatureScope::Local)
	FAIL("GED-routed qtcad edit preview should preserve identity, intent, revisions, and local scope");

    controller->clearRenderRequest();
    (void)bv_context_refresh_complete(static_cast<struct bv_context *>(view.viewContext()));
    if (edit_preview_clear(&view, gedPreviewId) != 1)
	FAIL("qtcad helper should clear GED-routed edit preview geometry");
    if (ged_view_feature_get_summary(view_ctx, gedPreviewId,
	    &gedSummary) && gedSummary.exists)
	FAIL("GED-routed edit preview clear should remove the transient feature");
    if (!bv_refresh_dirty_get(bv_context_view_const(static_cast<const struct bv_context *>(view.viewContext()))))
	FAIL("GED-routed edit preview clears should request a qtcad view refresh");

    (void)ged_view_context_obol_endpoint_set(second_view_ctx, NULL, 0);
    (void)ged_view_context_obol_endpoint_set(view_ctx, NULL, 0);
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
