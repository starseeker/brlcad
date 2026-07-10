/*        T E S T _ Q T C A D _ O B O L _ E D I T _ P R E V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/edit_preview.h"
#include "brlobol/export_action.h"
#include "brlobol/measure_action.h"
#include "brlobol/snap_action.h"
#include "brlobol/view_controller.h"
#include "brlobol/view_store.h"
#include "brlobol/vlist_shape.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "QgLegacyViewContext.h"
#include "QgObolEditPreviewPrivate.h"
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgPluginContext.h"
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
find_preview(BRLObolViewController *controller, const char *previewId)
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
	if (strcmp(preview->previewId.getValue().getString(), previewId) == 0)
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
    if (strcmp(shape->sourcePath.getValue().getString(), path) != 0 ||
	    strcmp(shape->sourceName.getValue().getString(), previewId) != 0 ||
	    strcmp(shape->sourceType.getValue().getString(), "edit-preview") != 0 ||
	    shape->sourceId.getValue() != sourceRevision ||
	    !shape->editEmphasis.getValue() ||
	    strcmp(shape->editIntentId.getValue().getString(), editIntentId) != 0 ||
	    strcmp(shape->editIntentRole.getValue().getString(), editIntentRole) != 0 ||
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

int
main(int argc, char **argv)
{
    QApplication app(argc, argv);

    QgView view(NULL, QgView_SW);
    view.resize(160, 120);

    QgPluginContext context;
    context.viewWidgetAccessor = [&view]() -> QgView * { return &view; };
    if (context.getViewWidget() != &view)
	FAIL("plugin context should expose the active QgView widget");
    if (context.activeView() != view.view())
	FAIL("plugin context should expose the active opaque legacy view");

    BRLObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");

    const char *previewId = "_test_edit_preview";
    const char *identity = "/box.s::edit-preview";
    SbVec3f points[4] = {
	SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 1.0f, 0.0f),
	SbVec3f(0.0f, 1.0f, 0.0f)
    };
    int32_t commands[4] = {
	QG_OBOL_EDIT_PREVIEW_MOVE,
	QG_OBOL_EDIT_PREVIEW_DRAW,
	QG_OBOL_EDIT_PREVIEW_DRAW,
	QG_OBOL_EDIT_PREVIEW_DRAW
    };

    controller->clearRenderRequest();
    if (qg_obol_edit_preview_update(&view, previewId, identity, points,
	    commands, 4, 7, 11) != 1)
	FAIL("qtcad helper should publish edit preview geometry into Obol");
    if (!controller->isRenderRequested() ||
	    strcmp(controller->getRenderReason().getString(), "edit-preview") != 0)
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

    SoBRLExportAction previewExport;
    previewExport.apply(controller->getSceneRoot());
    if (previewExport.getLineCount() != 3 ||
	    strcmp(previewExport.getLine(0).path.getString(), identity) != 0 ||
	    !previewExport.getLine(0).editEmphasis ||
	    strcmp(previewExport.getLine(0).editIntentId.getString(), previewId) != 0 ||
	    strcmp(previewExport.getLine(0).editIntentRole.getString(), "preview") != 0)
	FAIL("qtcad edit preview export should preserve edit-intent metadata");

    SoBRLSnapAction previewSnap;
    previewSnap.setEnabledKinds(SoBRLSnapAction::LINE_NEAREST);
    previewSnap.setQueryPoint(SbVec3f(1.0f, 0.25f, 0.0f));
    previewSnap.setTolerance(0.01f);
    previewSnap.apply(controller->getSceneRoot());
    if (!previewSnap.hasCandidate() ||
	    strcmp(previewSnap.getPath().getString(), identity) != 0 ||
	    strcmp(previewSnap.getEditIntentId().getString(), previewId) != 0 ||
	    strcmp(previewSnap.getEditIntentRole().getString(), "preview") != 0)
	FAIL("qtcad edit preview snap should preserve edit-intent metadata");

    SoBRLMeasureAction previewMeasure;
    previewMeasure.setQueryPoint(SbVec3f(1.0f, 0.25f, 0.0f));
    previewMeasure.apply(controller->getSceneRoot());
    if (!previewMeasure.hasNearestSegment() ||
	    strcmp(previewMeasure.getNearestPath().getString(), identity) != 0 ||
	    strcmp(previewMeasure.getNearestEditIntentId().getString(), previewId) != 0 ||
	    strcmp(previewMeasure.getNearestEditIntentRole().getString(), "preview") != 0)
	FAIL("qtcad edit preview measure should preserve edit-intent metadata");

    SbVec3f replacementPoints[2] = {
	SbVec3f(2.0f, 0.0f, 0.0f),
	SbVec3f(2.0f, 1.0f, 0.0f)
    };
    int32_t replacementCommands[2] = {
	QG_OBOL_EDIT_PREVIEW_MOVE,
	QG_OBOL_EDIT_PREVIEW_DRAW
    };

    controller->clearRenderRequest();
    if (qg_obol_edit_preview_update(&view, previewId, identity,
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
    replacementExport.apply(controller->getSceneRoot());
    if (replacementExport.getLineCount() != 1 ||
	    strcmp(replacementExport.getLine(0).editIntentId.getString(), previewId) != 0 ||
	    strcmp(replacementExport.getLine(0).editIntentRole.getString(), "preview") != 0)
	FAIL("replacement edit preview export should keep edit-intent metadata current");

    const char *customIntentId = "edit::box.s/translate";
    const char *customIntentRole = "move-handle";
    controller->clearRenderRequest();
    if (qg_obol_edit_preview_update_with_intent(&view, previewId, identity,
	    customIntentId, customIntentRole, replacementPoints,
	    replacementCommands, 2, 21, 22) != 1)
	FAIL("qtcad helper should publish explicit edit intent metadata");
    if (!controller->isRenderRequested() ||
	    strcmp(controller->getRenderReason().getString(), "edit-preview") != 0)
	FAIL("custom edit preview updates should request an Obol render");

    preview = find_preview(controller, previewId);
    if (!preview)
	FAIL("custom intent replacement should keep the edit preview node present");
    if (strcmp(preview->editIntentId.getValue().getString(), customIntentId) != 0 ||
	    strcmp(preview->editIntentRole.getValue().getString(), customIntentRole) != 0 ||
	    preview->sourceRevision.getValue() != 21 ||
	    preview->inputsRevision.getValue() != 22)
	FAIL("edit preview node should record explicit live edit intent fields");

    shape = preview_shape(preview);
    if (!check_shape_metadata(shape, identity, previewId, customIntentId,
	    customIntentRole, 21, 1))
	FAIL("custom intent replacement should update shape edit-intent metadata");

    SoBRLExportAction customIntentExport;
    customIntentExport.apply(controller->getSceneRoot());
    if (customIntentExport.getLineCount() != 1 ||
	    strcmp(customIntentExport.getLine(0).editIntentId.getString(),
		customIntentId) != 0 ||
	    strcmp(customIntentExport.getLine(0).editIntentRole.getString(),
		customIntentRole) != 0)
	FAIL("custom edit preview export should preserve explicit intent metadata");

    SoBRLSnapAction customIntentSnap;
    customIntentSnap.setEnabledKinds(SoBRLSnapAction::LINE_NEAREST);
    customIntentSnap.setQueryPoint(SbVec3f(2.0f, 0.25f, 0.0f));
    customIntentSnap.setTolerance(0.01f);
    customIntentSnap.apply(controller->getSceneRoot());
    if (!customIntentSnap.hasCandidate() ||
	    strcmp(customIntentSnap.getEditIntentId().getString(),
		customIntentId) != 0 ||
	    strcmp(customIntentSnap.getEditIntentRole().getString(),
		customIntentRole) != 0)
	FAIL("custom edit preview snap should preserve explicit intent metadata");

    SoBRLMeasureAction customIntentMeasure;
    customIntentMeasure.setQueryPoint(SbVec3f(2.0f, 0.25f, 0.0f));
    customIntentMeasure.apply(controller->getSceneRoot());
    if (!customIntentMeasure.hasNearestSegment() ||
	    strcmp(customIntentMeasure.getNearestEditIntentId().getString(),
		customIntentId) != 0 ||
	    strcmp(customIntentMeasure.getNearestEditIntentRole().getString(),
		customIntentRole) != 0)
	FAIL("custom edit preview measure should preserve explicit intent metadata");

    controller->clearRenderRequest();
    if (qg_obol_edit_preview_clear(&view, previewId) != 1)
	FAIL("qtcad helper should clear edit preview geometry");
    if (find_preview(controller, previewId))
	FAIL("Obol scene should remove cleared edit preview nodes");
    if (!controller->isRenderRequested() ||
	    strcmp(controller->getRenderReason().getString(), "edit-preview") != 0)
	FAIL("edit preview clears should request an Obol render");

    const char *dbpath = "qtcad_obol_edit_preview_tmp.g";
    if (!make_edit_preview_db(dbpath))
	FAIL("failed to create qtcad Obol edit-preview test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol edit-preview test database");

    void *view_ctx = qg_legacy_view_to_context(view.view());
    qg_legacy_view_ged_active_set(gedp, view.view());
    if (!ged_draw_obol_controller_attach_for_view(gedp, view_ctx,
	    controller, 0))
	FAIL("qtcad test should attach the QgView Obol controller to GED");

    const char *gedPreviewId = "_test_ged_edit_preview";
    const char *gedIdentity = "/box.s::ged-edit-preview";
    const char *gedIntentId = "edit::box.s/ged-preview";
    const char *gedIntentRole = "scale-handle";
    controller->clearRenderRequest();
    (void)bv_context_refresh_complete(qg_legacy_view_context(view.view()));
    if (qg_obol_edit_preview_update_with_intent(&view, gedPreviewId,
	    gedIdentity, gedIntentId, gedIntentRole, points, commands, 4,
	    41, 42) != 1)
	FAIL("qtcad helper should publish GED-owned edit preview geometry");
    if (!bv_refresh_dirty_get(qg_legacy_view_bv_const(view.view())))
	FAIL("GED-routed edit preview updates should request a qtcad view refresh");

    struct ged_draw_view_feature_summary gedSummary =
	GED_DRAW_VIEW_FEATURE_SUMMARY_INIT;
    if (!ged_draw_view_context_feature_summary(view_ctx, gedPreviewId,
	    &gedSummary) ||
	    !gedSummary.exists ||
	    gedSummary.kind != GED_DRAW_VIEW_FEATURE_KIND_EDIT_PREVIEW ||
	    gedSummary.scope != GED_DRAW_VIEW_FEATURE_SCOPE_LOCAL ||
	    gedSummary.overlay_class !=
	    GED_DRAW_VIEW_FEATURE_OVERLAY_CLASS_EDIT_HANDLE ||
	    gedSummary.lifecycle != GED_DRAW_VIEW_FEATURE_LIFECYCLE_PER_TOOL ||
	    !gedSummary.is_transient_preview ||
	    !gedSummary.owner_id[0])
	FAIL("GED-routed qtcad edit preview should be a local transient feature with owner metadata");

    BRLObolFeatureHandle gedHandle =
	controller->features().find(gedPreviewId);
    BRLObolFeatureRecord gedRecord;
    if (!gedHandle.isValid() ||
	    !controller->features().record(gedHandle, gedRecord) ||
	    gedRecord.kind != BRLObolFeatureKind::EditPreview ||
	    gedRecord.identity != gedIdentity ||
	    gedRecord.editIntentId != gedIntentId ||
	    gedRecord.editIntentRole != gedIntentRole ||
	    gedRecord.sourceRevision != 41 ||
	    gedRecord.inputsRevision != 42 ||
	    gedRecord.scope != BRLObolFeatureScope::Local)
	FAIL("GED-routed qtcad edit preview should preserve identity, intent, revisions, and local scope");

    controller->clearRenderRequest();
    (void)bv_context_refresh_complete(qg_legacy_view_context(view.view()));
    if (qg_obol_edit_preview_clear(&view, gedPreviewId) != 1)
	FAIL("qtcad helper should clear GED-routed edit preview geometry");
    if (ged_draw_view_context_feature_summary(view_ctx, gedPreviewId,
	    &gedSummary) && gedSummary.exists)
	FAIL("GED-routed edit preview clear should remove the transient feature");
    if (!bv_refresh_dirty_get(qg_legacy_view_bv_const(view.view())))
	FAIL("GED-routed edit preview clears should request a qtcad view refresh");

    ged_draw_obol_controller_detach_for_view(gedp, view_ctx);
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
