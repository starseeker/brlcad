/*        T E S T _ Q T C A D _ O B O L _ E D I T _ P R E V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/edit_preview.h"
#include "brlobol/view_controller.h"
#include "brlobol/vlist_shape.h"
#include "qtcad/QgObolEditPreview.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgView.h"

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
	    shape->getSegmentCount() != segmentCount)
	return 0;
    return 1;
}

int
main(int argc, char **argv)
{
    QApplication app(argc, argv);

    QgView view(NULL, QgView_SW, NULL);
    view.resize(160, 120);

    QgPluginContext context;
    context.viewWidgetAccessor = [&view]() -> QgView * { return &view; };
    if (context.getViewWidget() != &view)
	FAIL("plugin context should expose the active QgView widget");

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
    if (!check_shape_metadata(shape, identity, previewId, 7, 3))
	FAIL("edit preview shape should carry selection/export metadata");

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
    if (!check_shape_metadata(shape, identity, previewId, 8, 1))
	FAIL("replacement should drop stale geometry and update shape metadata");

    controller->clearRenderRequest();
    if (qg_obol_edit_preview_clear(&view, previewId) != 1)
	FAIL("qtcad helper should clear edit preview geometry");
    if (find_preview(controller, previewId))
	FAIL("Obol scene should remove cleared edit preview nodes");
    if (!controller->isRenderRequested() ||
	    strcmp(controller->getRenderReason().getString(), "edit-preview") != 0)
	FAIL("edit preview clears should request an Obol render");

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
