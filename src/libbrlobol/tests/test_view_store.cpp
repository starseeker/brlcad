/*                    T E S T _ V I E W _ S T O R E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol.h"
#include "brlobol/edit_preview.h"
#include "brlobol/hud_label_overlay.h"
#include "bu/app.h"
#include "bu/file.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/SoType.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoFont.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoText2.h>

#include <stdio.h>
#include <string.h>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static int
count_nodes_of_type(SoNode *node, SoType type)
{
    if (!node)
	return 0;

    int count = node->isOfType(type) ? 1 : 0;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    count += count_nodes_of_type(group->getChild(i), type);
    }
    return count;
}

static SoNode *
first_node_of_type(SoNode *node, SoType type)
{
    if (!node)
	return NULL;
    if (node->isOfType(type))
	return node;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SoNode *found = first_node_of_type(group->getChild(i), type);
	    if (found)
		return found;
	}
    }
    return NULL;
}

struct feature_visit_count {
    size_t count;
};

struct selection_visit_state {
    size_t count;
    const char *last;
};

struct edit_preview_callback_state {
    uint64_t revision;
    int updateCount;
};

static uint64_t
edit_preview_revision_cb(void *data)
{
    struct edit_preview_callback_state *ctx =
	(struct edit_preview_callback_state *)data;
    return ctx ? ctx->revision : 0;
}

static int
edit_preview_update_cb(void *data)
{
    struct edit_preview_callback_state *ctx =
	(struct edit_preview_callback_state *)data;
    if (!ctx)
	return 0;
    ctx->updateCount++;
    return 1;
}

static int
count_feature_visit_cb(const BRLObolFeatureRecord &record, void *userData)
{
    struct feature_visit_count *ctx = (struct feature_visit_count *)userData;
    if (!ctx || record.name.getLength() == 0)
	return 0;
    ctx->count++;
    return 1;
}

static int
count_selection_path_cb(const SbString &path, void *userData)
{
    struct selection_visit_state *ctx =
	(struct selection_visit_state *)userData;
    if (!ctx || path.getLength() == 0)
	return 0;
    ctx->count++;
    ctx->last = path.getString();
    return 1;
}

static int
test_feature_nodes(BRLObolViewController &view)
{
    std::vector<BRLObolLabel> labels;
    BRLObolLabel label;
    label.text = "store label";
    label.point = SbVec3f(1.0f, 2.0f, 0.0f);
    label.hasLeader = TRUE;
    label.target = SbVec3f(0.0f, 0.0f, 0.0f);
    label.fontSize = 17.0f;
    labels.push_back(label);

    BRLObolFeatureHandle labelHandle = view.features().publishLabels(
					   "labels",
					   BRLObolFeatureScope::Shared,
					   labels);
    if (!labelHandle.isValid())
	FAIL("publishLabels should return a valid handle");
    if (count_nodes_of_type(view.features().node(labelHandle),
			    SoText2::getClassTypeId()) != 1)
	FAIL("label feature should realize a SoText2 node");
    SoFont *font = static_cast<SoFont *>(first_node_of_type(
	    view.features().node(labelHandle), SoFont::getClassTypeId()));
    if (!font || font->size.getValue() != 17.0f)
	FAIL("label feature should preserve explicit font size");

    std::vector<BRLObolLabel> hudLabels;
    BRLObolLabel hudLabel;
    hudLabel.text = "store hud label";
    hudLabel.point = SbVec3f(8.0f, 58.0f, 0.0f);
    hudLabel.hasColor = TRUE;
    hudLabel.color = SbColor(1.0f, 1.0f, 0.0f);
    hudLabel.fontSize = 12.0f;
    hudLabel.sourceId = 42;
    hudLabels.push_back(hudLabel);

    BRLObolFeatureHandle hudHandle = view.features().publishHudLabels(
					 "hud-labels",
					 BRLObolFeatureScope::Shared,
					 hudLabels);
    if (!hudHandle.isValid())
	FAIL("publishHudLabels should return a valid handle");
    if (count_nodes_of_type(view.features().node(hudHandle),
			    SoBRLHUDLabelOverlay::getClassTypeId()) != 1)
	FAIL("HUD label feature should realize a SoBRLHUDLabelOverlay node");
    SoBRLHUDLabelOverlay *hud = static_cast<SoBRLHUDLabelOverlay *>(
				    first_node_of_type(view.features().node(hudHandle),
					SoBRLHUDLabelOverlay::getClassTypeId()));
    if (!hud || !hud->getHUDLabel() || hud->sourceId.getValue() != 42)
	FAIL("HUD label feature should preserve source and build HUD geometry");
    if (hud->position.getValue()[0] != 8.0f ||
	hud->position.getValue()[1] != 58.0f ||
	hud->fontSize.getValue() != 12.0f)
	FAIL("HUD label feature should preserve screen position and font size");

    std::vector<SbVec3f> points;
    points.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    points.push_back(SbVec3f(1.0f, 0.0f, 0.0f));
    points.push_back(SbVec3f(1.0f, 1.0f, 0.0f));
    points.push_back(SbVec3f(0.0f, 1.0f, 0.0f));
    std::vector<SbVec3f> normals;
    std::vector<int32_t> indices;
    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(2);
    indices.push_back(0);
    indices.push_back(2);
    indices.push_back(3);

    BRLObolFeatureHandle meshHandle = view.features().publishIndexedFaceSet(
					  "surface",
					  BRLObolFeatureScope::Shared,
					  points,
					  normals,
					  indices);
    if (!meshHandle.isValid())
	FAIL("publishIndexedFaceSet should return a valid handle");

    SoNode *node = view.features().node(meshHandle);
    if (!node || !node->isOfType(SoBRLMeshShape::getClassTypeId()))
	FAIL("indexed face feature should realize a SoBRLMeshShape");
    SoBRLMeshShape *mesh = static_cast<SoBRLMeshShape *>(node);
    if (mesh->getTriangleCount() != 2)
	FAIL("indexed face feature should triangulate to two triangles");

    BRLObolFeatureOwner ownerA;
    ownerA.ownerToken = (const void *)0x1;
    ownerA.ownerId = "view-A";
    ownerA.ownerRole = "view";
    BRLObolFeatureOwner ownerB;
    ownerB.ownerToken = (const void *)0x2;
    ownerB.ownerId = "view-B";
    ownerB.ownerRole = "view";

    std::vector<int32_t> lineCommands;
    lineCommands.push_back(static_cast<int32_t>(BRLObolLineCommand::Move));
    lineCommands.push_back(static_cast<int32_t>(BRLObolLineCommand::Draw));
    BRLObolFeatureHandle localA = view.features().publishLineSet(
				      "local-line",
				      BRLObolFeatureScope::Local,
				      points,
				      lineCommands,
				      NULL,
				      &ownerA);
    BRLObolFeatureHandle localB = view.features().publishLineSet(
				      "local-line",
				      BRLObolFeatureScope::Local,
				      points,
				      lineCommands,
				      NULL,
				      &ownerB);
    if (!localA.isValid() || !localB.isValid() || localA.id == localB.id)
	FAIL("owner-scoped local feature names should not collide");
    if (!view.features().existsOwned("local-line",
				     BRLOBOL_FEATURE_SCOPE_LOCAL, &ownerA))
	FAIL("owner A should find its local feature");
    if (!view.features().existsOwned("local-line",
				     BRLOBOL_FEATURE_SCOPE_LOCAL, &ownerB))
	FAIL("owner B should find its local feature");

    struct feature_visit_count visitA = {0};
    view.features().visitRecords(count_feature_visit_cb, &visitA,
				 BRLOBOL_FEATURE_SCOPE_LOCAL, &ownerA);
    if (visitA.count != 1)
	FAIL("owner-scoped feature visit should see only one local record");

    BRLObolFeatureRecord localRecord;
    if (!view.features().record(localA, localRecord) ||
	localRecord.points.size() != points.size() ||
	localRecord.commands.size() != lineCommands.size() ||
	localRecord.owner.ownerToken != ownerA.ownerToken)
	FAIL("feature record snapshot should preserve geometry and owner");

    BRLObolOverlayInfo overlay;
    overlay.isOverlay = TRUE;
    overlay.ownerToken = ownerA.ownerToken;
    overlay.role = BRLObolOverlayRole::Model;
    overlay.overlayClass = BRLObolOverlayClass::TclOverlay;
    overlay.lifecycle = BRLObolOverlayLifecycle::PerCommand;
    overlay.order = BRLObolOverlayOrder::PostTransparent;
    overlay.sortOrder = 7;
    overlay.sourcePath = "local-line";
    if (!view.features().setOverlayInfo(localA, overlay))
	FAIL("feature overlay metadata should be settable");

    BRLObolOverlayInfo overlayRead;
    if (!view.features().overlayInfo(localA, overlayRead) ||
	!overlayRead.isOverlay ||
	overlayRead.ownerToken != ownerA.ownerToken ||
	overlayRead.role != BRLObolOverlayRole::Model ||
	overlayRead.overlayClass != BRLObolOverlayClass::TclOverlay ||
	overlayRead.lifecycle != BRLObolOverlayLifecycle::PerCommand ||
	overlayRead.order != BRLObolOverlayOrder::PostTransparent ||
	overlayRead.sortOrder != 7 ||
	strcmp(overlayRead.sourcePath.getString(), "local-line") != 0)
	FAIL("feature overlay metadata should round-trip through overlayInfo");

    BRLObolFeatureSummary overlaySummary;
    if (!view.features().summaryOwned("local-line", overlaySummary,
				      BRLOBOL_FEATURE_SCOPE_LOCAL, &ownerA) ||
	!overlaySummary.exists ||
	!overlaySummary.overlay.isOverlay ||
	overlaySummary.overlay.overlayClass !=
	BRLObolOverlayClass::TclOverlay)
	FAIL("feature summary should preserve overlay metadata");

    BRLObolFeatureRecord overlayRecord;
    if (!view.features().record(localA, overlayRecord) ||
	!overlayRecord.overlay.isOverlay ||
	overlayRecord.overlay.lifecycle !=
	BRLObolOverlayLifecycle::PerCommand)
	FAIL("feature record should preserve overlay metadata");

    std::vector<BRLObolFeatureMetadata> metadata;
    BRLObolFeatureMetadata metadataItem;
    metadataItem.key = "result.kind";
    metadataItem.value = "diagnostic-line";
    metadata.push_back(metadataItem);
    metadataItem.key = "result.index";
    metadataItem.value = "7";
    metadata.push_back(metadataItem);
    if (!view.features().replaceMetadata(localA, metadata))
	FAIL("feature metadata should be settable");
    std::vector<BRLObolFeatureMetadata> metadataRead;
    if (!view.features().metadata(localA, metadataRead) ||
	metadataRead.size() != 2 ||
	strcmp(metadataRead[0].key.getString(), "result.kind") != 0 ||
	strcmp(metadataRead[0].value.getString(),
	       "diagnostic-line") != 0 ||
	strcmp(metadataRead[1].key.getString(), "result.index") != 0 ||
	strcmp(metadataRead[1].value.getString(), "7") != 0)
	FAIL("feature metadata should round-trip through metadata readback");

    if (!view.features().summaryOwned("local-line", overlaySummary,
				      BRLOBOL_FEATURE_SCOPE_LOCAL, &ownerA) ||
	overlaySummary.metadataCount != 2)
	FAIL("feature summary should report metadata count");

    if (!view.features().record(localA, overlayRecord) ||
	overlayRecord.metadata.size() != 2 ||
	strcmp(overlayRecord.metadata[0].key.getString(),
	       "result.kind") != 0)
	FAIL("feature record should preserve metadata");

    std::vector<BRLObolFeatureMetadata> primitiveMetadata;
    metadataItem.key = "overlap.objects";
    metadataItem.value = "box.s cone.s";
    primitiveMetadata.push_back(metadataItem);
    metadataItem.key = "overlap.depth";
    metadataItem.value = "0.25";
    primitiveMetadata.push_back(metadataItem);
    if (!view.features().replacePrimitiveMetadata(localA, 0,
	    primitiveMetadata))
	FAIL("feature primitive metadata should be settable");
    if (!view.features().primitiveMetadata(localA, 0, metadataRead) ||
	metadataRead.size() != 2 ||
	strcmp(metadataRead[0].key.getString(), "overlap.objects") != 0 ||
	strcmp(metadataRead[0].value.getString(),
	       "box.s cone.s") != 0 ||
	strcmp(metadataRead[1].key.getString(), "overlap.depth") != 0 ||
	strcmp(metadataRead[1].value.getString(), "0.25") != 0)
	FAIL("feature primitive metadata should round-trip through readback");
    if (!view.features().summaryOwned("local-line", overlaySummary,
				      BRLOBOL_FEATURE_SCOPE_LOCAL, &ownerA) ||
	overlaySummary.primitiveMetadataCount != 1)
	FAIL("feature summary should report primitive metadata count");
    if (!view.features().record(localA, overlayRecord) ||
	overlayRecord.primitiveMetadata.size() != 1 ||
	overlayRecord.primitiveMetadata[0].primitiveIndex != 0 ||
	overlayRecord.primitiveMetadata[0].metadata.size() != 2)
	FAIL("feature record should preserve primitive metadata");

    SoNode *localNode = view.features().node(localA);
    if (!localNode ||
	!localNode->isOfType(SoBRLVListShape::getClassTypeId()) ||
	strcmp(static_cast<SoBRLVListShape *>(localNode)->
	       recordRole.getValue().getString(), "view-feature") != 0)
	FAIL("feature overlay metadata should preserve feature-store node role");

    if (!view.features().clearOverlayInfo(localA) ||
	!view.features().overlayInfo(localA, overlayRead) ||
	overlayRead.isOverlay)
	FAIL("feature overlay metadata should be clearable");
    localNode = view.features().node(localA);
    if (!localNode ||
	!localNode->isOfType(SoBRLVListShape::getClassTypeId()) ||
	strcmp(static_cast<SoBRLVListShape *>(localNode)->
	       recordRole.getValue().getString(), "view-feature") != 0)
	FAIL("cleared feature overlay metadata should restore view-feature node role");

    localA = view.features().publishLineSet("local-line",
					    BRLObolFeatureScope::Local, points, lineCommands, NULL, &ownerA);
    if (!localA.isValid() ||
	!view.features().metadata(localA, metadataRead) ||
	!metadataRead.empty() ||
	!view.features().primitiveMetadata(localA, 0, metadataRead) ||
	!metadataRead.empty())
	FAIL("feature publish should clear stale metadata");

    BRLObolFeatureStyle unselectableStyle;
    unselectableStyle.hasSelectable = TRUE;
    unselectableStyle.selectable = FALSE;
    if (!view.features().applyStyle(localA, unselectableStyle))
	FAIL("feature selectable style should be applicable");
    localNode = view.features().node(localA);
    if (!localNode ||
	!localNode->isOfType(SoBRLVListShape::getClassTypeId()) ||
	static_cast<SoBRLVListShape *>(localNode)->selectable.getValue())
	FAIL("feature selectable style should control realized line shapes");

    BRLObolFeaturePrimitivePick primitivePick;
    if (!view.features().resolvePrimitivePick("local-line", 0,
	    primitivePick, BRLOBOL_FEATURE_SCOPE_LOCAL, &ownerA) ||
	primitivePick.handle.id != localA.id ||
	strcmp(primitivePick.featureName.getString(), "local-line") != 0 ||
	primitivePick.primitiveIndex != 0)
	FAIL("feature primitive pick resolver should find direct local features");

    BRLObolLineLayer layer;
    layer.name = "layered-line/red";
    layer.points = points;
    layer.commands = lineCommands;
    std::vector<BRLObolLineLayer> layers;
    layers.push_back(layer);
    BRLObolFeatureHandle layeredHandle = view.features().publishLineLayers(
	    "layered-line",
	    BRLObolFeatureScope::Shared,
	    layers,
	    &unselectableStyle);
    SoNode *layeredNode = view.features().node(layeredHandle);
    SoBRLVListShape *layeredShape = static_cast<SoBRLVListShape *>(
					first_node_of_type(layeredNode, SoBRLVListShape::getClassTypeId()));
    if (!layeredHandle.isValid() || !layeredShape ||
	layeredShape->selectable.getValue())
	FAIL("parent selectable style should reach realized line-layer children");

    std::vector<int32_t> selectedPrimitives;
    selectedPrimitives.push_back(0);
    std::vector<int32_t> highlightedPrimitives;
    highlightedPrimitives.push_back(0);
    if (!view.features().replaceSelectedPrimitives(layeredHandle,
	    selectedPrimitives) ||
	!view.features().replaceHighlightedPrimitives(layeredHandle,
		highlightedPrimitives))
	FAIL("line-layer primitive state should be settable");
    if (!view.features().selectedPrimitives(layeredHandle,
					    selectedPrimitives) ||
	selectedPrimitives.size() != 1 ||
	selectedPrimitives[0] != 0 ||
	!view.features().highlightedPrimitives(layeredHandle,
		highlightedPrimitives) ||
	highlightedPrimitives.size() != 1 ||
	highlightedPrimitives[0] != 0)
	FAIL("line-layer primitive state should round-trip through feature store");
    layeredNode = view.features().node(layeredHandle);
    layeredShape = static_cast<SoBRLVListShape *>(
		       first_node_of_type(layeredNode, SoBRLVListShape::getClassTypeId()));
    if (!layeredShape ||
	layeredShape->selectedPrimitive.getNum() != 1 ||
	layeredShape->selectedPrimitive[0] != 0 ||
	layeredShape->highlightedPrimitive.getNum() != 1 ||
	layeredShape->highlightedPrimitive[0] != 0)
	FAIL("line-layer primitive state should reach realized VLIST children");
    if (!view.features().record(layeredHandle, overlayRecord) ||
	overlayRecord.selectedPrimitives.size() != 1 ||
	overlayRecord.highlightedPrimitives.size() != 1)
	FAIL("feature record should preserve primitive state");

    if (!view.features().replacePrimitiveMetadata(layeredHandle, 0,
	    primitiveMetadata) ||
	!view.features().resolvePrimitivePick("layered-line/red", 0,
		primitivePick) ||
	primitivePick.handle.id != layeredHandle.id ||
	strcmp(primitivePick.featureName.getString(),
	       "layered-line") != 0 ||
	primitivePick.primitiveIndex != 0 ||
	primitivePick.metadata.size() != 2 ||
	strcmp(primitivePick.metadata[0].key.getString(),
	       "overlap.objects") != 0 ||
	strcmp(primitivePick.metadata[0].value.getString(),
	       "box.s cone.s") != 0)
	FAIL("line-layer child primitive picks should resolve to parent primitive metadata");
    if (view.features().resolvePrimitivePick("layered-line/red", 1,
	    primitivePick))
	FAIL("line-layer child primitive resolver should reject out-of-range picks");

    SoSeparator *customNode = new SoSeparator;
    BRLObolFeatureHandle customHandle = view.features().publishCustomNode(
					    "custom-node",
					    BRLObolFeatureScope::Shared,
					    customNode,
					    &unselectableStyle,
					    &ownerA);
    if (!customHandle.isValid() ||
	view.features().node(customHandle) != customNode)
	FAIL("custom feature should publish the provider Coin node directly");
    BRLObolFeatureSummary customSummary;
    if (!view.features().summary("custom-node", customSummary) ||
	!customSummary.exists ||
	customSummary.kind != BRLObolFeatureKind::CustomNode ||
	!customSummary.realized)
	FAIL("custom feature summary should preserve custom-node kind");
    BRLObolFeatureRecord customRecord;
    if (!view.features().record(customHandle, customRecord) ||
	customRecord.kind != BRLObolFeatureKind::CustomNode ||
	customRecord.points.size() != 0 ||
	customRecord.commands.size() != 0 ||
	customRecord.owner.ownerToken != ownerA.ownerToken)
	FAIL("custom feature record should preserve owner and avoid fake geometry");
    if (!view.features().replaceMetadata(customHandle, metadata) ||
	!view.features().replaceSelectedPrimitives(customHandle,
		selectedPrimitives) ||
	!view.features().replaceHighlightedPrimitives(customHandle,
		highlightedPrimitives) ||
	!view.features().applyStyle(customHandle, unselectableStyle) ||
	view.features().node(customHandle) != customNode)
	FAIL("custom feature metadata/style/primitive updates should preserve the provider node");
    if (view.features().appendLinePoint(customHandle, SbVec3f(1.0f, 1.0f,
					1.0f)) ||
	view.features().replaceLabels(customHandle, labels))
	FAIL("custom feature should reject typed geometry mutation helpers");

    SoSeparator *replacementCustomNode = new SoSeparator;
    BRLObolFeatureHandle replacementCustomHandle =
	view.features().publishCustomNode("custom-node",
					  BRLObolFeatureScope::Shared, replacementCustomNode,
					  &unselectableStyle, &ownerA);
    if (!replacementCustomHandle.isValid() ||
	replacementCustomHandle.id != customHandle.id ||
	replacementCustomHandle.revision <= customHandle.revision ||
	view.features().node(replacementCustomHandle) !=
	replacementCustomNode)
	FAIL("custom feature replacement should preserve identity and replace the Coin node");
    if (!view.features().remove(replacementCustomHandle) ||
	view.features().exists("custom-node"))
	FAIL("custom feature removal should use normal feature-store teardown");

    struct edit_preview_callback_state previewState;
    previewState.revision = 42;
    previewState.updateCount = 0;
    BRLObolEditPreviewCallbacks previewCallbacks;
    previewCallbacks.previewContext = &previewState;
    previewCallbacks.revisionCallback = edit_preview_revision_cb;
    previewCallbacks.updateCallback = edit_preview_update_cb;

    std::vector<SbVec3f> previewPoints;
    previewPoints.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    previewPoints.push_back(SbVec3f(2.0f, 0.0f, 0.0f));
    BRLObolFeatureHandle previewHandle = view.features().publishEditPreview(
	    "edit-preview",
	    "edit-preview-source",
	    "intent-id",
	    "intent-role",
	    previewPoints,
	    lineCommands,
	    0,
	    0,
	    &previewCallbacks,
	    &ownerA);
    SoNode *previewNode = view.features().node(previewHandle);
    if (!previewHandle.isValid() || !previewNode ||
	!previewNode->isOfType(SoBRLEditPreview::getClassTypeId()))
	FAIL("publishEditPreview should realize a SoBRLEditPreview node");
    if (view.features().editPreviewRevision(previewHandle) != 42)
	FAIL("edit preview revision callback should be preserved");

    previewPoints.push_back(SbVec3f(2.0f, 2.0f, 0.0f));
    lineCommands.push_back(static_cast<int32_t>(BRLObolLineCommand::Draw));
    if (!view.features().replaceEditPreviewGeometry(previewHandle,
	    "edit-preview-replaced", previewPoints, lineCommands))
	FAIL("replaceEditPreviewGeometry should update preview geometry");
    BRLObolFeatureRecord previewRecord;
    if (!view.features().record(previewHandle, previewRecord) ||
	previewRecord.points.size() != 3 ||
	previewRecord.commands.size() != 3 ||
	previewRecord.kind != BRLObolFeatureKind::EditPreview)
	FAIL("replaceEditPreviewGeometry should preserve edit preview records");
    previewNode = view.features().node(previewHandle);
    if (!previewNode ||
	!previewNode->isOfType(SoBRLEditPreview::getClassTypeId()))
	FAIL("replaceEditPreviewGeometry should preserve the edit preview node");
    if (view.features().updateEditPreview(previewHandle) != 1 ||
	previewState.updateCount != 1)
	FAIL("replaceEditPreviewGeometry should preserve edit preview callbacks");

    return 0;
}

static int
test_selection_store(BRLObolViewController &view)
{
    BRLObolFeatureOwner ownerA;
    ownerA.ownerToken = (const void *)0x11;
    ownerA.ownerId = "selection-view-A";
    ownerA.ownerRole = "view";
    BRLObolFeatureOwner ownerB;
    ownerB.ownerToken = (const void *)0x22;
    ownerB.ownerId = "selection-view-B";
    ownerB.ownerRole = "view";

    if (!view.selection().addPath("all.g/box.s",
				  BRLOBOL_SELECTION_SELECTED_PATH, &ownerA))
	FAIL("selection addPath should accept owner A selected path");
    if (!view.selection().addPath("all.g/box.s",
				  BRLOBOL_SELECTION_HIGHLIGHTED_REF, &ownerA))
	FAIL("selection addPath should allow highlighted record with same path");
    if (!view.selection().addPath("all.g/cone.s",
				  BRLOBOL_SELECTION_SELECTED_PATH, &ownerB))
	FAIL("selection addPath should accept owner B selected path");

    if (view.selection().count() != 3 ||
	view.selection().count(&ownerA) != 2 ||
	view.selection().count(&ownerA,
			       BRLOBOL_SELECTION_SELECTED_PATH) != 1 ||
	view.selection().count(&ownerA,
			       BRLOBOL_SELECTION_HIGHLIGHTED_REF) != 1)
	FAIL("selection count should filter by owner and kind");
    if (!view.selection().containsPath("all.g/box.s",
				       BRLOBOL_SELECTION_SELECTED_PATH, &ownerA))
	FAIL("selection containsPath should find owner A selected path");
    if (view.selection().containsPath("all.g/box.s",
				      BRLOBOL_SELECTION_SELECTED_PATH, &ownerB))
	FAIL("selection containsPath should not cross owner scopes");

    struct selection_visit_state visit = {0, NULL};
    view.selection().visitPaths(count_selection_path_cb, &visit, &ownerA,
				BRLOBOL_SELECTION_SELECTED_PATH);
    if (visit.count != 1 || !visit.last ||
	strcmp(visit.last, "all.g/box.s") != 0)
	FAIL("selection visitPaths should filter to owner A selected path");

    view.selection().clear(&ownerA, BRLOBOL_SELECTION_HIGHLIGHTED_REF);
    if (view.selection().count(&ownerA) != 1 ||
	view.selection().containsPath("all.g/box.s",
				      BRLOBOL_SELECTION_HIGHLIGHTED_REF, &ownerA))
	FAIL("selection clear should remove only matching highlighted records");

    std::vector<BRLObolSelectionRecord> picks;
    BRLObolSelectionRecord pick;
    pick.path = "all.g/sphere.s";
    pick.kind = BRLOBOL_SELECTION_SELECTED_PATH;
    pick.hitDistance = 3.0;
    picks.push_back(pick);
    if (!view.selection().applyPickResults(picks, NULL, NULL, &ownerA))
	FAIL("selection applyPickResults should succeed");
    if (view.selection().count(&ownerA) != 1 ||
	!view.selection().containsPath("all.g/sphere.s",
				       BRLOBOL_SELECTION_SELECTED_PATH, &ownerA) ||
	view.selection().count(&ownerB) != 1)
	FAIL("selection applyPickResults should replace only owner A records");

    return 0;
}

static int
test_polygon_nodes_and_sketch(BRLObolViewController &view)
{
    plane_t viewPlane;
    HSET(viewPlane, 0.0, 0.0, 1.0, 0.0);

    BRLObolPolygonHandle handle = view.polygons().create(
				      "poly",
				      BRLObolFeatureScope::Shared,
				      BRLObolPolygonType::Square,
				      SbVec3f(0.0f, 0.0f, 0.0f),
				      viewPlane,
				      0.0f);
    if (!handle.isValid())
	FAIL("polygon create should return a valid handle");
    if (!view.polygons().updateScreenPoint(handle, 10, 8,
					   BRLObolPolygonUpdate::Default))
	FAIL("square polygon update should succeed");
    if (!view.polygons().setFill(handle, TRUE, SbVec2f(1.0f, 0.0f), 2.0f))
	FAIL("polygon setFill should succeed");
    if (!view.polygons().setFillColor(handle, SbColor(0.2f, 0.4f, 0.8f)))
	FAIL("polygon setFillColor should succeed");
    if (!view.polygons().setEdgeColor(handle, SbColor(1.0f, 0.0f, 0.0f)))
	FAIL("polygon setEdgeColor should succeed");

    SoNode *node = view.polygons().node(handle);
    if (!node || !node->isOfType(SoGroup::getClassTypeId()))
	FAIL("filled polygon should realize as a group");
    if (count_nodes_of_type(node, SoBRLMeshShape::getClassTypeId()) != 0)
	FAIL("legacy polygon setFill should not imply mesh fill");
    if (count_nodes_of_type(node, SoBRLVListShape::getClassTypeId()) != 2)
	FAIL("legacy polygon setFill should realize hatch and outline children");
    if (!view.polygons().setFillFlags(handle,
				      BRLOBOL_POLYGON_FILL_HATCH | BRLOBOL_POLYGON_FILL_MESH))
	FAIL("polygon setFillFlags should succeed");
    node = view.polygons().node(handle);
    if (count_nodes_of_type(node, SoBRLMeshShape::getClassTypeId()) != 1)
	FAIL("mesh fill flag should add one mesh fill child");
    if (count_nodes_of_type(node, SoBRLVListShape::getClassTypeId()) != 2)
	FAIL("hatch plus mesh fill should keep hatch and outline children");

    BRLObolPolygonHandle csgA = view.polygons().create(
				    "csg-a",
				    BRLObolFeatureScope::Shared,
				    BRLObolPolygonType::Circle,
				    SbVec3f(0.0f, 0.0f, 0.0f),
				    viewPlane,
				    0.0f);
    BRLObolPolygonHandle csgB = view.polygons().create(
				    "csg-b",
				    BRLObolFeatureScope::Shared,
				    BRLObolPolygonType::Ellipse,
				    SbVec3f(0.5f, 0.0f, 0.0f),
				    viewPlane,
				    0.0f);
    if (!csgA.isValid() || !csgB.isValid() ||
	!view.polygons().updateModelPoint(csgA,
					  SbVec3f(1.0f, 0.0f, 0.0f),
					  BRLObolPolygonUpdate::Default) ||
	!view.polygons().updateModelPoint(csgB,
					  SbVec3f(1.5f, 0.75f, 0.0f),
					  BRLObolPolygonUpdate::Default) ||
	!view.polygons().csg(csgA, csgB, bg_Union))
	FAIL("polygon CSG union should succeed for overlapping Obol polygons");
    BRLObolPolygonRecord csgRecord;
    if (!view.polygons().record(csgA, csgRecord) ||
	csgRecord.pointCount == 0 ||
	csgRecord.type != BRLObolPolygonType::General)
	FAIL("polygon CSG union should preserve non-empty general geometry");

    char dbpath[MAXPATHLEN];
    bu_dir(dbpath, MAXPATHLEN, BU_DIR_CURR,
	   "brlobol_view_store_test.g", NULL);
    bu_file_delete(dbpath);

    struct rt_wdb *wdbp = wdb_fopen_v(dbpath, 5);
    if (!wdbp)
	FAIL("wdb_fopen_v should succeed");
    struct db_i *dbip = wdbp->dbip;

    if (!view.polygons().exportSketch(handle, dbip, "poly.s")) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("polygon exportSketch should succeed");
    }

    struct directory *dp = db_lookup(dbip, "poly.s", LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("exported sketch should be in the database directory");
    }

    BRLObolPolygonStore imported;
    BRLObolPolygonHandle importedHandle = imported.importSketch(
	    "imported",
	    BRLObolFeatureScope::Shared,
	    dbip,
	    dp);
    if (!importedHandle.isValid()) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("polygon importSketch should succeed");
    }

    BRLObolPolygonRecord record;
    if (!imported.record(importedHandle, record)) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("imported polygon record lookup should succeed");
    }
    if (record.pointCount != 4 || record.contourCount != 1) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("imported polygon should preserve square contour geometry");
    }
    if (!record.fill ||
	record.fillFlags != BRLOBOL_POLYGON_FILL_HATCH ||
	record.fillSpacing != 2.0f ||
	record.type != BRLObolPolygonType::Square) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("imported polygon should preserve visual/type attributes");
    }

    wdb_close(wdbp);
    bu_file_delete(dbpath);
    return 0;
}

int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);
    if (argc != 1)
	FAIL("unexpected arguments");

    brlobol_init(NULL);

    SoSeparator *root = new SoSeparator;
    root->ref();
    {
	BRLObolViewController view(root);
	if (test_feature_nodes(view) != 0) {
	    root->unref();
	    return 1;
	}
	if (test_selection_store(view) != 0) {
	    root->unref();
	    return 1;
	}
	if (test_polygon_nodes_and_sketch(view) != 0) {
	    root->unref();
	    return 1;
	}
    }
    root->unref();

    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
