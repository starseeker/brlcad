/*                    T E S T _ V I E W _ S T O R E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol.h"
#include "BObol/BEditPreview.h"
#include "BObol/BHUDLabelOverlay.h"
#include "BObol/BMeshShape.h"
#include "BObol/BVListShape.h"
#include "BObol/BEditPreview.h"
#include "BObol/BHUDLabelOverlay.h"
#include "bu/app.h"
#include "bu/file.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/SoType.h>
#include <Inventor/SoPath.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
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

static int
count_feature_visit_cb(const BObolFeatureRecord &record, void *userData)
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
test_feature_nodes(BObolViewController &view)
{
    std::vector<BObolLabel> labels;
    BObolLabel label;
    label.text = "store label";
    label.point = SbVec3f(1.0f, 2.0f, 0.0f);
    label.hasLeader = TRUE;
    label.target = SbVec3f(0.0f, 0.0f, 0.0f);
    label.fontSize = 17.0f;
    labels.push_back(label);

    BObolFeatureHandle labelHandle = view.features().publishLabels(
					   "labels",
					   BObolFeatureScope::Shared,
					   labels);
    if (!labelHandle.isValid())
	FAIL("publishLabels should return a valid handle");
    if (count_nodes_of_type(view.features().node(labelHandle),
			    SoText2::getClassTypeId()) != 1)
	FAIL("label feature should realize a SoText2 node");
    SoFont *font = static_cast<SoFont *>(first_node_of_type(
	    view.features().node(labelHandle), SoFont::getClassTypeId()));
    if (!font || !NEAR_EQUAL(font->size.getValue(), 17.0f, SMALL_FASTF))
	FAIL("label feature should preserve explicit font size");

    std::vector<BObolLabel> hudLabels;
    BObolLabel hudLabel;
    hudLabel.text = "store hud label";
    hudLabel.point = SbVec3f(8.0f, 58.0f, 0.0f);
    hudLabel.hasColor = TRUE;
    hudLabel.color = SbColor(1.0f, 1.0f, 0.0f);
    hudLabel.fontSize = 12.0f;
    hudLabel.sourceId = 42;
    hudLabels.push_back(hudLabel);

    BObolFeatureHandle hudHandle = view.features().publishHudLabels(
					 "hud-labels",
					 BObolFeatureScope::Shared,
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
    if (!NEAR_EQUAL(hud->position.getValue()[0], 8.0f, SMALL_FASTF) ||
	!NEAR_EQUAL(hud->position.getValue()[1], 58.0f, SMALL_FASTF) ||
	!NEAR_EQUAL(hud->fontSize.getValue(), 12.0f, SMALL_FASTF))
	FAIL("HUD label feature should preserve screen position and font size");

    /* A frame-rate/HUD synchronizer republishes its current value on a
     * periodic cadence.  Equal typed content must retain both node and handle
     * identity and must not request another frame; otherwise the faceplate can
     * starve a capacity-relevant LoD confirmation indefinitely. */
    view.clearRenderRequest();
    SoNode *hudNode = view.features().node(hudHandle);
    BObolFeatureHandle unchangedHud = view.features().publishHudLabels(
	"hud-labels", BObolFeatureScope::Shared, hudLabels);
    if (unchangedHud != hudHandle || view.features().node(unchangedHud) !=
	    hudNode || view.isRenderRequested())
	FAIL("equal HUD publication should be an idempotent retained update");
    hudLabels[0].text = "changed store hud label";
    BObolFeatureHandle changedHud = view.features().publishHudLabels(
	"hud-labels", BObolFeatureScope::Shared, hudLabels);
    if (changedHud.id != hudHandle.id ||
	changedHud.revision <= hudHandle.revision ||
	view.features().node(changedHud) == hudNode ||
	!view.isRenderRequested())
	FAIL("changed HUD publication should replace the retained presentation");
    view.clearRenderRequest();

    /* SoHUDKit updates viewportSize while rendering.  That runtime projection
     * state must not turn an otherwise identical retained line publication
     * into a new feature after a canvas resize. */
    BObolFeatureStyle hudLineStyle;
    hudLineStyle.hud = TRUE;
    std::vector<SbVec3f> hudLinePoints;
    hudLinePoints.push_back(SbVec3f(4.0f, 8.0f, 0.0f));
    hudLinePoints.push_back(SbVec3f(64.0f, 8.0f, 0.0f));
    std::vector<int32_t> hudLineCommands;
    hudLineCommands.push_back(static_cast<int32_t>(BObolLineCommand::Move));
    hudLineCommands.push_back(static_cast<int32_t>(BObolLineCommand::Draw));
    BObolFeatureHandle hudLineHandle = view.features().publishLineSet(
	"hud-line", BObolFeatureScope::Shared, hudLinePoints,
	hudLineCommands, &hudLineStyle);
    SoNode *hudLineNode = view.features().node(hudLineHandle);
    SoHUDKit *hudLine = static_cast<SoHUDKit *>(
	first_node_of_type(hudLineNode, SoHUDKit::getClassTypeId()));
    if (!hudLine)
	FAIL("HUD line feature should realize a SoHUDKit");
    hudLine->viewportSize = SbVec2f(947.0f, 693.0f);
    view.clearRenderRequest();
    BObolFeatureHandle unchangedHudLine = view.features().publishLineSet(
	"hud-line", BObolFeatureScope::Shared, hudLinePoints,
	hudLineCommands, &hudLineStyle);
    if (unchangedHudLine != hudLineHandle ||
	view.features().node(unchangedHudLine) != hudLineNode ||
	view.isRenderRequested())
	FAIL("HUD runtime viewport state should not invalidate equal publication");

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

    BObolFeatureHandle meshHandle = view.features().publishIndexedFaceSet(
					  "surface",
					  BObolFeatureScope::Shared,
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

    BObolFeatureOwner ownerA;
    ownerA.ownerToken = (const void *)0x1;
    ownerA.ownerId = "view-A";
    ownerA.ownerRole = "view";
    BObolFeatureOwner ownerB;
    ownerB.ownerToken = (const void *)0x2;
    ownerB.ownerId = "view-B";
    ownerB.ownerRole = "view";

    std::vector<int32_t> lineCommands;
    lineCommands.push_back(static_cast<int32_t>(BObolLineCommand::Move));
    lineCommands.push_back(static_cast<int32_t>(BObolLineCommand::Draw));
    BObolFeatureHandle localA = view.features().publishLineSet(
				      "local-line",
				      BObolFeatureScope::Local,
				      points,
				      lineCommands,
				      NULL,
				      &ownerA);
    BObolFeatureHandle localB = view.features().publishLineSet(
				      "local-line",
				      BObolFeatureScope::Local,
				      points,
				      lineCommands,
				      NULL,
				      &ownerB);
    if (!localA.isValid() || !localB.isValid() || localA.id == localB.id)
	FAIL("owner-scoped local feature names should not collide");
    if (!view.features().existsOwned("local-line",
				     BOBOL_FEATURE_SCOPE_LOCAL, &ownerA))
	FAIL("owner A should find its local feature");
    if (!view.features().existsOwned("local-line",
				     BOBOL_FEATURE_SCOPE_LOCAL, &ownerB))
	FAIL("owner B should find its local feature");

    struct feature_visit_count visitA = {0};
    view.features().visitRecords(count_feature_visit_cb, &visitA,
				 BOBOL_FEATURE_SCOPE_LOCAL, &ownerA);
    if (visitA.count != 1)
	FAIL("owner-scoped feature visit should see only one local record");

    BObolFeatureRecord localRecord;
    if (!view.features().record(localA, localRecord) ||
	localRecord.points.size() != points.size() ||
	localRecord.commands.size() != lineCommands.size() ||
	localRecord.owner.ownerToken != ownerA.ownerToken)
	FAIL("feature record snapshot should preserve geometry and owner");

    BObolOverlayInfo overlay;
    overlay.isOverlay = TRUE;
    overlay.ownerToken = ownerA.ownerToken;
    overlay.role = BObolOverlayRole::Model;
    overlay.overlayClass = BObolOverlayClass::TclOverlay;
    overlay.lifecycle = BObolOverlayLifecycle::PerCommand;
    overlay.order = BObolOverlayOrder::PostTransparent;
    overlay.sortOrder = 7;
    overlay.sourcePath = "local-line";
    if (!view.features().setOverlayInfo(localA, overlay))
	FAIL("feature overlay metadata should be settable");

    BObolOverlayInfo overlayRead;
    if (!view.features().overlayInfo(localA, overlayRead) ||
	!overlayRead.isOverlay ||
	overlayRead.ownerToken != ownerA.ownerToken ||
	overlayRead.role != BObolOverlayRole::Model ||
	overlayRead.overlayClass != BObolOverlayClass::TclOverlay ||
	overlayRead.lifecycle != BObolOverlayLifecycle::PerCommand ||
	overlayRead.order != BObolOverlayOrder::PostTransparent ||
	overlayRead.sortOrder != 7 ||
	bu_strcmp(overlayRead.sourcePath.getString(), "local-line") != 0)
	FAIL("feature overlay metadata should round-trip through overlayInfo");

    BObolFeatureSummary overlaySummary;
    if (!view.features().summaryOwned("local-line", overlaySummary,
				      BOBOL_FEATURE_SCOPE_LOCAL, &ownerA) ||
	!overlaySummary.exists ||
	!overlaySummary.overlay.isOverlay ||
	overlaySummary.overlay.overlayClass !=
	BObolOverlayClass::TclOverlay)
	FAIL("feature summary should preserve overlay metadata");

    BObolFeatureRecord overlayRecord;
    if (!view.features().record(localA, overlayRecord) ||
	!overlayRecord.overlay.isOverlay ||
	overlayRecord.overlay.lifecycle !=
	BObolOverlayLifecycle::PerCommand)
	FAIL("feature record should preserve overlay metadata");

    std::vector<BObolFeatureMetadata> metadata;
    BObolFeatureMetadata metadataItem;
    metadataItem.key = "result.kind";
    metadataItem.value = "diagnostic-line";
    metadata.push_back(metadataItem);
    metadataItem.key = "result.index";
    metadataItem.value = "7";
    metadata.push_back(metadataItem);
    if (!view.features().replaceMetadata(localA, metadata))
	FAIL("feature metadata should be settable");
    std::vector<BObolFeatureMetadata> metadataRead;
    if (!view.features().metadata(localA, metadataRead) ||
	metadataRead.size() != 2 ||
	bu_strcmp(metadataRead[0].key.getString(), "result.kind") != 0 ||
	bu_strcmp(metadataRead[0].value.getString(),
	       "diagnostic-line") != 0 ||
	bu_strcmp(metadataRead[1].key.getString(), "result.index") != 0 ||
	bu_strcmp(metadataRead[1].value.getString(), "7") != 0)
	FAIL("feature metadata should round-trip through metadata readback");

    if (!view.features().summaryOwned("local-line", overlaySummary,
				      BOBOL_FEATURE_SCOPE_LOCAL, &ownerA) ||
	overlaySummary.metadataCount != 2)
	FAIL("feature summary should report metadata count");

    if (!view.features().record(localA, overlayRecord) ||
	overlayRecord.metadata.size() != 2 ||
	bu_strcmp(overlayRecord.metadata[0].key.getString(),
	       "result.kind") != 0)
	FAIL("feature record should preserve metadata");

    std::vector<BObolFeatureMetadata> primitiveMetadata;
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
	bu_strcmp(metadataRead[0].key.getString(), "overlap.objects") != 0 ||
	bu_strcmp(metadataRead[0].value.getString(),
	       "box.s cone.s") != 0 ||
	bu_strcmp(metadataRead[1].key.getString(), "overlap.depth") != 0 ||
	bu_strcmp(metadataRead[1].value.getString(), "0.25") != 0)
	FAIL("feature primitive metadata should round-trip through readback");
    if (!view.features().summaryOwned("local-line", overlaySummary,
				      BOBOL_FEATURE_SCOPE_LOCAL, &ownerA) ||
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
	bu_strcmp(static_cast<SoBRLVListShape *>(localNode)->
	       recordRole.getValue().getString(), "view-feature") != 0)
	FAIL("feature overlay metadata should preserve feature-store node role");

    if (!view.features().clearOverlayInfo(localA) ||
	!view.features().overlayInfo(localA, overlayRead) ||
	overlayRead.isOverlay)
	FAIL("feature overlay metadata should be clearable");
    localNode = view.features().node(localA);
    if (!localNode ||
	!localNode->isOfType(SoBRLVListShape::getClassTypeId()) ||
	bu_strcmp(static_cast<SoBRLVListShape *>(localNode)->
	       recordRole.getValue().getString(), "view-feature") != 0)
	FAIL("cleared feature overlay metadata should restore view-feature node role");

    localA = view.features().publishLineSet("local-line",
					    BObolFeatureScope::Local, points, lineCommands, NULL, &ownerA);
    if (!localA.isValid() ||
	!view.features().metadata(localA, metadataRead) ||
	!metadataRead.empty() ||
	!view.features().primitiveMetadata(localA, 0, metadataRead) ||
	!metadataRead.empty())
	FAIL("feature publish should clear stale metadata");

    BObolFeatureStyle unselectableStyle;
    unselectableStyle.hasSelectable = TRUE;
    unselectableStyle.selectable = FALSE;
    if (!view.features().applyStyle(localA, unselectableStyle))
	FAIL("feature selectable style should be applicable");
    localNode = view.features().node(localA);
    if (!localNode ||
	!localNode->isOfType(SoBRLVListShape::getClassTypeId()) ||
	static_cast<SoBRLVListShape *>(localNode)->selectable.getValue())
	FAIL("feature selectable style should control realized line shapes");

    BObolFeaturePrimitivePick primitivePick;
    if (!view.features().resolvePrimitivePick("local-line", 0,
	    primitivePick, BOBOL_FEATURE_SCOPE_LOCAL, &ownerA) ||
	primitivePick.handle.id != localA.id ||
	bu_strcmp(primitivePick.featureName.getString(), "local-line") != 0 ||
	primitivePick.primitiveIndex != 0)
	FAIL("feature primitive pick resolver should find direct local features");

    BObolLineLayer layer;
    layer.name = "layered-line/red";
    layer.points = points;
    layer.commands = lineCommands;
    std::vector<BObolLineLayer> layers;
    layers.push_back(layer);
    BObolFeatureHandle layeredHandle = view.features().publishLineLayers(
	    "layered-line",
	    BObolFeatureScope::Shared,
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
	bu_strcmp(primitivePick.featureName.getString(),
	       "layered-line") != 0 ||
	primitivePick.primitiveIndex != 0 ||
	primitivePick.metadata.size() != 2 ||
	bu_strcmp(primitivePick.metadata[0].key.getString(),
	       "overlap.objects") != 0 ||
	bu_strcmp(primitivePick.metadata[0].value.getString(),
	       "box.s cone.s") != 0)
	FAIL("line-layer child primitive picks should resolve to parent primitive metadata");
    if (view.features().resolvePrimitivePick("layered-line/red", 1,
	    primitivePick))
	FAIL("line-layer child primitive resolver should reject out-of-range picks");

    SoSeparator *customNode = new SoSeparator;
    BObolFeatureHandle customHandle = view.features().publishCustomNode(
					    "custom-node",
					    BObolFeatureScope::Shared,
					    customNode,
					    &unselectableStyle,
					    &ownerA);
    if (!customHandle.isValid() ||
	view.features().node(customHandle) != customNode)
	FAIL("custom feature should publish the provider Coin node directly");
    BObolFeatureSummary customSummary;
    if (!view.features().summary("custom-node", customSummary) ||
	!customSummary.exists ||
	customSummary.kind != BObolFeatureKind::CustomNode ||
	!customSummary.realized)
	FAIL("custom feature summary should preserve custom-node kind");
    BObolFeatureRecord customRecord;
    if (!view.features().record(customHandle, customRecord) ||
	customRecord.kind != BObolFeatureKind::CustomNode ||
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
    BObolFeatureHandle replacementCustomHandle =
	view.features().publishCustomNode("custom-node",
					  BObolFeatureScope::Shared, replacementCustomNode,
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

    std::vector<SbVec3f> previewPoints;
    previewPoints.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    previewPoints.push_back(SbVec3f(2.0f, 0.0f, 0.0f));
    BObolFeatureHandle previewHandle = view.features().publishEditPreview(
	    "edit-preview",
	    "edit-preview-source",
	    "intent-id",
	    "intent-role",
	    previewPoints,
	    lineCommands,
	    0,
	    0,
	    &ownerA);
    SoNode *previewNode = view.features().node(previewHandle);
    if (!previewHandle.isValid() || !previewNode ||
	!previewNode->isOfType(SoBRLEditPreview::getClassTypeId()))
	FAIL("publishEditPreview should realize a SoBRLEditPreview node");

    previewPoints.push_back(SbVec3f(2.0f, 2.0f, 0.0f));
    lineCommands.push_back(static_cast<int32_t>(BObolLineCommand::Draw));
    if (!view.features().replaceEditPreviewGeometry(previewHandle,
	    "edit-preview-replaced", previewPoints, lineCommands))
	FAIL("replaceEditPreviewGeometry should update preview geometry");
    BObolFeatureRecord previewRecord;
    if (!view.features().record(previewHandle, previewRecord) ||
	previewRecord.points.size() != 3 ||
	previewRecord.commands.size() != 3 ||
	previewRecord.kind != BObolFeatureKind::EditPreview)
	FAIL("replaceEditPreviewGeometry should preserve edit preview records");
    previewNode = view.features().node(previewHandle);
    if (!previewNode ||
	!previewNode->isOfType(SoBRLEditPreview::getClassTypeId()))
	FAIL("replaceEditPreviewGeometry should preserve the edit preview node");
    return 0;
}

static int
test_selection_store(BObolViewController &view)
{
    BObolFeatureOwner ownerA;
    ownerA.ownerToken = (const void *)0x11;
    ownerA.ownerId = "selection-view-A";
    ownerA.ownerRole = "view";
    BObolFeatureOwner ownerB;
    ownerB.ownerToken = (const void *)0x22;
    ownerB.ownerId = "selection-view-B";
    ownerB.ownerRole = "view";

    if (!view.selection().addPath("all.g/box.s",
				  BOBOL_SELECTION_SELECTED_PATH, &ownerA))
	FAIL("selection addPath should accept owner A selected path");
    if (!view.selection().addPath("all.g/box.s",
				  BOBOL_SELECTION_HIGHLIGHTED_REF, &ownerA))
	FAIL("selection addPath should allow highlighted record with same path");
    if (!view.selection().addPath("all.g/cone.s",
				  BOBOL_SELECTION_SELECTED_PATH, &ownerB))
	FAIL("selection addPath should accept owner B selected path");

    if (view.selection().count() != 3 ||
	view.selection().count(&ownerA) != 2 ||
	view.selection().count(&ownerA,
			       BOBOL_SELECTION_SELECTED_PATH) != 1 ||
	view.selection().count(&ownerA,
			       BOBOL_SELECTION_HIGHLIGHTED_REF) != 1)
	FAIL("selection count should filter by owner and kind");
    if (!view.selection().containsPath("all.g/box.s",
				       BOBOL_SELECTION_SELECTED_PATH, &ownerA))
	FAIL("selection containsPath should find owner A selected path");
    if (view.selection().containsPath("all.g/box.s",
				      BOBOL_SELECTION_SELECTED_PATH, &ownerB))
	FAIL("selection containsPath should not cross owner scopes");

    struct selection_visit_state visit = {0, NULL};
    view.selection().visitPaths(count_selection_path_cb, &visit, &ownerA,
				BOBOL_SELECTION_SELECTED_PATH);
    if (visit.count != 1 || !visit.last ||
	bu_strcmp(visit.last, "all.g/box.s") != 0)
	FAIL("selection visitPaths should filter to owner A selected path");

    view.selection().clear(&ownerA, BOBOL_SELECTION_HIGHLIGHTED_REF);
    if (view.selection().count(&ownerA) != 1 ||
	view.selection().containsPath("all.g/box.s",
				      BOBOL_SELECTION_HIGHLIGHTED_REF, &ownerA))
	FAIL("selection clear should remove only matching highlighted records");

    std::vector<BObolSelectionRecord> picks;
    BObolSelectionRecord pick;
    pick.path = "all.g/sphere.s";
    pick.kind = BOBOL_SELECTION_SELECTED_PATH;
    pick.hitDistance = 3.0;
    picks.push_back(pick);
    if (!view.selection().applyPickResults(picks, NULL, NULL, &ownerA))
	FAIL("selection applyPickResults should succeed");
    if (view.selection().count(&ownerA) != 1 ||
	!view.selection().containsPath("all.g/sphere.s",
				       BOBOL_SELECTION_SELECTED_PATH, &ownerA) ||
	view.selection().count(&ownerB) != 1)
	FAIL("selection applyPickResults should replace only owner A records");

    return 0;
}

static int
test_polygon_nodes_and_sketch(BObolViewController &view)
{
    plane_t viewPlane;
    HSET(viewPlane, 0.0, 0.0, 1.0, 0.0);

    BObolPolygonHandle handle = view.polygons().create(
				      "poly",
				      BObolFeatureScope::Shared,
				      BObolPolygonType::Square,
				      SbVec3f(0.0f, 0.0f, 0.0f),
				      viewPlane,
				      0.0f);
    if (!handle.isValid())
	FAIL("polygon create should return a valid handle");
    if (!view.polygons().updateScreenPoint(handle, 10, 8,
					   BObolPolygonUpdate::Default))
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

    SoSearchAction outlineSearch;
    outlineSearch.setType(SoBRLVListShape::getClassTypeId());
    outlineSearch.setInterest(SoSearchAction::LAST);
    outlineSearch.apply(node);
    SoPath *outlinePath = outlineSearch.getPath();
    SoBRLVListShape *outline = outlinePath ?
	static_cast<SoBRLVListShape *>(outlinePath->getTail()) : NULL;
    if (!outline || outline->getSegmentCount() != 4)
	FAIL("closed square polygon outline should realize all four edges");

    if (!view.polygons().move(handle, SbVec3f(2.0f, 3.0f, 0.0f),
		SbVec3f(0.0f, 0.0f, 0.0f)))
	FAIL("square polygon move should succeed");
    node = view.polygons().node(handle);
    SoSearchAction movedOutlineSearch;
    movedOutlineSearch.setType(SoBRLVListShape::getClassTypeId());
    movedOutlineSearch.setInterest(SoSearchAction::LAST);
    movedOutlineSearch.apply(node);
    SoPath *movedOutlinePath = movedOutlineSearch.getPath();
    SoBRLVListShape *movedOutline = movedOutlinePath ?
	static_cast<SoBRLVListShape *>(movedOutlinePath->getTail()) : NULL;
    if (!movedOutline || movedOutline->getSegmentCount() != 4)
	FAIL("moving a constrained polygon should preserve all four edges");

    if (!view.polygons().setAllContoursOpen(handle, TRUE))
	FAIL("square polygon compatibility open update should succeed");
    node = view.polygons().node(handle);
    SoSearchAction constrainedOutlineSearch;
    constrainedOutlineSearch.setType(SoBRLVListShape::getClassTypeId());
    constrainedOutlineSearch.setInterest(SoSearchAction::LAST);
    constrainedOutlineSearch.apply(node);
    SoPath *constrainedOutlinePath = constrainedOutlineSearch.getPath();
    SoBRLVListShape *constrainedOutline = constrainedOutlinePath ?
	static_cast<SoBRLVListShape *>(constrainedOutlinePath->getTail()) : NULL;
    if (!constrainedOutline || constrainedOutline->getSegmentCount() != 4)
	FAIL("constrained square outline should remain closed while editing");
    if (!view.polygons().setAllContoursOpen(handle, FALSE))
	FAIL("square polygon close restore should succeed");
    if (!view.polygons().setFillFlags(handle,
				      BOBOL_POLYGON_FILL_HATCH | BOBOL_POLYGON_FILL_MESH))
	FAIL("polygon setFillFlags should succeed");
    node = view.polygons().node(handle);
    if (count_nodes_of_type(node, SoBRLMeshShape::getClassTypeId()) != 1)
	FAIL("mesh fill flag should add one mesh fill child");
    if (count_nodes_of_type(node, SoBRLVListShape::getClassTypeId()) != 2)
	FAIL("hatch plus mesh fill should keep hatch and outline children");

    BObolPolygonHandle csgA = view.polygons().create(
				    "csg-a",
				    BObolFeatureScope::Shared,
				    BObolPolygonType::Circle,
				    SbVec3f(0.0f, 0.0f, 0.0f),
				    viewPlane,
				    0.0f);
    BObolPolygonHandle csgB = view.polygons().create(
				    "csg-b",
				    BObolFeatureScope::Shared,
				    BObolPolygonType::Ellipse,
				    SbVec3f(0.5f, 0.0f, 0.0f),
				    viewPlane,
				    0.0f);
    if (!csgA.isValid() || !csgB.isValid() ||
	!view.polygons().updateModelPoint(csgA,
					  SbVec3f(1.0f, 0.0f, 0.0f),
					  BObolPolygonUpdate::Default) ||
	!view.polygons().updateModelPoint(csgB,
					  SbVec3f(1.5f, 0.75f, 0.0f),
					  BObolPolygonUpdate::Default) ||
	!view.polygons().csg(csgA, csgB, BG_POLYGON_BOOLEAN_UNION))
	FAIL("polygon CSG union should succeed for overlapping Obol polygons");
    BObolPolygonRecord csgRecord;
    if (!view.polygons().record(csgA, csgRecord) ||
	csgRecord.pointCount == 0 ||
	csgRecord.type != BObolPolygonType::General)
	FAIL("polygon CSG union should preserve non-empty general geometry");

    char dbpath[MAXPATHLEN];
    bu_dir(dbpath, MAXPATHLEN, BU_DIR_CURR,
	   "bobol_view_store_test.g", NULL);
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

    BObolPolygonStore imported;
    BObolPolygonHandle importedHandle = imported.importSketch(
	    "imported",
	    BObolFeatureScope::Shared,
	    dbip,
	    dp);
    if (!importedHandle.isValid()) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("polygon importSketch should succeed");
    }

    BObolPolygonRecord record;
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
	record.fillFlags != BOBOL_POLYGON_FILL_HATCH ||
	!NEAR_EQUAL(record.fillSpacing, 2.0f, SMALL_FASTF) ||
	record.type != BObolPolygonType::Square) {
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

    bobol_init(NULL);

    SoSeparator *root = new SoSeparator;
    root->ref();
    {
	BObolViewController view(root);
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
