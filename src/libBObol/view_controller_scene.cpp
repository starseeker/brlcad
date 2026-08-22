/*             V I E W _ C O N T R O L L E R _ S C E N E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file view_controller_scene.cpp
 *
 * Retained scene, feature-store, and edit-presentation forwarding.  View/LoD
 * policy remains in view_controller.cpp; this unit performs bounded semantic
 * mutations through the existing owners.
 */

#include "common.h"

#include "BObol/BEditPreview.h"
#include "BObol/BHUDLabelOverlay.h"
#include "BObol/BLineLayerOverlay.h"
#include "BObol/BSceneGroup.h"
#include "bu/str.h"
#include "view_controller_private.h"

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

BObolSceneController *
BObolViewController::getSceneController(void)
{
    return &this->d->sceneController;
}

const BObolSceneController *
BObolViewController::getSceneController(void) const
{
    return &this->d->sceneController;
}

BObolFeatureStore &
BObolViewController::features(void)
{
    return *this->d->featureStore;
}

const BObolFeatureStore &
BObolViewController::features(void) const
{
    return *this->d->featureStore;
}

BObolPolygonStore &
BObolViewController::polygons(void)
{
    return *this->d->polygonStore;
}

const BObolPolygonStore &
BObolViewController::polygons(void) const
{
    return *this->d->polygonStore;
}

BObolSelectionStore &
BObolViewController::selection(void)
{
    return *this->d->selectionStore;
}

const BObolSelectionStore &
BObolViewController::selection(void) const
{
    return *this->d->selectionStore;
}

SoViewport *
BObolViewController::getViewport(void)
{
    return this->d->viewport;
}

const SoViewport *
BObolViewController::getViewport(void) const
{
    return this->d->viewport;
}

void
BObolViewController::setRenderContextManager(SoDB::ContextManager *manager)
{
    SoGLRenderAction *action = this->d->renderManager ?
	this->d->renderManager->getGLRenderAction() : NULL;
    if (!action || action->getContextManager() == manager)
	return;
    /* The cached image renderer owns provider-specific context state.  Drop
     * it while the old provider is still alive rather than retaining a stale
     * provider through a host switch. */
    delete this->d->imageRenderer;
    this->d->imageRenderer = NULL;
    this->d->imageRendererManager = NULL;
    this->invalidateRendererPerformanceHistory();
    action->setContextManager(manager);
}

SoDB::ContextManager *
BObolViewController::getRenderContextManager(void) const
{
    SoGLRenderAction *action = this->d->renderManager ?
	this->d->renderManager->getGLRenderAction() : NULL;
    return action ? action->getContextManager() : NULL;
}

SoRenderManager *
BObolViewController::getRenderManager(void)
{
    return this->d->renderManager;
}

const SoRenderManager *
BObolViewController::getRenderManager(void) const
{
    return this->d->renderManager;
}

static int
find_edit_preview_child(SoGroup *group, const char *previewId)
{
    if (!group || !previewId)
	return -1;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(SoBRLEditPreview::getClassTypeId()))
	    continue;
	SoBRLEditPreview *preview = static_cast<SoBRLEditPreview *>(node);
	if (bu_strcmp(preview->previewId.getValue().getString(), previewId) == 0)
	    return i;
    }
    return -1;
}

static int
find_hud_label_overlay_child(SoGroup *group, const char *labelId)
{
    if (!group || !labelId)
	return -1;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(SoBRLHUDLabelOverlay::getClassTypeId()))
	    continue;
	SoBRLHUDLabelOverlay *label = static_cast<SoBRLHUDLabelOverlay *>(node);
	if (bu_strcmp(label->labelId.getValue().getString(), labelId) == 0)
	    return i;
    }
    return -1;
}

static int
find_line_layer_overlay_child(SoGroup *group, const char *overlayId)
{
    if (!group || !overlayId)
	return -1;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(SoBRLLineLayerOverlay::getClassTypeId()))
	    continue;
	SoBRLLineLayerOverlay *overlay = static_cast<SoBRLLineLayerOverlay *>(node);
	if (bu_strcmp(overlay->overlayId.getValue().getString(), overlayId) == 0)
	    return i;
    }
    return -1;
}

int
BObolViewController::replaceEditPreview(const char *previewId,
	const char *identity,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	uint32_t sourceRevision,
	uint32_t inputsRevision)
{
    return this->replaceEditPreviewWithIntent(previewId, identity, NULL, NULL,
	    points, commands, count, sourceRevision, inputsRevision);
}

int
BObolViewController::replaceEditPreviewWithIntent(const char *previewId,
	const char *identity,
	const char *editIntentId,
	const char *editIntentRole,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	uint32_t sourceRevision,
	uint32_t inputsRevision)
{
    if (!previewId || !previewId[0] || !points || !commands || count <= 0)
	return -1;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_edit_preview_child(group, previewId);
    SoBRLEditPreview *preview = NULL;
    if (childIndex >= 0)
	preview = static_cast<SoBRLEditPreview *>(group->getChild(childIndex));
    else
	preview = new SoBRLEditPreview;

    if (sourceRevision == 0)
	sourceRevision = preview->sourceRevision.getValue() + 1;
    if (inputsRevision == 0)
	inputsRevision = preview->inputsRevision.getValue() + 1;

    preview->previewId = previewId;
    preview->setEditIntent(editIntentId ? editIntentId : "",
			   editIntentRole ? editIntentRole : "preview");
    preview->sourceRevision = sourceRevision;
    preview->inputsRevision = inputsRevision;
    SoBRLVListShape *shape = preview->setLineSet(
				 (identity && identity[0]) ? identity : previewId,
				 points, commands, count);
    if (!shape)
	return -1;

    if (childIndex < 0)
	group->addChild(preview);

    this->requestRender("edit-preview");
    return 1;
}

int
BObolViewController::removeEditPreview(const char *previewId)
{
    if (!previewId || !previewId[0])
	return 0;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_edit_preview_child(group, previewId);
    if (childIndex < 0)
	return 0;

    group->removeChild(childIndex);
    this->requestRender("edit-preview");
    return 1;
}

int
BObolViewController::replaceLineLayerOverlay(const char *overlayId,
	const struct bg_line_layer_builder *builder,
	uint32_t sourceId,
	SbBool selectable)
{
    if (!overlayId || !overlayId[0])
	return -1;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_line_layer_overlay_child(group, overlayId);
    if (!builder) {
	if (childIndex >= 0) {
	    group->removeChild(childIndex);
	    this->requestRender("line-layer-overlay");
	}
	return 0;
    }

    SoBRLLineLayerOverlay *overlay = new SoBRLLineLayerOverlay;
    overlay->overlayId = overlayId;
    overlay->sourceId = sourceId;
    overlay->selectable = selectable;
    const int realized = overlay->rebuildGeometry(builder);

    if (childIndex >= 0)
	group->replaceChild(childIndex, overlay);
    else
	group->addChild(overlay);

    this->requestRender("line-layer-overlay");
    return realized;
}

int
BObolViewController::replaceHUDLabelOverlay(const char *labelId,
	const char *text,
	const SbVec2f &position,
	const SbColor &color,
	float fontSize,
	uint32_t sourceId)
{
    if (!labelId || !labelId[0])
	return -1;
    if (!text || !text[0])
	return this->removeHUDLabelOverlay(labelId);

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_hud_label_overlay_child(group, labelId);
    SoBRLHUDLabelOverlay *label = childIndex >= 0 ?
				  static_cast<SoBRLHUDLabelOverlay *>(group->getChild(childIndex)) :
				  new SoBRLHUDLabelOverlay;

    label->labelId = labelId;
    label->sourceId = sourceId;
    label->text = text;
    label->position = position;
    label->color = color;
    label->fontSize = fontSize;
    label->visible = TRUE;
    label->rebuildGeometry();

    if (childIndex < 0)
	group->addChild(label);

    this->requestRender("hud-label-overlay");
    return 1;
}

int
BObolViewController::removeHUDLabelOverlay(const char *labelId)
{
    if (!labelId || !labelId[0])
	return 0;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_hud_label_overlay_child(group, labelId);
    if (childIndex < 0)
	return 0;

    group->removeChild(childIndex);
    this->requestRender("hud-label-overlay");
    return 1;
}

SoGroup *
BObolViewController::findGroup(const char *groupPath) const
{
    return this->d->sceneController.findGroup(groupPath);
}

SoGroup *
BObolViewController::ensureGroup(const char *groupPath)
{
    const uint64_t revision = this->d->sceneController.getStructuralRevision();
    SoGroup *group = this->d->sceneController.ensureGroup(groupPath);
    if (group && this->d->sceneController.getStructuralRevision() != revision)
	this->requestRender("scene-group");
    return group;
}

int
BObolViewController::setGroupDrawIntent(const char *groupPath,
	const char *intentPath,
	int drawMode,
	int fallbackDrawMode,
	SbBool overlayIntent,
	uint32_t revalidationRevision)
{
    const int changed = this->d->sceneController.setGroupDrawIntent(groupPath,
			intentPath, drawMode, fallbackDrawMode, overlayIntent,
			revalidationRevision);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::setGroupDisplayState(const char *groupPath,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    const int changed = this->d->sceneController.setGroupDisplayState(
			    groupPath, visible, selected, highlighted, lineStyle,
			    lineWidth, transparency, colorOverride, color,
			    materialColorValid, materialColor, materialRevision);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::renameGroup(const char *groupPath,
				   const char *newLeafName)
{
    const int changed =
	this->d->sceneController.renameGroup(groupPath, newLeafName);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::appendChildToGroup(const char *groupPath,
	SoNode *child)
{
    const int changed =
	this->d->sceneController.appendChildToGroup(groupPath, child);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::removeChildFromGroup(const char *groupPath,
	SoNode *child)
{
    const int changed =
	this->d->sceneController.removeChildFromGroup(groupPath, child);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::eraseGroupSubpath(const char *parentGroupPath,
	const char *subpath)
{
    const int changed =
	this->d->sceneController.eraseGroupSubpath(parentGroupPath, subpath);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BObolViewController::removeGroup(const char *groupPath)
{
    const int removed = this->d->sceneController.removeGroup(groupPath);
    if (removed > 0)
	this->requestRender("scene-group");
    return removed;
}

int
BObolViewController::clearGroup(const char *groupPath)
{
    const int removed = this->d->sceneController.clearGroup(groupPath);
    if (removed > 0)
	this->requestRender("scene-group");
    return removed;
}

int
BObolViewController::getGroupChildCount(const char *groupPath) const
{
    return this->d->sceneController.getGroupChildCount(groupPath);
}

int
BObolViewController::getGroupDescendantGroupCount(
    const char *groupPath) const
{
    return this->d->sceneController.getGroupDescendantGroupCount(groupPath);
}

int
BObolViewController::getGroupDatabaseSourceCount(
    const char *groupPath) const
{
    return this->d->sceneController.getGroupDatabaseSourceCount(groupPath);
}

SoNode *
BObolViewController::findShape(const char *shapePath) const
{
    return this->d->sceneController.findShape(shapePath);
}

SoGroup *
BObolViewController::findShapeParent(const char *shapePath) const
{
    return this->d->sceneController.findShapeParent(shapePath);
}

int
BObolViewController::moveShapeToGroup(const char *shapePath,
					const char *groupPath)
{
    const int changed =
	this->d->sceneController.moveShapeToGroup(shapePath, groupPath);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BObolViewController::removeShape(const char *shapePath)
{
    const int removed = this->d->sceneController.removeShape(shapePath);
    if (removed > 0)
	this->requestRender("scene-shape");
    return removed;
}

int
BObolViewController::setShapeDrawState(const char *shapePath,
	int drawMode,
	SbBool databaseIntent,
	SbBool overlayIntent,
	SbBool hudIntent)
{
    const int changed = this->d->sceneController.setShapeDrawState(shapePath,
			drawMode, databaseIntent, overlayIntent, hudIntent);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BObolViewController::setShapeDisplayState(const char *shapePath,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    const int changed = this->d->sceneController.setShapeDisplayState(
			    shapePath, visible, selected, highlighted, lineStyle, lineWidth,
			    transparency, colorOverride, color, materialColorValid,
			    materialColor, materialRevision);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BObolViewController::setShapePlacementState(const char *shapePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize)
{
    const int changed = this->d->sceneController.setShapePlacementState(
			    shapePath, drawMatrixValid, drawMatrix, drawCenterValid,
			    drawCenter, drawSizeValid, drawSize);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BObolViewController::setShapeSourceState(const char *shapePath,
	const char *ownerSourcePath,
	uint32_t ownerSourceRevision,
	uint32_t ownerInputsRevision,
	uint32_t ownerViewRevision,
	uint32_t ownerRealizedRevision,
	uint32_t ownerRealizedSourceRevision,
	uint32_t ownerRealizedInputsRevision,
	uint32_t ownerRealizedViewRevision,
	int ownerRealizationStatus,
	const char *ownerRealizationDiagnostic,
	const char *ownerRealizationIdentity,
	SbBool ownerSourceStale,
	uint32_t ownerStaleReason)
{
    const int changed = this->d->sceneController.setShapeSourceState(
			    shapePath, ownerSourcePath, ownerSourceRevision,
			    ownerInputsRevision, ownerViewRevision, ownerRealizedRevision,
			    ownerRealizedSourceRevision, ownerRealizedInputsRevision,
			    ownerRealizedViewRevision, ownerRealizationStatus,
			    ownerRealizationDiagnostic, ownerRealizationIdentity,
			    ownerSourceStale, ownerStaleReason);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BObolViewController::replaceDatabaseSource(const char *sourcePath,
	struct db_i *dbip,
	int drawMode,
	uint32_t sourceRevision)
{
    int changed = this->d->sceneController.replaceDatabaseSource(sourcePath,
		  dbip, drawMode, sourceRevision);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::replaceDatabaseSourceInstance(
    const char *sourceInstanceKey,
    const char *sourcePath,
    struct db_i *dbip,
    int drawMode,
    uint32_t sourceRevision)
{
    int changed = this->d->sceneController.replaceDatabaseSourceInstance(
		      sourceInstanceKey, sourcePath, dbip, drawMode, sourceRevision);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceState(const char *sourcePath,
	SbBool sourceRevisionValid,
	uint32_t sourceRevision,
	uint32_t inputsRevision,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    SoBRLDatabaseSource *source =
	this->d->sceneController.findDatabaseSource(sourcePath);
    const uint32_t previousSourceRevision = source ?
	source->sourceRevision.getValue() : 0;
    const uint32_t previousInputsRevision = source ?
	source->inputsRevision.getValue() : 0;
    const SbBool previousStale = source ? source->stale.getValue() : FALSE;
    const uint32_t previousStaleReason = source ?
	source->staleReason.getValue() : 0;
    const int changed = this->d->sceneController.setDatabaseSourceState(
			    sourcePath, sourceRevisionValid, sourceRevision, inputsRevision,
			    visible, selected, highlighted, lineStyle, lineWidth, transparency,
			    colorOverride, color, materialColorValid, materialColor,
			    materialRevision);
    if (changed > 0) {
	source = this->d->sceneController.findDatabaseSource(sourcePath);
	if (!source ||
	    source->sourceRevision.getValue() != previousSourceRevision ||
	    source->inputsRevision.getValue() != previousInputsRevision ||
	    source->stale.getValue() != previousStale ||
	    source->staleReason.getValue() != previousStaleReason)
	    this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceInstanceState(
    const char *sourceInstanceKey,
    SbBool sourceRevisionValid,
    uint32_t sourceRevision,
    uint32_t inputsRevision,
    SbBool visible,
    SbBool selected,
    SbBool highlighted,
    int lineStyle,
    int lineWidth,
    float transparency,
    SbBool colorOverride,
    const SbColor &color,
    SbBool materialColorValid,
    const SbColor &materialColor,
    uint32_t materialRevision)
{
    SoBRLDatabaseSource *source =
	this->d->sceneController.findDatabaseSourceInstance(sourceInstanceKey);
    const uint32_t previousSourceRevision = source ?
	source->sourceRevision.getValue() : 0;
    const uint32_t previousInputsRevision = source ?
	source->inputsRevision.getValue() : 0;
    const SbBool previousStale = source ? source->stale.getValue() : FALSE;
    const uint32_t previousStaleReason = source ?
	source->staleReason.getValue() : 0;
    const int changed = this->d->sceneController.setDatabaseSourceInstanceState(
			    sourceInstanceKey, sourceRevisionValid, sourceRevision,
			    inputsRevision, visible, selected, highlighted, lineStyle, lineWidth,
			    transparency, colorOverride, color, materialColorValid,
			    materialColor, materialRevision);
    if (changed > 0) {
	source = this->d->sceneController.findDatabaseSourceInstance(
		 sourceInstanceKey);
	if (!source ||
	    source->sourceRevision.getValue() != previousSourceRevision ||
	    source->inputsRevision.getValue() != previousInputsRevision ||
	    source->stale.getValue() != previousStale ||
	    source->staleReason.getValue() != previousStaleReason)
	    this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceDisplayPatch(const char *sourcePath,
	const BObolDatabaseSourceDisplayPatch &patch)
{
    const int changed = this->d->sceneController.setDatabaseSourceDisplayPatch(
			    sourcePath, patch);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceInstanceDisplayPatch(
    const char *sourceInstanceKey,
    const BObolDatabaseSourceDisplayPatch &patch)
{
    const int changed =
	this->d->sceneController.setDatabaseSourceInstanceDisplayPatch(
	    sourceInstanceKey, patch);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceDisplayName(const char *sourcePath,
	const char *displayName)
{
    const int changed = this->d->sceneController.setDatabaseSourceDisplayName(
			    sourcePath, displayName);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceInstanceDisplayName(
    const char *sourceInstanceKey,
    const char *displayName)
{
    const int changed =
	this->d->sceneController.setDatabaseSourceInstanceDisplayName(
	    sourceInstanceKey, displayName);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceBoundsState(const char *sourcePath,
	SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax,
	SbBool boundsExact)
{
    const int changed = this->d->sceneController.setDatabaseSourceBoundsState(
			    sourcePath, boundsValid, boundsMin, boundsMax,
			    boundsExact);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceInstanceBoundsState(
    const char *sourceInstanceKey,
    SbBool boundsValid,
    const SbVec3f &boundsMin,
    const SbVec3f &boundsMax,
    SbBool boundsExact)
{
    const int changed =
	this->d->sceneController.setDatabaseSourceInstanceBoundsState(
	    sourceInstanceKey, boundsValid, boundsMin, boundsMax,
	    boundsExact);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceMaterialPolicy(
    const char *sourcePath,
    int materialPolicy)
{
    const int changed = this->d->sceneController.setDatabaseSourceMaterialPolicy(
			    sourcePath, materialPolicy);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::setDatabaseSourceInstanceMaterialPolicy(
    const char *sourceInstanceKey,
    int materialPolicy)
{
    const int changed =
	this->d->sceneController.setDatabaseSourceInstanceMaterialPolicy(
	    sourceInstanceKey, materialPolicy);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::markDatabaseSourceStale(const char *sourcePath,
	uint32_t staleReason)
{
    const int changed = this->d->sceneController.markDatabaseSourceStale(
			    sourcePath, staleReason);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::markDatabaseSourceInstanceStale(
    const char *sourceInstanceKey,
    uint32_t staleReason)
{
    const int changed = this->d->sceneController.markDatabaseSourceInstanceStale(
			    sourceInstanceKey, staleReason);
    if (changed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BObolViewController::removeDatabaseSource(const char *sourcePath)
{
    int removed = this->d->sceneController.removeDatabaseSource(sourcePath);
    if (removed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return removed;
}

int
BObolViewController::removeDatabaseSourceInstance(
    const char *sourceInstanceKey)
{
    int removed = this->d->sceneController.removeDatabaseSourceInstance(
		      sourceInstanceKey);
    if (removed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return removed;
}

int
BObolViewController::moveDatabaseSourceToGroup(const char *sourcePath,
	const char *groupPath)
{
    int moved = this->d->sceneController.moveDatabaseSourceToGroup(sourcePath,
		groupPath);
    if (moved > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return moved;
}

int
BObolViewController::moveDatabaseSourceInstanceToGroup(
    const char *sourceInstanceKey,
    const char *groupPath)
{
    int moved = this->d->sceneController.moveDatabaseSourceInstanceToGroup(
		    sourceInstanceKey, groupPath);
    if (moved > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return moved;
}

int
BObolViewController::clearDatabaseSources(void)
{
    int removed = this->d->sceneController.clearDatabaseSources();
    if (removed > 0) {
	this->invalidateDatabaseSourceLodState();
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return removed;
}

SoBRLDatabaseSource *
BObolViewController::getDatabaseSource(int index) const
{
    return this->d->sceneController.getDatabaseSource(index);
}

int
BObolViewController::getDatabaseSourceCount(void) const
{
    return this->d->sceneController.getDatabaseSourceCount();
}

std::vector<SoBRLDatabaseSource *>
BObolViewController::getRenderDatabaseSources(void) const
{
    return controller_render_database_sources(this);
}

SoBRLDatabaseSource *
BObolViewController::findDatabaseSourceInstance(
    const char *sourceInstanceKey) const
{
    return this->d->sceneController.findDatabaseSourceInstance(
	       sourceInstanceKey);
}

SbBool
BObolViewController::getDatabaseSourceSummary(int index,
	BObolDatabaseSourceSummary &summary) const
{
    return this->d->sceneController.getDatabaseSourceSummary(index, summary);
}

