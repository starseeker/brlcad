/*         D A T A B A S E _ S O U R C E _ S U M M A R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file database_source_summary.cpp
 *
 * Read-only inspection of realized database-source presentation.  Keeping
 * traversal/report construction separate from realization and compact
 * mutation makes the source ownership boundary visible without adding an
 * accessor layer to the occurrence hot path.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BMaterialObject.h"
#include "BObol/BMeshShape.h"
#include "BObol/BVListShape.h"
#include "database_source_private.h"

#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoMatrixTransform.h>

static int
count_material_objects_in_node(SoNode *node)
{
    if (!node)
	return 0;

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return 1;

    int ret = 0;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    ret += count_material_objects_in_node(group->getChild(i));
    }

    return ret;
}

int
SoBRLDatabaseSource::getRealizedMaterialObjectCount(void) const
{
    int ret = 0;
    for (int i = 0; i < this->getNumChildren(); i++)
	ret += count_material_objects_in_node(this->getChild(i));
    return ret;
}

static void
realized_material_summary_owner(const SoBRLDatabaseSource *source,
				int sourceIndex, BObolRealizedMaterialSummary &summary);

static void
realized_material_summary(const SoBRLMaterialObject *object,
			  BObolRealizedMaterialSummary &summary)
{
    summary = BObolRealizedMaterialSummary();
    if (!object)
	return;

    summary.valid = TRUE;
    summary.sourcePath = object->sourcePath.getValue();
    summary.sourceName = object->sourceName.getValue();
    summary.sourceType = object->sourceType.getValue();
    summary.sourceId = object->sourceId.getValue();
    summary.materialName = object->materialName.getValue();
    summary.parentName = object->parentName.getValue();
    summary.materialSource = object->materialSource.getValue();
    summary.propertyCount = object->getPropertyCount();
}

int
SoBRLDatabaseSource::getRealizedMaterialSummaryCount(void) const
{
    return this->getRealizedMaterialObjectCount();
}

SbBool
SoBRLDatabaseSource::getRealizedMaterialSummary(int index,
	BObolRealizedMaterialSummary &summary) const
{
    summary = BObolRealizedMaterialSummary();
    if (index < 0)
	return FALSE;

    SoBRLMaterialObject *object = this->getRealizedMaterialObject(index);
    if (!object)
	return FALSE;

    realized_material_summary(object, summary);
    realized_material_summary_owner(this, -1, summary);
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getRealizedMaterialProperty(int materialIndex,
	int propertyIndex, SbString &groupOut, SbString &nameOut,
	SbString &valueOut) const
{
    if (materialIndex < 0 || propertyIndex < 0)
	return FALSE;

    SoBRLMaterialObject *object =
	this->getRealizedMaterialObject(materialIndex);
    if (!object)
	return FALSE;

    return object->getProperty(propertyIndex, groupOut, nameOut, valueOut);
}

static void
realized_summary_owner_instance_from_node(const SoNode *node,
	SbString &ownerSourceInstanceKey)
{
    if (!node || ownerSourceInstanceKey.getLength() > 0)
	return;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	ownerSourceInstanceKey = source_effective_instance_key(source);
	return;
    }
    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	ownerSourceInstanceKey =
	    static_cast<const SoBRLVListShape *>(node)->
	    ownerSourceInstanceKey.getValue();
	return;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	ownerSourceInstanceKey =
	    static_cast<const SoBRLMeshShape *>(node)->
	    ownerSourceInstanceKey.getValue();
	return;
    }
}

static void
realized_tree_summary_fill(const SoNode *node, int depth, SbBool hasParent,
			   int ownerSourceIndex, const SbString &ownerSourcePath,
			   BObolSceneTreeSummary &summary)
{
    summary = BObolSceneTreeSummary();
    if (!node)
	return;

    summary.valid = TRUE;
    summary.hasParent = hasParent;
    summary.drawTreeDepth = depth;
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    realized_summary_owner_instance_from_node(node,
	    summary.ownerSourceInstanceKey);
    summary.isGroup =
	node->isOfType(SoGroup::getClassTypeId()) ? TRUE : FALSE;
    summary.isDatabaseSource =
	node->isOfType(SoBRLDatabaseSource::getClassTypeId()) ? TRUE : FALSE;
    summary.isMaterialObject =
	node->isOfType(SoBRLMaterialObject::getClassTypeId()) ? TRUE : FALSE;
    summary.isShape =
	(node->isOfType(SoBRLVListShape::getClassTypeId()) ||
	 node->isOfType(SoBRLMeshShape::getClassTypeId())) ? TRUE : FALSE;

    if (summary.isGroup)
	summary.childCount = static_cast<const SoGroup *>(node)->getNumChildren();

    if (summary.isDatabaseSource) {
	const SoBRLDatabaseSource *sourceNode =
	    static_cast<const SoBRLDatabaseSource *>(node);
	summary.nodeKind = BObolSceneTreeSummary::NODE_DATABASE_SOURCE;
	summary.path = sourceNode->path.getValue();
	summary.displayName = sourceNode->displayName.getValue();
	if (summary.ownerSourcePath.getLength() == 0)
	    summary.ownerSourcePath = sourceNode->path.getValue();
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		source_effective_instance_key(sourceNode);
	return;
    }

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	const SoBRLVListShape *shape =
	    static_cast<const SoBRLVListShape *>(node);
	summary.nodeKind = BObolSceneTreeSummary::NODE_VLIST_SHAPE;
	summary.path = shape->sourcePath.getValue();
	summary.sourceName = shape->sourceName.getValue();
	summary.sourceType = shape->sourceType.getValue();
	summary.sourceId = shape->sourceId.getValue();
	summary.displayName = shape->displayName.getValue();
	summary.geometryName = shape->geometryName.getValue();
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	const SoBRLMeshShape *shape =
	    static_cast<const SoBRLMeshShape *>(node);
	summary.nodeKind = BObolSceneTreeSummary::NODE_MESH_SHAPE;
	summary.path = shape->sourcePath.getValue();
	summary.sourceName = shape->sourceName.getValue();
	summary.sourceType = shape->sourceType.getValue();
	summary.sourceId = shape->sourceId.getValue();
	summary.displayName = shape->displayName.getValue();
	summary.geometryName = shape->geometryName.getValue();
	return;
    }

    if (summary.isMaterialObject) {
	const SoBRLMaterialObject *object =
	    static_cast<const SoBRLMaterialObject *>(node);
	summary.nodeKind = BObolSceneTreeSummary::NODE_MATERIAL_OBJECT;
	summary.path = object->sourcePath.getValue();
	summary.sourceName = object->sourceName.getValue();
	summary.sourceType = object->sourceType.getValue();
	summary.sourceId = object->sourceId.getValue();
	summary.displayName = object->materialName.getValue();
	return;
    }

    summary.nodeKind = summary.isGroup ?
		       BObolSceneTreeSummary::NODE_GROUP :
		       BObolSceneTreeSummary::NODE_OTHER;
}

static int
realized_tree_summary_node_count(const SoNode *node)
{
    if (!node)
	return 0;
    if (node_is_source_placement_transform(node))
	return 0;

    int count = 1;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    count += realized_tree_summary_node_count(group->getChild(i));
    }
    return count;
}

static SbBool
find_realized_tree_summary_in_node(const SoNode *node, int &index, int depth,
				   SbBool hasParent, int ownerSourceIndex,
				   const SbString &ownerSourcePath, BObolSceneTreeSummary &summary)
{
    if (!node)
	return FALSE;
    if (node_is_source_placement_transform(node))
	return FALSE;

    if (index == 0) {
	realized_tree_summary_fill(node, depth, hasParent, ownerSourceIndex,
				   ownerSourcePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_tree_summary_in_node(group->getChild(i),
						   index, depth + 1, TRUE, ownerSourceIndex,
						   ownerSourcePath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedTreeSummaryCount(void) const
{
    return realized_tree_summary_node_count(this);
}

SbBool
SoBRLDatabaseSource::getRealizedTreeSummary(int index,
	BObolSceneTreeSummary &summary) const
{
    summary = BObolSceneTreeSummary();
    if (index < 0)
	return FALSE;

    const SbString ownerPath = this->path.getValue();
    const SbBool ret = find_realized_tree_summary_in_node(this, index, 0,
		       FALSE, -1, ownerPath, summary);
    if (ret && summary.ownerSourceInstanceKey.getLength() == 0)
	summary.ownerSourceInstanceKey = source_effective_instance_key(this);
    return ret;
}

static void
realized_display_summary_fill_common(BObolSceneDisplaySummary &summary,
				     int nodeKind, int ownerSourceIndex, const SbString &ownerSourcePath,
				     const SbString &nodePath)
{
    summary = BObolSceneDisplaySummary();
    summary.valid = TRUE;
    summary.nodeKind = nodeKind;
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    summary.path = nodePath;
}

template <typename ShapeT>
static void
realized_display_summary_fill_shape(const ShapeT *shape,
				    int nodeKind, int ownerSourceIndex, const SbString &ownerSourcePath,
				    BObolSceneDisplaySummary &summary)
{
    realized_display_summary_fill_common(summary, nodeKind, ownerSourceIndex,
					 ownerSourcePath, shape ? shape->sourcePath.getValue() : SbString(""));
    if (!shape)
	return;

    summary.ownerSourceInstanceKey = shape->ownerSourceInstanceKey.getValue();
    summary.hasDrawIntent = shape->sourcePath.getValue().getLength() > 0;
    summary.intentPath = shape->sourcePath.getValue();
    summary.intentDrawMode = shape->drawMode.getValue();
    summary.visible = shape->visible.getValue();
    summary.highlighted = shape->highlighted.getValue();
    summary.lineStyle = shape->lineStyle.getValue();
    summary.lineWidth = shape->lineWidth.getValue();
    summary.transparency = shape->transparency.getValue();
    summary.drawMode = shape->drawMode.getValue();
    summary.materialValid = TRUE;
    summary.materialRevision = shape->materialRevision.getValue();
    if (shape->materialColorValid.getValue())
	summary.materialColor = shape->materialColor.getValue();
    else if (shape->colorOverride.getValue())
	summary.materialColor = shape->color.getValue();
    summary.drawMatrixValid = shape->drawMatrixValid.getValue();
    summary.drawMatrix = shape->drawMatrix.getValue();
    summary.drawCenterValid = shape->drawCenterValid.getValue();
    summary.drawCenter = shape->drawCenter.getValue();
    summary.drawSizeValid = shape->drawSizeValid.getValue();
    summary.drawSize = shape->drawSize.getValue();
}

static void
realized_display_summary_fill(const SoNode *node, int ownerSourceIndex,
			      const SbString &ownerSourcePath, BObolSceneDisplaySummary &summary)
{
    summary = BObolSceneDisplaySummary();
    if (!node)
	return;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	realized_display_summary_fill_common(summary,
					     BObolSceneTreeSummary::NODE_DATABASE_SOURCE,
					     ownerSourceIndex, ownerSourcePath, source->path.getValue());
	summary.ownerSourceInstanceKey = source_effective_instance_key(source);
	summary.isDatabaseSource = TRUE;
	summary.hasDrawIntent = source->path.getValue().getLength() > 0;
	summary.intentPath = source->path.getValue();
	summary.intentDrawMode = source->drawMode.getValue();
	summary.visible = source->visible.getValue();
	summary.selected = source->selected.getValue();
	summary.highlighted = source->highlighted.getValue();
	summary.lineStyle = source->lineStyle.getValue();
	summary.lineWidth = source->lineWidth.getValue();
	summary.transparency = source->transparency.getValue();
	summary.drawMode = source->drawMode.getValue();
	summary.materialValid = TRUE;
	summary.materialRevision = source->materialRevision.getValue();
	if (source->materialColorValid.getValue())
	    summary.materialColor = source->materialColor.getValue();
	else if (source->colorOverride.getValue())
	    summary.materialColor = source->color.getValue();
	return;
    }

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	realized_display_summary_fill_shape(
	    static_cast<const SoBRLVListShape *>(node),
	    BObolSceneTreeSummary::NODE_VLIST_SHAPE,
	    ownerSourceIndex, ownerSourcePath, summary);
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	realized_display_summary_fill_shape(
	    static_cast<const SoBRLMeshShape *>(node),
	    BObolSceneTreeSummary::NODE_MESH_SHAPE,
	    ownerSourceIndex, ownerSourcePath, summary);
	return;
    }

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId())) {
	const SoBRLMaterialObject *object =
	    static_cast<const SoBRLMaterialObject *>(node);
	realized_display_summary_fill_common(summary,
					     BObolSceneTreeSummary::NODE_MATERIAL_OBJECT,
					     ownerSourceIndex, ownerSourcePath,
					     object->sourcePath.getValue());
	return;
    }

    const int nodeKind = node->isOfType(SoGroup::getClassTypeId()) ?
			 BObolSceneTreeSummary::NODE_GROUP :
			 BObolSceneTreeSummary::NODE_OTHER;
    realized_display_summary_fill_common(summary, nodeKind, ownerSourceIndex,
					 ownerSourcePath, "");
}

static SbBool
find_realized_display_summary_in_node(const SoNode *node, int &index,
				      int ownerSourceIndex, const SbString &ownerSourcePath,
				      BObolSceneDisplaySummary &summary)
{
    if (!node)
	return FALSE;
    if (node_is_source_placement_transform(node))
	return FALSE;

    if (index == 0) {
	realized_display_summary_fill(node, ownerSourceIndex,
				      ownerSourcePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_display_summary_in_node(group->getChild(i),
		index, ownerSourceIndex, ownerSourcePath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedDisplaySummaryCount(void) const
{
    return this->getRealizedTreeSummaryCount();
}

SbBool
SoBRLDatabaseSource::getRealizedDisplaySummary(int index,
	BObolSceneDisplaySummary &summary) const
{
    summary = BObolSceneDisplaySummary();
    if (index < 0)
	return FALSE;

    const SbString ownerPath = this->path.getValue();
    const SbBool ret = find_realized_display_summary_in_node(this, index,
		       -1, ownerPath, summary);
    if (ret && summary.ownerSourceInstanceKey.getLength() == 0)
	summary.ownerSourceInstanceKey = source_effective_instance_key(this);
    return ret;
}

static void
scene_material_summary_from_display(const BObolSceneDisplaySummary &display,
				    BObolSceneMaterialSummary &summary)
{
    summary = BObolSceneMaterialSummary();
    if (!display.valid)
	return;

    summary.valid = TRUE;
    summary.nodeKind = display.nodeKind;
    summary.materialValid =
	(display.nodeKind == BObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	 display.nodeKind == BObolSceneTreeSummary::NODE_MESH_SHAPE) ?
	display.materialValid : FALSE;
    summary.materialRevision = display.materialRevision;
    summary.materialColor = display.materialColor;
    summary.ownerSourceIndex = display.ownerSourceIndex;
    summary.ownerSourcePath = display.ownerSourcePath;
    summary.ownerSourceInstanceKey = display.ownerSourceInstanceKey;
    summary.path = display.path;
}

int
SoBRLDatabaseSource::getRealizedSceneMaterialSummaryCount(void) const
{
    return this->getRealizedDisplaySummaryCount();
}

SbBool
SoBRLDatabaseSource::getRealizedSceneMaterialSummary(int index,
	BObolSceneMaterialSummary &summary) const
{
    summary = BObolSceneMaterialSummary();
    BObolSceneDisplaySummary display;
    if (!this->getRealizedDisplaySummary(index, display))
	return FALSE;

    scene_material_summary_from_display(display, summary);
    return TRUE;
}

static int
realized_bounds_node_kind(const SoNode *node)
{
    if (!node)
	return BObolSceneTreeSummary::NODE_UNKNOWN;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return BObolSceneTreeSummary::NODE_DATABASE_SOURCE;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return BObolSceneTreeSummary::NODE_VLIST_SHAPE;
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return BObolSceneTreeSummary::NODE_MESH_SHAPE;
    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return BObolSceneTreeSummary::NODE_MATERIAL_OBJECT;
    if (node->isOfType(SoGroup::getClassTypeId()))
	return BObolSceneTreeSummary::NODE_GROUP;
    return BObolSceneTreeSummary::NODE_OTHER;
}

static SbString
realized_bounds_node_path(const SoNode *node)
{
    if (!node)
	return "";
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return static_cast<const SoBRLDatabaseSource *>(node)->path.getValue();
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->sourcePath.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->sourcePath.getValue();
    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return static_cast<const SoBRLMaterialObject *>(node)->sourcePath.getValue();
    return "";
}

static SbBool
realized_bounds_for_vlist_shape(const SoBRLVListShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    const SoBRLVListShape *geom = shape->getGeometrySource();
    for (int i = 0; i < geom->point.getNum(); i++)
	bounds.extendBy(geom->point[i]);
    return geom->point.getNum() > 0;
}

static SbBool
realized_bounds_for_mesh_shape(const SoBRLMeshShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    for (int i = 0; i < geom->point.getNum(); i++)
	bounds.extendBy(geom->point[i]);
    return geom->point.getNum() > 0;
}

static SbBool
realized_bounds_for_node(const SoNode *node, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!node)
	return FALSE;

    SoGetBoundingBoxAction bboxAction(SbViewportRegion(1, 1));
    bboxAction.apply(const_cast<SoNode *>(node));
    bounds = bboxAction.getBoundingBox();
    if (!bounds.isEmpty())
	return TRUE;

    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return realized_bounds_for_vlist_shape(
		   static_cast<const SoBRLVListShape *>(node), bounds);

    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return realized_bounds_for_mesh_shape(
		   static_cast<const SoBRLMeshShape *>(node), bounds);

    SbBool valid = FALSE;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SbBox3f childBounds;
	    if (realized_bounds_for_node(group->getChild(i),
					 childBounds)) {
		bounds.extendBy(childBounds);
		valid = TRUE;
	    }
	}
    }

    if (valid)
	return TRUE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()) &&
	static_cast<const SoBRLDatabaseSource *>(node)->
	getEffectiveSourceBounds(bounds))
	return TRUE;

    return valid;
}

static void
realized_bounds_summary_fill(const SoNode *node, int ownerSourceIndex,
			     const SbString &ownerSourcePath, BObolSceneBoundsSummary &summary)
{
    summary = BObolSceneBoundsSummary();
    if (!node)
	return;

    summary.valid = TRUE;
    summary.nodeKind = realized_bounds_node_kind(node);
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    realized_summary_owner_instance_from_node(node,
	    summary.ownerSourceInstanceKey);
    summary.path = realized_bounds_node_path(node);
    summary.boundsValid = realized_bounds_for_node(node, summary.bounds);
}

static SbBool
find_realized_bounds_summary_in_node(const SoNode *node, int &index,
				     int ownerSourceIndex, const SbString &ownerSourcePath,
				     BObolSceneBoundsSummary &summary)
{
    if (!node)
	return FALSE;
    if (node_is_source_placement_transform(node))
	return FALSE;

    if (index == 0) {
	realized_bounds_summary_fill(node, ownerSourceIndex, ownerSourcePath,
				     summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_bounds_summary_in_node(group->getChild(i),
		index, ownerSourceIndex, ownerSourcePath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedBoundsSummaryCount(void) const
{
    return this->getRealizedTreeSummaryCount();
}

SbBool
SoBRLDatabaseSource::getRealizedBoundsSummary(int index,
	BObolSceneBoundsSummary &summary) const
{
    summary = BObolSceneBoundsSummary();
    if (index < 0)
	return FALSE;

    const SbString ownerPath = this->path.getValue();
    const SbBool ret = find_realized_bounds_summary_in_node(this, index,
		       -1, ownerPath, summary);
    if (ret && summary.ownerSourceInstanceKey.getLength() == 0)
	summary.ownerSourceInstanceKey = source_effective_instance_key(this);
    return ret;
}

template <typename ShapeT>
static void
realized_shape_summary_common(const ShapeT *shape,
			      BObolRealizedShapeSummary &summary)
{
    summary.path = shape->sourcePath.getValue();
    summary.sourceName = shape->sourceName.getValue();
    summary.sourceType = shape->sourceType.getValue();
    summary.sourceId = shape->sourceId.getValue();
    summary.displayName = shape->displayName.getValue();
    summary.geometryName = shape->geometryName.getValue();
    summary.cacheIdentity = shape->cacheIdentity.getValue();
    summary.sourceIdentity = shape->sourceIdentity.getValue();
    summary.databaseIntent = shape->databaseIntent.getValue();
    summary.overlayIntent = shape->overlayIntent.getValue();
    summary.hudIntent = shape->hudIntent.getValue();
    summary.localSource = shape->localSource.getValue();
    summary.sharedSource = shape->sharedSource.getValue();
    summary.nonDatabaseSource = shape->nonDatabaseSource.getValue();
    summary.drawMode = shape->drawMode.getValue();
    summary.recordRole = shape->recordRole.getValue();
    summary.geometryKind = shape->geometryKind.getValue();
    summary.regionId = shape->regionId.getValue();
    summary.airCode = shape->airCode.getValue();
    summary.materialId = shape->materialId.getValue();
    summary.los = shape->los.getValue();
    summary.materialColorValid = shape->materialColorValid.getValue();
    summary.materialColor = shape->materialColor.getValue();
    summary.materialShader = shape->materialShader.getValue();
    summary.materialRevision = shape->materialRevision.getValue();
    summary.visible = shape->visible.getValue();
    summary.selectable = shape->selectable.getValue();
    summary.selected = shape->selected.getValue();
    summary.highlighted = shape->highlighted.getValue();
    summary.ghosted = shape->ghosted.getValue();
    summary.hiddenLine = shape->hiddenLine.getValue();
    summary.editEmphasis = shape->editEmphasis.getValue();
    summary.lineStyle = shape->lineStyle.getValue();
    summary.lineWidth = shape->lineWidth.getValue();
    summary.transparency = shape->transparency.getValue();
    summary.editIntentId = shape->editIntentId.getValue();
    summary.editIntentRole = shape->editIntentRole.getValue();
    summary.lodPolicy = shape->lodPolicy.getValue();
    summary.colorOverride = shape->colorOverride.getValue();
    summary.color = shape->color.getValue();
}

static void
realized_shape_summary_owner(const SoBRLDatabaseSource *source,
			     int sourceIndex, BObolRealizedShapeSummary &summary)
{
    if (!source)
	return;

    summary.ownerSourceIndex = sourceIndex;
    summary.ownerSourcePath = source->path.getValue();
    summary.ownerSourceInstanceKey = source_effective_instance_key(source);
    summary.ownerDrawMode = source->drawMode.getValue();
    summary.ownerSourceRevision = source->sourceRevision.getValue();
    summary.ownerInputsRevision = source->inputsRevision.getValue();
    summary.ownerViewRevision = source->viewRevision.getValue();
    summary.ownerRealizedRevision = source->realizedRevision.getValue();
    summary.ownerRealizedSourceRevision =
	source->realizedSourceRevision.getValue();
    summary.ownerRealizedInputsRevision =
	source->realizedInputsRevision.getValue();
    summary.ownerRealizedViewRevision =
	source->realizedViewRevision.getValue();
    summary.ownerRealizationStatus = source->realizationStatus.getValue();
    summary.ownerRealizationDiagnostic =
	source->realizationDiagnostic.getValue();
    summary.ownerRealizationIdentity =
	source->realizationIdentity.getValue();
    summary.ownerSourceStale = source->stale.getValue();
    summary.ownerStaleReason = source->staleReason.getValue();
}

static void
realized_material_summary_owner(const SoBRLDatabaseSource *source,
				int sourceIndex, BObolRealizedMaterialSummary &summary)
{
    if (!source)
	return;

    summary.ownerSourceIndex = sourceIndex;
    summary.ownerSourcePath = source->path.getValue();
    summary.ownerSourceInstanceKey = source_effective_instance_key(source);
    summary.ownerDrawMode = source->drawMode.getValue();
    summary.ownerSourceRevision = source->sourceRevision.getValue();
    summary.ownerInputsRevision = source->inputsRevision.getValue();
    summary.ownerViewRevision = source->viewRevision.getValue();
    summary.ownerRealizedRevision = source->realizedRevision.getValue();
    summary.ownerRealizedSourceRevision =
	source->realizedSourceRevision.getValue();
    summary.ownerRealizedInputsRevision =
	source->realizedInputsRevision.getValue();
    summary.ownerRealizedViewRevision =
	source->realizedViewRevision.getValue();
    summary.ownerRealizationStatus = source->realizationStatus.getValue();
    summary.ownerRealizationDiagnostic =
	source->realizationDiagnostic.getValue();
    summary.ownerRealizationIdentity =
	source->realizationIdentity.getValue();
    summary.ownerSourceStale = source->stale.getValue();
    summary.ownerStaleReason = source->staleReason.getValue();
}

static void
realized_shape_summary_bounds(const SoMFVec3f &points,
			      BObolRealizedShapeSummary &summary)
{
    summary.bounds.makeEmpty();
    summary.boundsValid = FALSE;
    for (int i = 0; i < points.getNum(); i++) {
	summary.bounds.extendBy(points[i]);
	summary.boundsValid = TRUE;
    }
}

void
realized_vlist_shape_summary(const SoBRLVListShape *shape,
			     BObolRealizedShapeSummary &summary)
{
    summary = BObolRealizedShapeSummary();
    if (!shape)
	return;

    summary.valid = TRUE;
    summary.shapeKind = BObolRealizedShapeSummary::SHAPE_VLIST;
    realized_shape_summary_common(shape, summary);
    const SoBRLVListShape *geom = shape->getGeometrySource();
    summary.pointCount = geom->point.getNum();
    summary.commandCount = geom->command.getNum();
    summary.segmentCount = shape->getSegmentCount();
    summary.pointPrimitiveCount = shape->getPointPrimitiveCount();
    realized_shape_summary_bounds(geom->point, summary);
}

void
realized_mesh_shape_summary(const SoBRLMeshShape *shape,
			    BObolRealizedShapeSummary &summary)
{
    summary = BObolRealizedShapeSummary();
    if (!shape)
	return;

    summary.valid = TRUE;
    summary.shapeKind = BObolRealizedShapeSummary::SHAPE_MESH;
    realized_shape_summary_common(shape, summary);
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    summary.pointCount = geom->point.getNum();
    summary.indexCount = geom->coordIndex.getNum();
    summary.triangleCount = shape->getTriangleCount();
    realized_shape_summary_bounds(geom->point, summary);
    summary.lodAvailable = shape->lodAvailable.getValue();
    summary.lodActiveCut = shape->lodActiveCut.getValue();
    summary.lodFaceCount = shape->lodFaceCount.getValue();
    summary.lodPointCount = shape->lodPointCount.getValue();
    summary.lodOriginalPointCount = shape->lodOriginalPointCount.getValue();
    summary.lodNormalCount = shape->lodNormalCount.getValue();
    summary.lodHasSnappedPoints = shape->lodHasSnappedPoints.getValue();
    summary.lodHasNormals = shape->lodHasNormals.getValue();
    summary.lodBoundsMin = shape->lodBoundsMin.getValue();
    summary.lodBoundsMax = shape->lodBoundsMax.getValue();
}

static SbBool
find_realized_shape_summary_in_node(SoNode *node, int &index,
				    BObolRealizedShapeSummary &summary)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	if (index == 0) {
	    realized_vlist_shape_summary(
		static_cast<SoBRLVListShape *>(node), summary);
	    return TRUE;
	}
	index--;
	return FALSE;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	if (index == 0) {
	    realized_mesh_shape_summary(
		static_cast<SoBRLMeshShape *>(node), summary);
	    return TRUE;
	}
	index--;
	return FALSE;
    }

    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_shape_summary_in_node(group->getChild(i),
						    index, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedShapeSummaryCount(void) const
{
    int count = this->d->compactIndexActive && this->d->compactIndex ?
	static_cast<int>(this->d->compactIndex->entries.size()) : 0;
    count += this->getRealizedShapeCount();
    count += this->getRealizedMeshCount();
    return count;
}

SbBool
SoBRLDatabaseSource::getRealizedShapeSummary(int index,
	BObolRealizedShapeSummary &summary) const
{
    summary = BObolRealizedShapeSummary();
    if (index < 0)
	return FALSE;

    if (this->d->compactIndexActive && this->d->compactIndex) {
	const size_t compactCount = this->d->compactIndex->entries.size();
	if (static_cast<size_t>(index) < compactCount) {
	    summary = this->d->compactIndex->entries[
		static_cast<size_t>(index)].shapeSummary;
	    realized_shape_summary_owner(this, -1, summary);
	    return summary.valid;
	}
	index -= static_cast<int>(compactCount);
    }

    int remaining = index;
    for (int i = 0; i < this->getNumChildren(); i++) {
	if (find_realized_shape_summary_in_node(this->getChild(i),
						remaining, summary)) {
	    realized_shape_summary_owner(this, -1, summary);
	    return TRUE;
	}
    }

    return FALSE;
}

SbBool
SoBRLDatabaseSource::getSummary(BObolDatabaseSourceSummary &summary) const
{
    summary = BObolDatabaseSourceSummary();
    summary.valid = TRUE;
    summary.path = this->path.getValue();
    summary.instanceKey = source_effective_instance_key(this);
    summary.parentInstanceKey = this->parentInstanceKey.getValue();
    summary.occurrenceIndex = this->occurrenceIndex.getValue();
    summary.booleanOperation = this->booleanOperation.getValue();
    summary.displayName = this->displayName.getValue();
    summary.representationKey =
	source_effective_representation_key(this);
    summary.representationMode = this->representationMode.getValue();
    summary.drawMode = this->drawMode.getValue();
    summary.sourceRevision = this->sourceRevision.getValue();
    summary.inputsRevision = this->inputsRevision.getValue();
    summary.viewRevision = this->viewRevision.getValue();
    summary.realizedRevision = this->realizedRevision.getValue();
    summary.realizedSourceRevision = this->realizedSourceRevision.getValue();
    summary.realizedInputsRevision = this->realizedInputsRevision.getValue();
    summary.realizedViewRevision = this->realizedViewRevision.getValue();
    summary.realizationStatus = this->realizationStatus.getValue();
    summary.realizationDiagnostic = this->realizationDiagnostic.getValue();
    summary.realizationIdentity = this->realizationIdentity.getValue();
    summary.realizationRoleFlags = this->realizationRoleFlags.getValue();
    summary.realizationViewDependent =
	this->realizationViewDependent.getValue();
    summary.realizationCsgLodEnabled =
	this->realizationCsgLodEnabled.getValue();
    summary.realizationMeshLodEnabled =
	this->realizationMeshLodEnabled.getValue();
    summary.realizationViewScale = this->realizationViewScale.getValue();
    summary.realizationLodScale = this->realizationLodScale.getValue();
    summary.realizationViewWidth = this->realizationViewWidth.getValue();
    summary.realizationViewHeight = this->realizationViewHeight.getValue();
    summary.realizationBotThreshold =
	this->realizationBotThreshold.getValue();
    summary.realizationCurveScale = this->realizationCurveScale.getValue();
    summary.realizationPointScale = this->realizationPointScale.getValue();
    summary.visible = this->visible.getValue();
    summary.selected = this->selected.getValue();
    summary.highlighted = this->highlighted.getValue();
    summary.lineStyle = this->lineStyle.getValue();
    summary.lineWidth = this->lineWidth.getValue();
    summary.transparency = this->transparency.getValue();
    summary.materialColorValid = this->materialColorValid.getValue();
    summary.materialColor = this->materialColor.getValue();
    summary.materialRevision = this->materialRevision.getValue();
    summary.materialPolicy = this->materialPolicy.getValue();
    summary.databaseMetadataValid = this->databaseMetadataValid.getValue();
    summary.databaseRegionId = this->databaseRegionId.getValue();
    summary.databaseAirCode = this->databaseAirCode.getValue();
    summary.databaseMaterialId = this->databaseMaterialId.getValue();
    summary.databaseLos = this->databaseLos.getValue();
    summary.databaseMaterialColorValid =
	this->databaseMaterialColorValid.getValue();
    summary.databaseMaterialColor = this->databaseMaterialColor.getValue();
    summary.databaseMaterialShader = this->databaseMaterialShader.getValue();
    summary.colorOverride = this->colorOverride.getValue();
    summary.color = this->color.getValue();
    summary.selectedColor = this->selectedColor.getValue();
    summary.highlightedColor = this->highlightedColor.getValue();
    summary.ghostedColor = this->ghostedColor.getValue();
    summary.drawMatrixValid = this->drawMatrixValid.getValue();
    summary.drawMatrix = this->drawMatrix.getValue();
    summary.drawCenterValid = this->drawCenterValid.getValue();
    summary.drawCenter = this->drawCenter.getValue();
    summary.drawSizeValid = this->drawSizeValid.getValue();
    summary.drawSize = this->drawSize.getValue();
    summary.sourceBoundsValid = this->getEffectiveSourceBounds(
				    summary.sourceBounds);
    summary.sourceBoundsExact = this->hasExactSourceBounds();
    summary.hasViewDependentCsgGeometry = this->d->compactIndex &&
	this->d->compactIndex->viewDependentCsgGeometryCount > 0 ? TRUE : FALSE;
    summary.stale = this->stale.getValue();
    summary.staleReason = this->staleReason.getValue();
    summary.realizedShapeCount = this->getRealizedShapeCount();
    summary.realizedMeshCount = this->getRealizedMeshCount();
    summary.realizedMaterialObjectCount =
	this->getRealizedMaterialObjectCount();
    return TRUE;
}
