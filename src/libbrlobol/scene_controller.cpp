/*                S C E N E _ C O N T R O L L E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/line_layer_overlay.h"
#include "brlobol/material_object.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/realize_action.h"
#include "brlobol/scene_controller.h"
#include "brlobol/scene_group.h"
#include "brlobol/vlist_shape.h"
#include "performance_private.h"

#include <Inventor/SbName.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>

#include <algorithm>
#include <string.h>
#include <string>
#include <unordered_map>
#include <vector>

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
scene_path_equal(const char *stored, const char *path)
{
    if (!stored || !path)
	return 0;
    if (strcmp(stored, path) == 0)
	return 1;
    return strcmp(skip_leading_slash(stored), skip_leading_slash(path)) == 0;
}

static int
database_source_path_equal(const SoBRLDatabaseSource *source,
			   const char *path)
{
    if (!source || !path)
	return 0;
    const char *sourcePath = source->path.getValue().getString();
    return scene_path_equal(sourcePath, path);
}

static SbString
database_source_effective_instance_key(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";
    const SbString key = source->instanceKey.getValue();
    if (key.getLength() > 0)
	return key;
    return source->path.getValue();
}

static SbString
database_source_effective_representation_key(
    const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";
    const SbString key = source->representationKey.getValue();
    if (key.getLength() > 0)
	return key;
    return database_source_effective_instance_key(source);
}

static int
database_source_instance_key_equal(const SoBRLDatabaseSource *source,
				   const char *instanceKey)
{
    if (!source || !instanceKey)
	return 0;
    const SbString key = database_source_effective_instance_key(source);
    const char *stored = key.getString();
    return scene_path_equal(stored, instanceKey);
}

static void
index_put(std::unordered_map<std::string, SoBRLDatabaseSource *> &index,
	  const char *key,
	  SoBRLDatabaseSource *source)
{
    if (!key || !key[0] || !source)
	return;

    const char *normalized = skip_leading_slash(key);
    if (normalized && normalized[0])
	index[std::string(normalized)] = source;
}

static void
index_erase(std::unordered_map<std::string, SoBRLDatabaseSource *> &index,
	    const char *key)
{
    if (!key || !key[0])
	return;

    const char *normalized = skip_leading_slash(key);
    if (normalized && normalized[0])
	index.erase(std::string(normalized));
}

static void
parent_index_put(std::unordered_map<std::string, SoGroup *> &index,
		 const char *key,
		 SoGroup *parent)
{
    if (!key || !key[0] || !parent)
	return;

    const char *normalized = skip_leading_slash(key);
    if (normalized && normalized[0])
	index[std::string(normalized)] = parent;
}

static void
parent_index_erase(std::unordered_map<std::string, SoGroup *> &index,
		   const char *key)
{
    if (!key || !key[0])
	return;

    const char *normalized = skip_leading_slash(key);
    if (normalized && normalized[0])
	index.erase(std::string(normalized));
}

static void
group_index_put(std::unordered_map<std::string, SoGroup *> &index,
		const char *key,
		SoGroup *group)
{
    if (!key || !key[0] || !group)
	return;

    const char *normalized = skip_leading_slash(key);
    if (normalized && normalized[0])
	index[std::string(normalized)] = group;
}

static SbString
scene_group_index_path(const SoGroup *group)
{
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return "";

    const SoBRLSceneGroup *sceneGroup =
	static_cast<const SoBRLSceneGroup *>(group);
    return sceneGroup->groupPath.getValue();
}

static void
source_order_put(std::vector<SoBRLDatabaseSource *> &order,
		 std::unordered_map<SoBRLDatabaseSource *, size_t> &orderIndex,
		 SoBRLDatabaseSource *source)
{
    if (!source || orderIndex.find(source) != orderIndex.end())
	return;

    orderIndex[source] = order.size();
    order.push_back(source);
}

static void
source_order_erase(std::vector<SoBRLDatabaseSource *> &order,
		   std::unordered_map<SoBRLDatabaseSource *, size_t> &orderIndex,
		   SoBRLDatabaseSource *source)
{
    if (!source)
	return;

    auto it = orderIndex.find(source);
    if (it == orderIndex.end())
	return;

    const size_t pos = it->second;
    orderIndex.erase(it);
    if (pos >= order.size())
	return;

    order.erase(order.begin() + static_cast<std::ptrdiff_t>(pos));
    for (size_t i = pos; i < order.size(); i++)
	orderIndex[order[i]] = i;
}

static SoBRLDatabaseSource *
find_database_source_instance_recursive(SoGroup *group,
					const char *sourceInstanceKey,
					SoGroup **parentOut = NULL, int *childIndexOut = NULL)
{
    if (parentOut)
	*parentOut = NULL;
    if (childIndexOut)
	*childIndexOut = -1;
    if (!group || !sourceInstanceKey)
	return NULL;

    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    SoBRLDatabaseSource *source =
		static_cast<SoBRLDatabaseSource *>(node);
	    if (database_source_instance_key_equal(source,
						   sourceInstanceKey)) {
		if (parentOut)
		    *parentOut = group;
		if (childIndexOut)
		    *childIndexOut = i;
		return source;
	    }
	}
	if (node->isOfType(SoGroup::getClassTypeId())) {
	    SoBRLDatabaseSource *found =
		find_database_source_instance_recursive(
		    static_cast<SoGroup *>(node), sourceInstanceKey,
		    parentOut, childIndexOut);
	    if (found)
		return found;
	}
    }

    return NULL;
}

static void
index_database_sources_recursive(
    SoGroup *group,
    std::unordered_map<std::string, SoGroup *> &groupPathIndex,
    std::unordered_map<std::string, SoBRLDatabaseSource *> &pathIndex,
    std::unordered_map<std::string, SoBRLDatabaseSource *> &instanceIndex,
    std::unordered_map<std::string, SoGroup *> &instanceParentIndex,
    std::vector<SoBRLDatabaseSource *> &sourceOrder,
    std::unordered_map<SoBRLDatabaseSource *, size_t> &sourceOrderIndex)
{
    if (!group)
	return;

    const SbString groupPath = scene_group_index_path(group);
    group_index_put(groupPathIndex, groupPath.getString(), group);

    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    SoBRLDatabaseSource *source =
		static_cast<SoBRLDatabaseSource *>(node);
	    const SbString path = source->path.getValue();
	    const SbString instanceKey =
		database_source_effective_instance_key(source);
	    index_put(pathIndex, path.getString(), source);
	    index_put(instanceIndex, instanceKey.getString(), source);
	    parent_index_put(instanceParentIndex, instanceKey.getString(), group);
	    source_order_put(sourceOrder, sourceOrderIndex, source);
	}
	if (node->isOfType(SoGroup::getClassTypeId()))
	    index_database_sources_recursive(static_cast<SoGroup *>(node),
		groupPathIndex, pathIndex, instanceIndex, instanceParentIndex,
		sourceOrder, sourceOrderIndex);
    }
}

static int
count_database_sources_recursive(const SoGroup *group)
{
    if (!group)
	return 0;

    int count = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    count++;
	}
	if (node->isOfType(SoGroup::getClassTypeId()))
	    count += count_database_sources_recursive(
			 static_cast<const SoGroup *>(node));
    }
    return count;
}

static int
count_scene_groups_recursive(const SoGroup *group)
{
    if (!group)
	return 0;

    int count = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	    continue;
	if (node->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    count++;
	if (node->isOfType(SoGroup::getClassTypeId()))
	    count += count_scene_groups_recursive(
			 static_cast<const SoGroup *>(node));
    }
    return count;
}

static int
clear_database_sources_recursive(SoGroup *group)
{
    if (!group)
	return 0;

    int removed = 0;
    for (int i = group->getNumChildren() - 1; i >= 0; i--) {
	SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    SoBRLDatabaseSource *source =
		static_cast<SoBRLDatabaseSource *>(node);
	    if (source->auxiliarySource.getValue())
		continue;
	    group->removeChild(i);
	    removed++;
	    continue;
	}
	if (node->isOfType(SoGroup::getClassTypeId()))
	    removed += clear_database_sources_recursive(
			   static_cast<SoGroup *>(node));
    }
    return removed;
}

static void
scene_group_path_components(const char *groupPath,
			    std::vector<std::string> &components)
{
    components.clear();
    if (!groupPath)
	return;

    const char *cursor = groupPath;
    while (*cursor == '/')
	cursor++;

    while (*cursor) {
	const char *start = cursor;
	while (*cursor && *cursor != '/')
	    cursor++;
	if (cursor > start)
	    components.push_back(std::string(start,
					     static_cast<size_t>(cursor - start)));
	while (*cursor == '/')
	    cursor++;
    }
}

static SoGroup *
scene_root_group(SoNode *node)
{
    if (!node || !node->isOfType(SoGroup::getClassTypeId()))
	return NULL;
    return static_cast<SoGroup *>(node);
}

static const SoGroup *
scene_root_group_const(const SoNode *node)
{
    if (!node || !node->isOfType(SoGroup::getClassTypeId()))
	return NULL;
    return static_cast<const SoGroup *>(node);
}

static int
scene_group_node_name_equal(const SoNode *node, const char *leafName)
{
    if (!node || !leafName)
	return 0;
    const SbName nodeName = node->getName();
    if (strcmp(nodeName.getString(), leafName) == 0)
	return 1;

    if (node->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	const SoBRLSceneGroup *group =
	    static_cast<const SoBRLSceneGroup *>(node);
	const char *groupPath = group->groupPath.getValue().getString();
	const char *groupLeaf = strrchr(groupPath, '/');
	groupLeaf = groupLeaf ? groupLeaf + 1 : groupPath;
	if (groupLeaf && strcmp(groupLeaf, leafName) == 0)
	    return 1;
    }

    return 0;
}

static SoGroup *
scene_group_find_child(SoGroup *parent, const char *leafName)
{
    if (!parent || !leafName || !leafName[0])
	return NULL;

    for (int i = 0; i < parent->getNumChildren(); i++) {
	SoNode *child = parent->getChild(i);
	if (!child || !child->isOfType(SoBRLSceneGroup::getClassTypeId()) ||
	    !scene_group_node_name_equal(child, leafName))
	    continue;
	return static_cast<SoGroup *>(child);
    }

    return NULL;
}

static int
scene_group_find_child_index(SoGroup *parent, const SoNode *child)
{
    if (!parent || !child)
	return -1;
    for (int i = 0; i < parent->getNumChildren(); i++) {
	if (parent->getChild(i) == child)
	    return i;
    }
    return -1;
}

static const SoGroup *
scene_group_find_child_const(const SoGroup *parent, const char *leafName)
{
    if (!parent || !leafName || !leafName[0])
	return NULL;

    for (int i = 0; i < parent->getNumChildren(); i++) {
	const SoNode *child = parent->getChild(i);
	if (!child || !child->isOfType(SoBRLSceneGroup::getClassTypeId()) ||
	    !scene_group_node_name_equal(child, leafName))
	    continue;
	return static_cast<const SoGroup *>(child);
    }

    return NULL;
}

static const SoGroup *
scene_group_find_path_const(const SoNode *sceneRoot, const char *groupPath)
{
    if (!groupPath)
	return NULL;

    const SoGroup *current = scene_root_group_const(sceneRoot);
    if (!current)
	return NULL;

    std::vector<std::string> components;
    scene_group_path_components(groupPath, components);
    for (size_t i = 0; i < components.size(); i++) {
	current = scene_group_find_child_const(current,
					       components[i].c_str());
	if (!current)
	    return NULL;
    }

    return current;
}

static SbString
scene_child_summary_path(const SbString &parentPath, const SoNode *child)
{
    if (!child)
	return "";

    const SbName childName = child->getName();
    if (childName.getLength() == 0)
	return "";

    if (parentPath.getLength() == 0 ||
	strcmp(parentPath.getString(), "/") == 0)
	return SbString(childName.getString());

    SbString childPath = parentPath;
    childPath += "/";
    childPath += childName.getString();
    return childPath;
}

static SbString
scene_group_append_path(const SbString &parentPath, const char *leafName)
{
    if (!leafName || !leafName[0])
	return parentPath;

    if (parentPath.getLength() == 0 ||
	strcmp(parentPath.getString(), "/") == 0)
	return SbString(leafName);

    SbString childPath = parentPath;
    childPath += "/";
    childPath += leafName;
    return childPath;
}

static SbString
scene_group_component_path(const std::vector<std::string> &components,
			   size_t count)
{
    SbString path("");
    const size_t limit = count < components.size() ?
			 count : components.size();
    for (size_t i = 0; i < limit; i++)
	path = scene_group_append_path(path, components[i].c_str());
    return path;
}

static SbString
scene_group_summary_path(const SoNode *node, const SbString &fallbackPath)
{
    if (node && node->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	const SoBRLSceneGroup *group =
	    static_cast<const SoBRLSceneGroup *>(node);
	if (group->groupPath.getValue().getLength() > 0)
	    return group->groupPath.getValue();
    }

    if (node && node->isOfType(SoBRLLineLayerOverlay::getClassTypeId())) {
	const SoBRLLineLayerOverlay *overlay =
	    static_cast<const SoBRLLineLayerOverlay *>(node);
	if (overlay->overlayId.getValue().getLength() > 0)
	    return overlay->overlayId.getValue();
    }

    return fallbackPath;
}

static void
scene_group_update_path_recursive(SoGroup *group, const SbString &groupPath)
{
    if (!group)
	return;

    if (group->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	SoBRLSceneGroup *sceneGroup =
	    static_cast<SoBRLSceneGroup *>(group);
	sceneGroup->groupPath = groupPath;
    }

    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *child = group->getChild(i);
	if (!child || !child->isOfType(SoGroup::getClassTypeId()))
	    continue;

	const SbString childPath = scene_child_summary_path(groupPath, child);
	scene_group_update_path_recursive(static_cast<SoGroup *>(child),
					  childPath);
    }
}

static int
scene_group_float_different(float a, float b)
{
    const float delta = a - b;
    return delta < -0.000001f || delta > 0.000001f;
}

static int
scene_group_color_equal(const SbColor &a, const SbColor &b)
{
    return !scene_group_float_different(a[0], b[0]) &&
	   !scene_group_float_different(a[1], b[1]) &&
	   !scene_group_float_different(a[2], b[2]);
}

static int
scene_group_vec3f_equal(const SbVec3f &a, const SbVec3f &b)
{
    return !scene_group_float_different(a[0], b[0]) &&
	   !scene_group_float_different(a[1], b[1]) &&
	   !scene_group_float_different(a[2], b[2]);
}

static int
scene_group_set_string(SoSFString &field, const SbString &value)
{
    if (strcmp(field.getValue().getString(), value.getString()) == 0)
	return 0;
    field = value;
    return 1;
}

static int
scene_node_is_shape(const SoNode *node)
{
    return node &&
	   (node->isOfType(SoBRLVListShape::getClassTypeId()) ||
	    node->isOfType(SoBRLMeshShape::getClassTypeId()));
}

static SbString
scene_shape_node_path(const SoNode *node)
{
    if (!node)
	return "";
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->sourcePath.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->sourcePath.getValue();
    return "";
}

static int
scene_shape_path_equal(const SoNode *node, const char *shapePath)
{
    if (!node || !shapePath || !shapePath[0])
	return 0;

    const SbString nodePath = scene_shape_node_path(node);
    if (nodePath.getLength() == 0)
	return 0;

    const char *stored = nodePath.getString();
    if (strcmp(stored, shapePath) == 0)
	return 1;
    return strcmp(skip_leading_slash(stored),
		  skip_leading_slash(shapePath)) == 0;
}

static SoNode *
scene_shape_find_in_group(SoGroup *group, const char *shapePath,
			  SoGroup **parentOut)
{
    if (!group || !shapePath)
	return NULL;

    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *child = group->getChild(i);
	if (!child)
	    continue;

	if (scene_node_is_shape(child) &&
	    scene_shape_path_equal(child, shapePath)) {
	    if (parentOut)
		*parentOut = group;
	    return child;
	}

	if (!child->isOfType(SoGroup::getClassTypeId()) ||
	    child->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	    continue;

	SoNode *found = scene_shape_find_in_group(
			    static_cast<SoGroup *>(child), shapePath, parentOut);
	if (found)
	    return found;
    }

    return NULL;
}

static SoNode *
scene_shape_find_path(SoNode *sceneRoot, const char *shapePath,
		      SoGroup **parentOut)
{
    if (parentOut)
	*parentOut = NULL;
    if (!shapePath || !shapePath[0])
	return NULL;

    SoGroup *rootGroup = scene_root_group(sceneRoot);
    if (!rootGroup)
	return NULL;
    return scene_shape_find_in_group(rootGroup, shapePath, parentOut);
}

template <typename ShapeT>
static int
scene_shape_set_draw_state_typed(ShapeT *shape,
				 int drawMode,
				 SbBool databaseIntent,
				 SbBool overlayIntent,
				 SbBool hudIntent)
{
    if (!shape)
	return 0;

    int changed = 0;
    if (shape->drawMode.getValue() != drawMode) {
	shape->drawMode = drawMode;
	changed = 1;
    }
    if (shape->databaseIntent.getValue() != databaseIntent) {
	shape->databaseIntent = databaseIntent;
	changed = 1;
    }
    if (shape->overlayIntent.getValue() != overlayIntent) {
	shape->overlayIntent = overlayIntent;
	changed = 1;
    }
    if (shape->hudIntent.getValue() != hudIntent) {
	shape->hudIntent = hudIntent;
	changed = 1;
    }
    return changed;
}

static int
scene_shape_set_draw_state(SoNode *node,
			   int drawMode,
			   SbBool databaseIntent,
			   SbBool overlayIntent,
			   SbBool hudIntent)
{
    if (!node)
	return 0;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return scene_shape_set_draw_state_typed(
		   static_cast<SoBRLVListShape *>(node), drawMode,
		   databaseIntent, overlayIntent, hudIntent);
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return scene_shape_set_draw_state_typed(
		   static_cast<SoBRLMeshShape *>(node), drawMode,
		   databaseIntent, overlayIntent, hudIntent);
    return 0;
}

template <typename ShapeT>
static int
scene_shape_set_display_state_typed(ShapeT *shape,
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
    if (!shape)
	return 0;

    int changed = 0;
    if (shape->visible.getValue() != visible) {
	shape->visible = visible;
	changed = 1;
    }
    if (shape->selected.getValue() != selected) {
	shape->selected = selected;
	changed = 1;
    }
    if (shape->highlighted.getValue() != highlighted) {
	shape->highlighted = highlighted;
	changed = 1;
    }
    if (shape->lineStyle.getValue() != lineStyle) {
	shape->lineStyle = lineStyle;
	changed = 1;
    }
    if (shape->lineWidth.getValue() != lineWidth) {
	shape->lineWidth = lineWidth;
	changed = 1;
    }
    if (scene_group_float_different(shape->transparency.getValue(),
				    transparency)) {
	shape->transparency = transparency;
	changed = 1;
    }
    if (shape->colorOverride.getValue() != colorOverride) {
	shape->colorOverride = colorOverride;
	changed = 1;
    }
    if (!scene_group_color_equal(shape->color.getValue(), color)) {
	shape->color = color;
	changed = 1;
    }
    if (shape->materialColorValid.getValue() != materialColorValid) {
	shape->materialColorValid = materialColorValid;
	changed = 1;
    }
    if (!scene_group_color_equal(shape->materialColor.getValue(),
				 materialColor)) {
	shape->materialColor = materialColor;
	changed = 1;
    }
    if (shape->materialRevision.getValue() != materialRevision) {
	shape->materialRevision = materialRevision;
	changed = 1;
    }
    return changed;
}

static int
scene_shape_set_display_state(SoNode *node,
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
    if (!node)
	return 0;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return scene_shape_set_display_state_typed(
		   static_cast<SoBRLVListShape *>(node), visible, selected,
		   highlighted, lineStyle, lineWidth, transparency,
		   colorOverride, color, materialColorValid, materialColor,
		   materialRevision);
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return scene_shape_set_display_state_typed(
		   static_cast<SoBRLMeshShape *>(node), visible, selected,
		   highlighted, lineStyle, lineWidth, transparency,
		   colorOverride, color, materialColorValid, materialColor,
		   materialRevision);
    return 0;
}

template <typename ShapeT>
static int
scene_shape_set_placement_state_typed(ShapeT *shape,
				      SbBool drawMatrixValid,
				      const SbMatrix &drawMatrix,
				      SbBool drawCenterValid,
				      const SbVec3f &drawCenter,
				      SbBool drawSizeValid,
				      float drawSize)
{
    if (!shape)
	return 0;

    int changed = 0;
    if (shape->drawMatrixValid.getValue() != drawMatrixValid) {
	shape->drawMatrixValid = drawMatrixValid;
	changed = 1;
    }
    if (!shape->drawMatrix.getValue().equals(drawMatrix, 0.000001f)) {
	shape->drawMatrix = drawMatrix;
	changed = 1;
    }
    if (shape->drawCenterValid.getValue() != drawCenterValid) {
	shape->drawCenterValid = drawCenterValid;
	changed = 1;
    }
    if (!scene_group_vec3f_equal(shape->drawCenter.getValue(), drawCenter)) {
	shape->drawCenter = drawCenter;
	changed = 1;
    }
    if (shape->drawSizeValid.getValue() != drawSizeValid) {
	shape->drawSizeValid = drawSizeValid;
	changed = 1;
    }
    if (scene_group_float_different(shape->drawSize.getValue(), drawSize)) {
	shape->drawSize = drawSize;
	changed = 1;
    }
    return changed;
}

static int
scene_shape_set_placement_state(SoNode *node,
				SbBool drawMatrixValid,
				const SbMatrix &drawMatrix,
				SbBool drawCenterValid,
				const SbVec3f &drawCenter,
				SbBool drawSizeValid,
				float drawSize)
{
    if (!node)
	return 0;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return scene_shape_set_placement_state_typed(
		   static_cast<SoBRLVListShape *>(node), drawMatrixValid,
		   drawMatrix, drawCenterValid, drawCenter, drawSizeValid,
		   drawSize);
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return scene_shape_set_placement_state_typed(
		   static_cast<SoBRLMeshShape *>(node), drawMatrixValid,
		   drawMatrix, drawCenterValid, drawCenter, drawSizeValid,
		   drawSize);
    return 0;
}

template <typename ShapeT>
static int
scene_shape_set_source_state_typed(ShapeT *shape,
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
    if (!shape)
	return 0;

    const SbString sourcePath(ownerSourcePath ? ownerSourcePath : "");
    const SbString diagnostic(ownerRealizationDiagnostic ?
			      ownerRealizationDiagnostic : "");
    const SbString identity(ownerRealizationIdentity ?
			    ownerRealizationIdentity : "");
    int changed = 0;
    changed |= scene_group_set_string(shape->ownerSourcePath, sourcePath);
    if (shape->ownerSourceRevision.getValue() != ownerSourceRevision) {
	shape->ownerSourceRevision = ownerSourceRevision;
	changed = 1;
    }
    if (shape->ownerInputsRevision.getValue() != ownerInputsRevision) {
	shape->ownerInputsRevision = ownerInputsRevision;
	changed = 1;
    }
    if (shape->ownerViewRevision.getValue() != ownerViewRevision) {
	shape->ownerViewRevision = ownerViewRevision;
	changed = 1;
    }
    if (shape->ownerRealizedRevision.getValue() != ownerRealizedRevision) {
	shape->ownerRealizedRevision = ownerRealizedRevision;
	changed = 1;
    }
    if (shape->ownerRealizedSourceRevision.getValue() !=
	ownerRealizedSourceRevision) {
	shape->ownerRealizedSourceRevision = ownerRealizedSourceRevision;
	changed = 1;
    }
    if (shape->ownerRealizedInputsRevision.getValue() !=
	ownerRealizedInputsRevision) {
	shape->ownerRealizedInputsRevision = ownerRealizedInputsRevision;
	changed = 1;
    }
    if (shape->ownerRealizedViewRevision.getValue() !=
	ownerRealizedViewRevision) {
	shape->ownerRealizedViewRevision = ownerRealizedViewRevision;
	changed = 1;
    }
    if (shape->ownerRealizationStatus.getValue() !=
	ownerRealizationStatus) {
	shape->ownerRealizationStatus = ownerRealizationStatus;
	changed = 1;
    }
    changed |= scene_group_set_string(shape->ownerRealizationDiagnostic,
				      diagnostic);
    changed |= scene_group_set_string(shape->ownerRealizationIdentity,
				      identity);
    if (shape->ownerSourceStale.getValue() != ownerSourceStale) {
	shape->ownerSourceStale = ownerSourceStale;
	changed = 1;
    }
    if (shape->ownerStaleReason.getValue() != ownerStaleReason) {
	shape->ownerStaleReason = ownerStaleReason;
	changed = 1;
    }
    return changed;
}

static int
scene_shape_set_source_state(SoNode *node,
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
    if (!node)
	return 0;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return scene_shape_set_source_state_typed(
		   static_cast<SoBRLVListShape *>(node), ownerSourcePath,
		   ownerSourceRevision, ownerInputsRevision, ownerViewRevision,
		   ownerRealizedRevision, ownerRealizedSourceRevision,
		   ownerRealizedInputsRevision, ownerRealizedViewRevision,
		   ownerRealizationStatus, ownerRealizationDiagnostic,
		   ownerRealizationIdentity, ownerSourceStale,
		   ownerStaleReason);
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return scene_shape_set_source_state_typed(
		   static_cast<SoBRLMeshShape *>(node), ownerSourcePath,
		   ownerSourceRevision, ownerInputsRevision, ownerViewRevision,
		   ownerRealizedRevision, ownerRealizedSourceRevision,
		   ownerRealizedInputsRevision, ownerRealizedViewRevision,
		   ownerRealizationStatus, ownerRealizationDiagnostic,
		   ownerRealizationIdentity, ownerSourceStale,
		   ownerStaleReason);
    return 0;
}

static SbString
scene_shape_owner_source_path(const SoNode *node)
{
    if (!node)
	return "";
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->
	       ownerSourcePath.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->
	       ownerSourcePath.getValue();
    return "";
}

static SbString
scene_shape_owner_source_instance_key(const SoNode *node)
{
    if (!node)
	return "";
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->
	       ownerSourceInstanceKey.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->
	       ownerSourceInstanceKey.getValue();
    return "";
}

static int
scene_group_erase_nested_subpath(SoGroup *parent,
				 const std::vector<std::string> &components)
{
    if (!parent || components.empty())
	return 0;

    SoGroup *current = parent;
    for (size_t i = 0; i + 1 < components.size(); i++) {
	current = scene_group_find_child(current, components[i].c_str());
	if (!current)
	    return 0;
    }

    SoGroup *target = scene_group_find_child(current,
		      components.back().c_str());
    if (!target)
	return 0;

    current->removeChild(target);
    return 1;
}

BRLObolSceneSummary::BRLObolSceneSummary(void) :
    valid(FALSE),
    hasRoot(FALSE),
    rootIsGroup(FALSE),
    rootChildCount(0),
    databaseSourceCount(0),
    nonDatabaseRootChildCount(0),
    structuralRevision(0),
    frameRevision(0),
    lastVisitedSourceCount(0),
    lastRealizedSourceCount(0),
    lastFailedSourceCount(0),
    lastDiagnostics("")
{
}

BRLObolDatabaseSourcePublishState::BRLObolDatabaseSourcePublishState(void) :
    sourceInstanceKey(NULL),
    sourcePath(NULL),
    sourceRepresentationKey(NULL),
    targetGroupPath(NULL),
    database(NULL),
    drawMode(SoBRLDatabaseSource::WIREFRAME),
    representationMode(SoBRLDatabaseSource::REPRESENTATION_DEFAULT),
    sourceRevisionValid(FALSE),
    sourceRevision(0),
    inputsRevision(0),
    visible(TRUE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    colorOverride(FALSE),
    color(1.0f, 1.0f, 1.0f),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0),
    roleFlagsValid(FALSE),
    roleFlags(SoBRLDatabaseSource::REALIZATION_ROLE_NONE),
    viewPolicyValid(FALSE),
    viewDependent(FALSE),
    viewScale(0.0f),
    botThreshold(0),
    curveScale(0.0f),
    pointScale(0.0f),
    placementValid(FALSE),
    drawMatrixValid(FALSE),
    drawMatrix(SbMatrix::identity()),
    drawCenterValid(FALSE),
    drawCenter(0.0f, 0.0f, 0.0f),
    drawSizeValid(FALSE),
    drawSize(0.0f)
{
}

SoBRLSceneController::SoBRLSceneController(void) :
    root(NULL),
    structuralRevision(0),
    frameRevision(0),
    mutationBatchDepth(0),
    mutationBatchStructuralRevisionPending(FALSE),
    mutationBatchFrameRevisionPending(FALSE),
    compactCadRealizationEnabled(FALSE),
    databaseSourceIndexValid(FALSE),
    lastVisitedSourceCount(0),
    lastRealizedSourceCount(0),
    lastFailedSourceCount(0),
    lastDiagnostics("")
{
}

SoBRLSceneController::SoBRLSceneController(SoNode *sceneRoot) :
    root(NULL),
    structuralRevision(0),
    frameRevision(0),
    mutationBatchDepth(0),
    mutationBatchStructuralRevisionPending(FALSE),
    mutationBatchFrameRevisionPending(FALSE),
    compactCadRealizationEnabled(FALSE),
    databaseSourceIndexValid(FALSE),
    lastVisitedSourceCount(0),
    lastRealizedSourceCount(0),
    lastFailedSourceCount(0),
    lastDiagnostics("")
{
    this->setSceneRoot(sceneRoot);
}

SoBRLSceneController::~SoBRLSceneController(void)
{
    this->setSceneRoot(NULL);
}

void
SoBRLSceneController::setSceneRoot(SoNode *sceneRoot)
{
    if (this->root == sceneRoot)
	return;

    if (sceneRoot)
	sceneRoot->ref();
    if (this->root)
	this->root->unref();
    this->root = sceneRoot;
    this->clearDatabaseSourceIndex();
    this->advanceStructuralRevision();
}

SoNode *
SoBRLSceneController::getSceneRoot(void) const
{
    return this->root;
}

rt_view_scene_ref
SoBRLSceneController::getSceneRef(void)
{
    return rt_view_scene_ref_make(this, RT_VIEW_SCENE_BACKEND_OBOL);
}

rt_view_scene_ref
SoBRLSceneController::getSceneRef(void) const
{
    return rt_view_scene_ref_make(const_cast<SoBRLSceneController *>(this),
				  RT_VIEW_SCENE_BACKEND_OBOL);
}

SbBool
SoBRLSceneController::sceneRefIsObol(rt_view_scene_ref ref)
{
    return (!rt_view_scene_ref_is_null(ref) &&
	    rt_view_scene_ref_backend(ref) == RT_VIEW_SCENE_BACKEND_OBOL) ?
	   TRUE : FALSE;
}

SoBRLSceneController *
SoBRLSceneController::fromSceneRef(rt_view_scene_ref ref)
{
    if (!SoBRLSceneController::sceneRefIsObol(ref))
	return NULL;
    return static_cast<SoBRLSceneController *>(
	       rt_view_scene_ref_context(ref));
}

const SoBRLSceneController *
SoBRLSceneController::fromConstSceneRef(rt_view_scene_ref ref)
{
    return SoBRLSceneController::fromSceneRef(ref);
}

uint64_t
SoBRLSceneController::getStructuralRevision(void) const
{
    return this->structuralRevision;
}

uint64_t
SoBRLSceneController::getFrameRevision(void) const
{
    return this->frameRevision;
}

void
SoBRLSceneController::setCompactCadRealizationEnabled(SbBool enabled)
{
    this->compactCadRealizationEnabled = enabled ? TRUE : FALSE;
}

SbBool
SoBRLSceneController::getCompactCadRealizationEnabled(void) const
{
    return this->compactCadRealizationEnabled;
}

SbBool
SoBRLSceneController::getSceneSummary(BRLObolSceneSummary &summary) const
{
    summary = BRLObolSceneSummary();
    summary.valid = TRUE;
    summary.hasRoot = this->root ? TRUE : FALSE;
    summary.structuralRevision = this->structuralRevision;
    summary.frameRevision = this->frameRevision;
    summary.lastVisitedSourceCount = this->lastVisitedSourceCount;
    summary.lastRealizedSourceCount = this->lastRealizedSourceCount;
    summary.lastFailedSourceCount = this->lastFailedSourceCount;
    summary.lastDiagnostics = this->lastDiagnostics;

    if (!this->root)
	return TRUE;

    summary.rootIsGroup =
	this->root->isOfType(SoGroup::getClassTypeId()) ? TRUE : FALSE;
    if (!summary.rootIsGroup)
	return TRUE;

    SoGroup *group = static_cast<SoGroup *>(this->root);
    summary.rootChildCount = group->getNumChildren();
    summary.databaseSourceCount = this->getDatabaseSourceCount();
    int rootDatabaseSourceCount = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *child = group->getChild(i);
	if (child && child->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	    rootDatabaseSourceCount++;
    }
    summary.nonDatabaseRootChildCount =
	summary.rootChildCount - rootDatabaseSourceCount;
    return TRUE;
}

void
SoBRLSceneController::clearDatabaseSourceIndex(void) const
{
    this->groupPathIndex.clear();
    this->databaseSourcePathIndex.clear();
    this->databaseSourceInstanceIndex.clear();
    this->databaseSourceInstanceParentIndex.clear();
    this->databaseSourceOrder.clear();
    this->databaseSourceOrderIndex.clear();
    this->databaseSourceIndexValid = FALSE;
}

void
SoBRLSceneController::indexSceneGroup(SoGroup *group) const
{
    const SbString groupPath = scene_group_index_path(group);
    group_index_put(this->groupPathIndex, groupPath.getString(), group);
}

void
SoBRLSceneController::indexDatabaseSource(SoBRLDatabaseSource *source,
	SoGroup *parent) const
{
    if (!source)
	return;

    const SbString path = source->path.getValue();
    const SbString instanceKey = database_source_effective_instance_key(source);
    index_put(this->databaseSourcePathIndex, path.getString(), source);
    index_put(this->databaseSourceInstanceIndex, instanceKey.getString(),
	source);
    parent_index_put(this->databaseSourceInstanceParentIndex,
	instanceKey.getString(), parent);
    source_order_put(this->databaseSourceOrder,
	this->databaseSourceOrderIndex, source);
}

void
SoBRLSceneController::unindexDatabaseSource(SoBRLDatabaseSource *source) const
{
    if (!source)
	return;

    const SbString path = source->path.getValue();
    const SbString instanceKey = database_source_effective_instance_key(source);
    index_erase(this->databaseSourcePathIndex, path.getString());
    index_erase(this->databaseSourceInstanceIndex, instanceKey.getString());
    parent_index_erase(this->databaseSourceInstanceParentIndex,
	instanceKey.getString());
    source_order_erase(this->databaseSourceOrder,
	this->databaseSourceOrderIndex, source);
}

void
SoBRLSceneController::rebuildDatabaseSourceIndex(void) const
{
    BRLObolPerformanceTimer timer(BRLOBOL_PERF_SOURCE_INDEX_REBUILD_US);
    if (timer.active())
	brlobol_performance_counter_add(
	    BRLOBOL_PERF_SOURCE_INDEX_REBUILD_CALLS, 1);

    this->groupPathIndex.clear();
    this->databaseSourcePathIndex.clear();
    this->databaseSourceInstanceIndex.clear();
    this->databaseSourceInstanceParentIndex.clear();
    this->databaseSourceOrder.clear();
    this->databaseSourceOrderIndex.clear();
    if (this->root && this->root->isOfType(SoGroup::getClassTypeId())) {
	index_database_sources_recursive(static_cast<SoGroup *>(this->root),
	    this->groupPathIndex,
	    this->databaseSourcePathIndex,
	    this->databaseSourceInstanceIndex,
	    this->databaseSourceInstanceParentIndex,
	    this->databaseSourceOrder,
	    this->databaseSourceOrderIndex);
    }
    this->databaseSourceIndexValid = TRUE;
}

SoGroup *
SoBRLSceneController::findIndexedGroup(const char *groupPath) const
{
    if (!groupPath)
	return NULL;

    std::vector<std::string> components;
    scene_group_path_components(groupPath, components);
    if (components.empty())
	return scene_root_group(this->root);

    if (!this->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    const char *normalized = skip_leading_slash(groupPath);
    auto it = this->groupPathIndex.find(
		  std::string(normalized ? normalized : groupPath));
    if (it != this->groupPathIndex.end())
	return it->second;
    return NULL;
}

SoBRLDatabaseSource *
SoBRLSceneController::findIndexedDatabaseSource(const char *sourcePath) const
{
    if (!sourcePath || !sourcePath[0])
	return NULL;
    if (!this->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    const char *normalized = skip_leading_slash(sourcePath);
    auto it = this->databaseSourcePathIndex.find(
		  std::string(normalized ? normalized : sourcePath));
    if (it != this->databaseSourcePathIndex.end())
	return it->second;
    return NULL;
}

SbString
SoBRLSceneController::databaseSourceInstanceKeyForPath(
    const char *sourcePath) const
{
    SoBRLDatabaseSource *source = this->findIndexedDatabaseSource(sourcePath);
    return database_source_effective_instance_key(source);
}

SoBRLDatabaseSource *
SoBRLSceneController::findIndexedDatabaseSourceInstance(
    const char *sourceInstanceKey) const
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return NULL;
    if (!this->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    const char *normalized = skip_leading_slash(sourceInstanceKey);
    auto it = this->databaseSourceInstanceIndex.find(
		  std::string(normalized ? normalized : sourceInstanceKey));
    if (it != this->databaseSourceInstanceIndex.end())
	return it->second;
    return NULL;
}

SoGroup *
SoBRLSceneController::findIndexedDatabaseSourceInstanceParent(
    const char *sourceInstanceKey) const
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return NULL;
    if (!this->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    const char *normalized = skip_leading_slash(sourceInstanceKey);
    auto it = this->databaseSourceInstanceParentIndex.find(
	      std::string(normalized ? normalized : sourceInstanceKey));
    if (it != this->databaseSourceInstanceParentIndex.end())
	return it->second;
    return NULL;
}

SbBool
SoBRLSceneController::realizePending(void)
{
    this->lastVisitedSourceCount = 0;
    this->lastRealizedSourceCount = 0;
    this->lastFailedSourceCount = 0;
    this->lastDiagnostics = "";

    if (!this->root)
	return FALSE;

    SoBRLRealizeAction action;
    action.setCompactCadRealizationEnabled(this->compactCadRealizationEnabled);
    action.apply(this->root);
    this->lastVisitedSourceCount = action.getVisitedSourceCount();
    this->lastRealizedSourceCount = action.getRealizedSourceCount();
    this->lastFailedSourceCount = action.getFailedSourceCount();
    this->lastDiagnostics = action.getDiagnostics();

    size_t compactedSourceCount = 0;
    if (this->compactCadRealizationEnabled &&
	this->lastFailedSourceCount == 0) {
	if (!this->databaseSourceIndexValid)
	    this->rebuildDatabaseSourceIndex();
	for (size_t i = 0; i < this->databaseSourceOrder.size(); i++) {
	    SoBRLDatabaseSource *source = this->databaseSourceOrder[i];
	    if (!source || source->hasCompactInstanceIndex())
		continue;
	    if (source->compactRealizedGeometry() > 0)
		compactedSourceCount++;
	}
    }

    if (this->lastRealizedSourceCount > 0 ||
	this->lastFailedSourceCount > 0 ||
	compactedSourceCount > 0)
	this->advanceFrameRevision();
    return this->lastFailedSourceCount == 0;
}

void
SoBRLSceneController::beginSceneMutationBatch(size_t expectedDatabaseSources,
	size_t expectedGroups)
{
    if (this->mutationBatchDepth == 0) {
	this->mutationBatchStructuralRevisionPending = FALSE;
	this->mutationBatchFrameRevisionPending = FALSE;
	if (this->root && this->root->isOfType(SoGroup::getClassTypeId()) &&
	    !this->databaseSourceIndexValid)
	    this->rebuildDatabaseSourceIndex();
	if (expectedDatabaseSources > 0) {
	    this->databaseSourcePathIndex.reserve(
		this->databaseSourcePathIndex.size() + expectedDatabaseSources);
	    this->databaseSourceInstanceIndex.reserve(
		this->databaseSourceInstanceIndex.size() +
		expectedDatabaseSources);
	    this->databaseSourceInstanceParentIndex.reserve(
		this->databaseSourceInstanceParentIndex.size() +
		expectedDatabaseSources);
	    this->databaseSourceOrder.reserve(
		this->databaseSourceOrder.size() + expectedDatabaseSources);
	    this->databaseSourceOrderIndex.reserve(
		this->databaseSourceOrderIndex.size() + expectedDatabaseSources);
	}
	if (expectedGroups > 0)
	    this->groupPathIndex.reserve(
		this->groupPathIndex.size() + expectedGroups);
    }

    this->mutationBatchDepth++;
}

void
SoBRLSceneController::endSceneMutationBatch(void)
{
    if (this->mutationBatchDepth <= 0)
	return;

    this->mutationBatchDepth--;
    if (this->mutationBatchDepth > 0)
	return;

    const SbBool structuralChanged =
	this->mutationBatchStructuralRevisionPending;
    const SbBool frameChanged =
	this->mutationBatchFrameRevisionPending;
    this->mutationBatchStructuralRevisionPending = FALSE;
    this->mutationBatchFrameRevisionPending = FALSE;

    if (structuralChanged) {
	this->structuralRevision++;
	if (this->structuralRevision == 0)
	    this->structuralRevision++;
	this->frameRevision++;
	if (this->frameRevision == 0)
	    this->frameRevision++;
    } else if (frameChanged) {
	this->frameRevision++;
	if (this->frameRevision == 0)
	    this->frameRevision++;
    }
}

SoGroup *
SoBRLSceneController::findGroup(const char *groupPath) const
{
    return this->findIndexedGroup(groupPath);
}

SoGroup *
SoBRLSceneController::ensureGroup(const char *groupPath)
{
    if (!groupPath)
	return NULL;

    SoGroup *current = scene_root_group(this->root);
    if (!current)
	return NULL;

    if (!this->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    std::vector<std::string> components;
    scene_group_path_components(groupPath, components);
    SbBool created = FALSE;
    SbString currentPath("");
    for (size_t i = 0; i < components.size(); i++) {
	const SbString childPath =
	    scene_group_append_path(currentPath, components[i].c_str());
	SoGroup *child = this->findIndexedGroup(childPath.getString());
	if (!child) {
	    SoBRLSceneGroup *newGroup = new SoBRLSceneGroup;
	    newGroup->setName(SbName(components[i].c_str()));
	    newGroup->groupPath = childPath;
	    current->addChild(newGroup);
	    child = newGroup;
	    this->indexSceneGroup(child);
	    created = TRUE;
	}
	current = child;
	currentPath = childPath;
    }

    if (created)
	this->advanceStructuralRevision();
    return current;
}

int
SoBRLSceneController::setGroupDrawIntent(const char *groupPath,
	const char *intentPath,
	int drawMode,
	int fallbackDrawMode,
	SbBool overlayIntent,
	uint32_t revalidationRevision)
{
    SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;
    if (!group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return -1;

    SoBRLSceneGroup *sceneGroup = static_cast<SoBRLSceneGroup *>(group);
    SbString retainedPath = sceneGroup->groupPath.getValue();
    if (retainedPath.getLength() == 0)
	retainedPath = skip_leading_slash(groupPath);
    const SbString nextIntentPath =
	(intentPath && intentPath[0]) ? SbString(intentPath) : retainedPath;

    int changed = 0;
    if (!sceneGroup->drawIntentValid.getValue()) {
	sceneGroup->drawIntentValid = TRUE;
	changed = 1;
    }
    changed |= scene_group_set_string(sceneGroup->drawIntentPath,
				      nextIntentPath);
    if (sceneGroup->drawMode.getValue() != drawMode) {
	sceneGroup->drawMode = drawMode;
	changed = 1;
    }
    if (sceneGroup->fallbackDrawMode.getValue() != fallbackDrawMode) {
	sceneGroup->fallbackDrawMode = fallbackDrawMode;
	changed = 1;
    }
    if (sceneGroup->overlayIntent.getValue() != overlayIntent) {
	sceneGroup->overlayIntent = overlayIntent;
	changed = 1;
    }
    if (sceneGroup->revalidationRevision.getValue() !=
	revalidationRevision) {
	sceneGroup->revalidationRevision = revalidationRevision;
	changed = 1;
    }

    if (changed)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setGroupDisplayState(const char *groupPath,
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
    SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;
    if (!group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return -1;

    SoBRLSceneGroup *sceneGroup = static_cast<SoBRLSceneGroup *>(group);
    int changed = 0;
    if (sceneGroup->visible.getValue() != visible) {
	sceneGroup->visible = visible;
	changed = 1;
    }
    if (sceneGroup->selected.getValue() != selected) {
	sceneGroup->selected = selected;
	changed = 1;
    }
    if (sceneGroup->highlighted.getValue() != highlighted) {
	sceneGroup->highlighted = highlighted;
	changed = 1;
    }
    if (sceneGroup->lineStyle.getValue() != lineStyle) {
	sceneGroup->lineStyle = lineStyle;
	changed = 1;
    }
    if (sceneGroup->lineWidth.getValue() != lineWidth) {
	sceneGroup->lineWidth = lineWidth;
	changed = 1;
    }
    if (scene_group_float_different(sceneGroup->transparency.getValue(),
				    transparency)) {
	sceneGroup->transparency = transparency;
	changed = 1;
    }
    if (sceneGroup->colorOverride.getValue() != colorOverride) {
	sceneGroup->colorOverride = colorOverride;
	changed = 1;
    }
    if (!scene_group_color_equal(sceneGroup->color.getValue(), color)) {
	sceneGroup->color = color;
	changed = 1;
    }
    if (sceneGroup->materialColorValid.getValue() !=
	materialColorValid) {
	sceneGroup->materialColorValid = materialColorValid;
	changed = 1;
    }
    if (!scene_group_color_equal(sceneGroup->materialColor.getValue(),
				 materialColor)) {
	sceneGroup->materialColor = materialColor;
	changed = 1;
    }
    if (sceneGroup->materialRevision.getValue() != materialRevision) {
	sceneGroup->materialRevision = materialRevision;
	changed = 1;
    }

    if (changed)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::renameGroup(const char *groupPath,
				  const char *newLeafName)
{
    if (!newLeafName || !newLeafName[0] || strchr(newLeafName, '/'))
	return 0;

    std::vector<std::string> components;
    scene_group_path_components(groupPath, components);
    if (components.empty())
	return 0;

    SoGroup *parent = scene_root_group(this->root);
    if (!parent)
	return -1;

    for (size_t i = 0; i + 1 < components.size(); i++) {
	parent = scene_group_find_child(parent, components[i].c_str());
	if (!parent)
	    return 0;
    }

    SoGroup *target = scene_group_find_child(parent,
		      components.back().c_str());
    if (!target)
	return 0;
    if (scene_group_node_name_equal(target, newLeafName))
	return 0;
    if (scene_group_find_child(parent, newLeafName))
	return 0;

    const SbString parentPath =
	scene_group_component_path(components, components.size() - 1);
    const SbString newGroupPath =
	scene_group_append_path(parentPath, newLeafName);
    target->setName(SbName(newLeafName));
    scene_group_update_path_recursive(target, newGroupPath);
    this->clearDatabaseSourceIndex();
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::appendChildToGroup(const char *groupPath,
	SoNode *child)
{
    if (!child)
	return -1;

    SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;
    if (scene_group_find_child_index(group, child) >= 0)
	return 0;

    group->addChild(child);
    if (child->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	this->indexDatabaseSource(static_cast<SoBRLDatabaseSource *>(child),
	    group);
    else if (child->isOfType(SoGroup::getClassTypeId()))
	this->clearDatabaseSourceIndex();
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::removeChildFromGroup(const char *groupPath,
	SoNode *child)
{
    if (!child)
	return -1;

    SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;

    const int childIndex = scene_group_find_child_index(group, child);
    if (childIndex < 0)
	return 0;

    if (child->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	this->unindexDatabaseSource(static_cast<SoBRLDatabaseSource *>(child));
    else if (child->isOfType(SoGroup::getClassTypeId()))
	this->clearDatabaseSourceIndex();
    group->removeChild(childIndex);
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::eraseGroupSubpath(const char *parentGroupPath,
					const char *subpath)
{
    SoGroup *parent = this->findIndexedGroup(parentGroupPath);
    if (!parent)
	return -1;

    std::vector<std::string> components;
    scene_group_path_components(subpath, components);
    if (components.empty())
	return 0;

    const int erased = scene_group_erase_nested_subpath(parent, components);
    if (erased > 0) {
	this->clearDatabaseSourceIndex();
	this->advanceStructuralRevision();
    }
    return erased;
}

int
SoBRLSceneController::removeGroup(const char *groupPath)
{
    std::vector<std::string> components;
    scene_group_path_components(groupPath, components);
    if (components.empty())
	return 0;

    SoGroup *parent = scene_root_group(this->root);
    if (!parent)
	return -1;

    for (size_t i = 0; i + 1 < components.size(); i++) {
	parent = scene_group_find_child(parent, components[i].c_str());
	if (!parent)
	    return 0;
    }

    SoGroup *target = scene_group_find_child(parent,
		      components.back().c_str());
    if (!target)
	return 0;

    this->clearDatabaseSourceIndex();
    parent->removeChild(target);
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::clearGroup(const char *groupPath)
{
    SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;

    const int removed = group->getNumChildren();
    if (removed <= 0)
	return 0;

    this->clearDatabaseSourceIndex();
    group->removeAllChildren();
    this->advanceStructuralRevision();
    return removed;
}

int
SoBRLSceneController::getGroupChildCount(const char *groupPath) const
{
    const SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;
    return group->getNumChildren();
}

int
SoBRLSceneController::getGroupDescendantGroupCount(
    const char *groupPath) const
{
    const SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;
    return count_scene_groups_recursive(group);
}

int
SoBRLSceneController::getGroupDatabaseSourceCount(
    const char *groupPath) const
{
    const SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;
    return count_database_sources_recursive(group);
}

SoNode *
SoBRLSceneController::findShape(const char *shapePath) const
{
    return scene_shape_find_path(this->root, shapePath, NULL);
}

SoGroup *
SoBRLSceneController::findShapeParent(const char *shapePath) const
{
    SoGroup *parent = NULL;
    (void)scene_shape_find_path(this->root, shapePath, &parent);
    return parent;
}

int
SoBRLSceneController::moveShapeToGroup(const char *shapePath,
				       const char *groupPath)
{
    SoGroup *currentParent = NULL;
    SoNode *shape = scene_shape_find_path(this->root, shapePath,
					  &currentParent);
    if (!shape)
	return 0;

    SoGroup *targetGroup = this->findIndexedGroup(groupPath);
    if (!targetGroup)
	return -1;
    if (targetGroup == currentParent)
	return 0;
    if (scene_group_find_child_index(targetGroup, shape) >= 0)
	return 0;

    const int currentIndex =
	scene_group_find_child_index(currentParent, shape);
    if (currentIndex < 0)
	return -1;

    shape->ref();
    currentParent->removeChild(currentIndex);
    targetGroup->addChild(shape);
    shape->unref();
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::removeShape(const char *shapePath)
{
    SoGroup *parent = NULL;
    SoNode *shape = scene_shape_find_path(this->root, shapePath, &parent);
    if (!shape)
	return 0;

    const int childIndex = scene_group_find_child_index(parent, shape);
    if (childIndex < 0)
	return -1;

    parent->removeChild(childIndex);
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::setShapeDrawState(const char *shapePath,
					int drawMode,
					SbBool databaseIntent,
					SbBool overlayIntent,
					SbBool hudIntent)
{
    SoNode *shape = scene_shape_find_path(this->root, shapePath, NULL);
    if (!shape)
	return -1;

    const int changed = scene_shape_set_draw_state(shape, drawMode,
			databaseIntent, overlayIntent, hudIntent);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setShapeDisplayState(const char *shapePath,
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
    SoNode *shape = scene_shape_find_path(this->root, shapePath, NULL);
    if (!shape)
	return -1;

    const int changed = scene_shape_set_display_state(shape, visible,
			selected, highlighted, lineStyle, lineWidth, transparency,
			colorOverride, color, materialColorValid, materialColor,
			materialRevision);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setShapePlacementState(const char *shapePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize)
{
    SoNode *shape = scene_shape_find_path(this->root, shapePath, NULL);
    if (!shape)
	return -1;

    const int changed = scene_shape_set_placement_state(shape,
			drawMatrixValid, drawMatrix, drawCenterValid, drawCenter,
			drawSizeValid, drawSize);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::publishDatabaseSourceAuxiliaryLineSet(
    const char *sourcePath,
    const char *name,
    const SbVec3f *points,
    const int32_t *commands,
    int count,
    const BRLObolAuxiliaryLineSetDisplayState *displayState)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->publishDatabaseSourceInstanceAuxiliaryLineSet(
	       sourceInstanceKey.getString(),
	       name, points, commands, count, displayState);
}

int
SoBRLSceneController::publishDatabaseSourceAuxiliarySourceLineSet(
    const char *sourcePath,
    const char *auxiliarySourcePath,
    const char *displayName,
    const SbVec3f *points,
    const int32_t *commands,
    int count,
    const BRLObolAuxiliaryLineSetDisplayState *displayState)
{
    return this->publishDatabaseSourceInstanceAuxiliarySourceLineSet(
	       sourcePath, auxiliarySourcePath, displayName, points, commands,
	       count, displayState);
}

int
SoBRLSceneController::publishDatabaseSourceInstanceAuxiliaryLineSet(
    const char *sourceInstanceKey,
    const char *name,
    const SbVec3f *points,
    const int32_t *commands,
    int count,
    const BRLObolAuxiliaryLineSetDisplayState *displayState)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] || !name || !name[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setAuxiliaryLineSet(name, points, commands,
			count, displayState);
    if (changed > 0) {
	if (count == 0)
	    this->advanceStructuralRevision();
	else
	    this->advanceFrameRevision();
    }
    return changed;
}

int
SoBRLSceneController::publishDatabaseSourceInstanceAuxiliarySourceLineSet(
    const char *sourceInstanceKey,
    const char *auxiliarySourcePath,
    const char *displayName,
    const SbVec3f *points,
    const int32_t *commands,
    int count,
    const BRLObolAuxiliaryLineSetDisplayState *displayState)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] ||
	!auxiliarySourcePath || !auxiliarySourcePath[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setAuxiliarySourceLineSet(
			    auxiliarySourcePath, displayName, points, commands, count,
			    displayState);
    if (changed > 0) {
	if (count == 0)
	    this->advanceStructuralRevision();
	else
	    this->advanceFrameRevision();
    }
    return changed;
}

int
SoBRLSceneController::publishDatabaseSourceExternalLineSet(
    const char *sourcePath,
    const BRLObolExternalLineSet &lineSet)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->publishDatabaseSourceInstanceExternalLineSet(
	       sourceInstanceKey.getString(), lineSet);
}

int
SoBRLSceneController::publishDatabaseSourceInstanceExternalLineSet(
    const char *sourceInstanceKey,
    const BRLObolExternalLineSet &lineSet)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int published = source->publishExternalLineSet(lineSet);
    if (published > 0)
	this->advanceStructuralRevision();
    return published;
}

int
SoBRLSceneController::publishDatabaseSourceExternalPointSet(
    const char *sourcePath,
    const BRLObolExternalPointSet &pointSet)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->publishDatabaseSourceInstanceExternalPointSet(
	       sourceInstanceKey.getString(), pointSet);
}

int
SoBRLSceneController::publishDatabaseSourceInstanceExternalPointSet(
    const char *sourceInstanceKey,
    const BRLObolExternalPointSet &pointSet)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int published = source->publishExternalPointSet(pointSet);
    if (published > 0)
	this->advanceStructuralRevision();
    return published;
}

int
SoBRLSceneController::publishDatabaseSourceExternalTriangleMesh(
    const char *sourcePath,
    const BRLObolExternalTriangleMesh &triangleMesh)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->publishDatabaseSourceInstanceExternalTriangleMesh(
	       sourceInstanceKey.getString(), triangleMesh);
}

int
SoBRLSceneController::publishDatabaseSourceInstanceExternalTriangleMesh(
    const char *sourceInstanceKey,
    const BRLObolExternalTriangleMesh &triangleMesh)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int published = source->publishExternalTriangleMesh(triangleMesh);
    if (published > 0)
	this->advanceStructuralRevision();
    return published;
}

int
SoBRLSceneController::publishDatabaseSourceExternalAnnotation(
    const char *sourcePath,
    const BRLObolExternalAnnotation &annotation)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->publishDatabaseSourceInstanceExternalAnnotation(
	       sourceInstanceKey.getString(), annotation);
}

int
SoBRLSceneController::publishDatabaseSourceInstanceExternalAnnotation(
    const char *sourceInstanceKey,
    const BRLObolExternalAnnotation &annotation)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int published = source->publishExternalAnnotation(annotation);
    if (published > 0)
	this->advanceStructuralRevision();
    return published;
}

int
SoBRLSceneController::clearDatabaseSourceExternalPrimaryGeometry(
    const char *sourcePath)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->clearDatabaseSourceInstanceExternalPrimaryGeometry(
	       sourceInstanceKey.getString());
}

int
SoBRLSceneController::clearDatabaseSourceInstanceExternalPrimaryGeometry(
    const char *sourceInstanceKey)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int cleared = source->clearExternalPrimaryGeometry();
    if (cleared > 0)
	this->advanceStructuralRevision();
    return cleared;
}

int
SoBRLSceneController::clearDatabaseSourceAuxiliaryShapes(const char *sourcePath)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->clearDatabaseSourceInstanceAuxiliaryShapes(
	       sourceInstanceKey.getString());
}

int
SoBRLSceneController::clearDatabaseSourceInstanceAuxiliaryShapes(
    const char *sourceInstanceKey)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int removed = source->clearAuxiliaryShapes();
    if (removed > 0)
	this->advanceStructuralRevision();
    return removed;
}

int
SoBRLSceneController::setShapeSourceState(const char *shapePath,
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
    SoNode *shape = scene_shape_find_path(this->root, shapePath, NULL);
    if (!shape)
	return -1;

    const int changed = scene_shape_set_source_state(shape,
			ownerSourcePath, ownerSourceRevision, ownerInputsRevision,
			ownerViewRevision, ownerRealizedRevision,
			ownerRealizedSourceRevision, ownerRealizedInputsRevision,
			ownerRealizedViewRevision, ownerRealizationStatus,
			ownerRealizationDiagnostic, ownerRealizationIdentity,
			ownerSourceStale, ownerStaleReason);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

SoBRLDatabaseSource *
SoBRLSceneController::getDatabaseSource(int index) const
{
    if (index < 0 || !this->root ||
	!this->root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    if (!this->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();
    const size_t pos = static_cast<size_t>(index);
    if (pos >= this->databaseSourceOrder.size())
	return NULL;
    return this->databaseSourceOrder[pos];
}

int
SoBRLSceneController::getDatabaseSourceCount(void) const
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return 0;

    if (!this->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();
    return static_cast<int>(this->databaseSourceOrder.size());
}

SoBRLDatabaseSource *
SoBRLSceneController::findDatabaseSource(const char *sourcePath) const
{
    if (!sourcePath || !sourcePath[0] || !this->root ||
	!this->root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    return this->findIndexedDatabaseSource(sourcePath);
}

SoBRLDatabaseSource *
SoBRLSceneController::findDatabaseSourceInstance(
    const char *sourceInstanceKey) const
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] || !this->root ||
	!this->root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    return this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
}

int
SoBRLSceneController::replaceDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int drawMode,
	uint32_t sourceRevision)
{
    return this->replaceDatabaseSourceInstance(sourcePath, sourcePath,
	    database, drawMode, sourceRevision);
}

int
SoBRLSceneController::replaceDatabaseSourceInstance(
    const char *sourceInstanceKey,
    const char *sourcePath,
    struct db_i *database,
    int drawMode,
    uint32_t sourceRevision)
{
    return this->replaceDatabaseSourceInstanceRepresentation(
	       sourceInstanceKey, sourcePath, NULL, -1, database, drawMode,
	       sourceRevision);
}

int
SoBRLSceneController::replaceDatabaseSourceInstanceRepresentation(
    const char *sourceInstanceKey,
    const char *sourcePath,
    const char *sourceRepresentationKey,
    int sourceRepresentationMode,
    struct db_i *database,
    int drawMode,
    uint32_t sourceRevision)
{
    BRLObolDatabaseSourcePublishState state;
    state.sourceInstanceKey = sourceInstanceKey;
    state.sourcePath = sourcePath;
    state.sourceRepresentationKey = sourceRepresentationKey;
    state.database = database;
    state.drawMode = drawMode;
    state.representationMode = sourceRepresentationMode;
    state.sourceRevisionValid = sourceRevision ? TRUE : FALSE;
    state.sourceRevision = sourceRevision;
    return this->publishDatabaseSourceInstance(state);
}

int
SoBRLSceneController::publishDatabaseSourceInstance(
    const BRLObolDatabaseSourcePublishState &state)
{
    BRLObolPerformanceTimer timer(BRLOBOL_PERF_SOURCE_REPLACE_US);
    if (timer.active())
	brlobol_performance_counter_add(BRLOBOL_PERF_SOURCE_REPLACE_CALLS, 1);

    const char *sourceInstanceKey = state.sourceInstanceKey;
    const char *sourcePath = state.sourcePath;
    if (!sourceInstanceKey || !sourceInstanceKey[0] ||
	!sourcePath || !sourcePath[0])
	return -1;
    if (!state.database)
	return this->removeDatabaseSourceInstance(sourceInstanceKey);
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *rootGroup = static_cast<SoGroup *>(this->root);
    const SbBool explicitTargetGroup =
	(state.targetGroupPath && state.targetGroupPath[0]) ? TRUE : FALSE;
    SoGroup *targetGroup = explicitTargetGroup ?
	this->ensureGroup(state.targetGroupPath) : rootGroup;
    if (!targetGroup)
	return -1;

    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    const int sourceExisted = source ? 1 : 0;
    SoGroup *sourceParent = NULL;
    int sourceIndex = -1;
    SbBool indexKeysUnchanged = FALSE;
    if (!source) {
	source = new SoBRLDatabaseSource;
    } else {
	sourceParent = this->findIndexedDatabaseSourceInstanceParent(
			   sourceInstanceKey);
	if (sourceParent)
	    sourceIndex = scene_group_find_child_index(sourceParent, source);
	if (!explicitTargetGroup && sourceParent)
	    targetGroup = sourceParent;
	indexKeysUnchanged =
	    database_source_instance_key_equal(source, sourceInstanceKey) &&
	    database_source_path_equal(source, sourcePath) ? TRUE : FALSE;
	if (!indexKeysUnchanged)
	    this->unindexDatabaseSource(source);
    }

    int drawMode = state.drawMode;
    if (drawMode != SoBRLDatabaseSource::SHADED)
	drawMode = SoBRLDatabaseSource::WIREFRAME;

    uint32_t sourceRevision = state.sourceRevision;
    if (!state.sourceRevisionValid || sourceRevision == 0)
	sourceRevision = source->sourceRevision.getValue() + 1;

    const char *effectiveInstanceKey =
	(sourceInstanceKey && sourceInstanceKey[0]) ?
	sourceInstanceKey : sourcePath;
    const char *effectiveRepresentationKey =
	(state.sourceRepresentationKey && state.sourceRepresentationKey[0]) ?
	state.sourceRepresentationKey : effectiveInstanceKey;

    int changed = 0;
    const int needsConfigure =
	!sourceExisted ||
	!indexKeysUnchanged ||
	source->stale.getValue() ||
	source->auxiliarySource.getValue() ||
	source->getDatabase() != state.database ||
	source->drawMode.getValue() != drawMode ||
	source->representationMode.getValue() != state.representationMode ||
	source->sourceRevision.getValue() != sourceRevision ||
	!database_source_path_equal(source, sourcePath) ||
	!database_source_instance_key_equal(source, effectiveInstanceKey) ||
	!scene_path_equal(
	    database_source_effective_representation_key(source).getString(),
	    effectiveRepresentationKey);

    if (needsConfigure) {
	source->configureDatabaseSourceInstanceRepresentation(
	    sourceInstanceKey, sourcePath, state.sourceRepresentationKey,
	    state.representationMode, state.database, drawMode,
	    sourceRevision);
	changed = 1;
    }

    if (state.sourceRevisionValid || state.inputsRevision != 0 ||
	source->visible.getValue() != state.visible ||
	source->highlighted.getValue() != state.highlighted ||
	source->lineStyle.getValue() != state.lineStyle ||
	source->lineWidth.getValue() != state.lineWidth ||
	scene_group_float_different(source->transparency.getValue(),
	    state.transparency) ||
	source->colorOverride.getValue() != state.colorOverride ||
	(state.colorOverride &&
	 !scene_group_color_equal(source->color.getValue(), state.color)) ||
	source->materialColorValid.getValue() != state.materialColorValid ||
	(state.materialColorValid &&
	 !scene_group_color_equal(source->materialColor.getValue(),
	     state.materialColor)) ||
	source->materialRevision.getValue() != state.materialRevision) {
	if (source->setDisplayState(state.sourceRevisionValid,
		sourceRevision, state.inputsRevision,
		state.visible, state.highlighted, state.lineStyle,
		state.lineWidth, state.transparency, state.colorOverride,
		state.color, state.materialColorValid, state.materialColor,
		state.materialRevision) > 0)
	    changed = 1;
    }

    if (state.roleFlagsValid &&
	source->setRealizationRoleFlags(state.roleFlags) > 0)
	changed = 1;

    if (state.viewPolicyValid &&
	source->setRealizationViewPolicy(state.viewDependent,
	    state.viewScale, state.botThreshold, state.curveScale,
	    state.pointScale) > 0)
	changed = 1;

    if (state.placementValid &&
	source->setPlacementState(state.drawMatrixValid,
	    state.drawMatrix, state.drawCenterValid, state.drawCenter,
	    state.drawSizeValid, state.drawSize) > 0)
	changed = 1;

    if (!sourceExisted) {
	targetGroup->addChild(source);
	this->indexDatabaseSource(source, targetGroup);
	this->advanceStructuralRevision();
	return 1;
    }

    if (!indexKeysUnchanged)
	this->indexDatabaseSource(source, sourceParent);

    if (sourceParent && sourceParent != targetGroup) {
	if (sourceIndex < 0)
	    sourceIndex = scene_group_find_child_index(sourceParent, source);
	if (sourceIndex < 0) {
	    SoGroup *foundParent = NULL;
	    int foundIndex = -1;
	    SoBRLDatabaseSource *found =
		find_database_source_instance_recursive(rootGroup,
		    sourceInstanceKey, &foundParent, &foundIndex);
	    if (found == source) {
		sourceParent = foundParent;
		sourceIndex = foundIndex;
	    }
	}
	if (sourceParent && sourceIndex >= 0) {
	    source->ref();
	    this->unindexDatabaseSource(source);
	    sourceParent->removeChild(sourceIndex);
	    targetGroup->addChild(source);
	    this->indexDatabaseSource(source, targetGroup);
	    source->unref();
	    this->advanceStructuralRevision();
	    return 1;
	}
    }

    if (!indexKeysUnchanged)
	this->indexDatabaseSource(source, sourceParent);

    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::renameDatabaseSource(const char *sourcePath,
	const char *newSourcePath,
	uint32_t sourceRevision)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoBRLDatabaseSource *source = this->findIndexedDatabaseSource(sourcePath);
    SbString sourceInstanceKey = database_source_effective_instance_key(source);
    if (sourceInstanceKey.getLength() == 0)
	return 0;

    return this->renameDatabaseSourceInstance(sourceInstanceKey.getString(),
	    newSourcePath,
	    newSourcePath, sourceRevision);
}

int
SoBRLSceneController::renameDatabaseSourceInstance(
    const char *sourceInstanceKey,
    const char *newSourceInstanceKey,
    const char *newSourcePath,
    uint32_t sourceRevision)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] ||
	!newSourceInstanceKey || !newSourceInstanceKey[0] ||
	!newSourcePath || !newSourcePath[0])
	return -1;
    if ((strcmp(sourceInstanceKey, newSourceInstanceKey) == 0 ||
	 strcmp(skip_leading_slash(sourceInstanceKey),
		skip_leading_slash(newSourceInstanceKey)) == 0) &&
	database_source_path_equal(this->findDatabaseSourceInstance(
				       sourceInstanceKey), newSourcePath))
	return 0;
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *sourceParent = NULL;
    int sourceIndex = -1;
    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (source) {
	sourceParent = this->findIndexedDatabaseSourceInstanceParent(
			   sourceInstanceKey);
	if (sourceParent)
	    sourceIndex = scene_group_find_child_index(sourceParent, source);
    }
    if (!source || !sourceParent || sourceIndex < 0) {
	SoGroup *rootGroup = static_cast<SoGroup *>(this->root);
	source = find_database_source_instance_recursive(rootGroup,
	    sourceInstanceKey, &sourceParent, &sourceIndex);
    }
    if (!source || !sourceParent || sourceIndex < 0)
	return 0;

    SoGroup *conflictParent = NULL;
    int conflictIndex = -1;
    SoBRLDatabaseSource *conflict =
	this->findIndexedDatabaseSourceInstance(newSourceInstanceKey);
    if (conflict) {
	conflictParent = this->findIndexedDatabaseSourceInstanceParent(
			     newSourceInstanceKey);
	if (conflictParent)
	    conflictIndex = scene_group_find_child_index(conflictParent,
		conflict);
    }
    if (conflict && conflict != source && conflictParent &&
	conflictIndex >= 0) {
	this->unindexDatabaseSource(conflict);
	conflictParent->removeChild(conflictIndex);
    }

    this->unindexDatabaseSource(source);
    const int changed = source->retargetDatabaseSourceInstance(
			    newSourceInstanceKey, newSourcePath, sourceRevision);
    this->indexDatabaseSource(source, sourceParent);
    if (changed > 0 || conflict)
	this->advanceStructuralRevision();
    return changed > 0 || conflict ? 1 : changed;
}

int
SoBRLSceneController::setDatabaseSourceState(const char *sourcePath,
	SbBool sourceRevisionValid,
	uint32_t sourceRevision,
	uint32_t inputsRevision,
	SbBool visible,
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
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceState(sourceInstanceKey.getString(),
	    sourceRevisionValid, sourceRevision, inputsRevision, visible,
	    highlighted, lineStyle, lineWidth, transparency, colorOverride,
	    color, materialColorValid, materialColor, materialRevision);
}

int
SoBRLSceneController::setDatabaseSourceInstanceState(
    const char *sourceInstanceKey,
    SbBool sourceRevisionValid,
    uint32_t sourceRevision,
    uint32_t inputsRevision,
    SbBool visible,
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
    BRLObolPerformanceTimer timer(BRLOBOL_PERF_SOURCE_STATE_US);
    if (timer.active())
	brlobol_performance_counter_add(BRLOBOL_PERF_SOURCE_STATE_CALLS, 1);

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setDisplayState(sourceRevisionValid,
			sourceRevision, inputsRevision, visible, highlighted, lineStyle,
			lineWidth, transparency, colorOverride, color, materialColorValid,
			materialColor, materialRevision);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setDatabaseSourceDisplayPatch(const char *sourcePath,
	const BRLObolDatabaseSourceDisplayPatch &patch)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceDisplayPatch(
	       sourceInstanceKey.getString(), patch);
}

int
SoBRLSceneController::setDatabaseSourceInstanceDisplayPatch(
    const char *sourceInstanceKey,
    const BRLObolDatabaseSourceDisplayPatch &patch)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->applyDisplayPatch(patch);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setDatabaseSourceDisplayName(const char *sourcePath,
	const char *displayName)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceDisplayName(
	       sourceInstanceKey.getString(),
	       displayName);
}

int
SoBRLSceneController::setDatabaseSourceInstanceDisplayName(
    const char *sourceInstanceKey,
    const char *displayName)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setDisplayNameState(displayName);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setDatabaseSourceDrawMode(const char *sourcePath,
	int drawMode)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceDrawMode(
	       sourceInstanceKey.getString(), drawMode);
}

int
SoBRLSceneController::setDatabaseSourceInstanceDrawMode(
    const char *sourceInstanceKey,
    int drawMode)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    int changed = source->setDrawModeState(drawMode);
    const int currentRepresentationMode =
	source->representationMode.getValue();
    if (currentRepresentationMode ==
	SoBRLDatabaseSource::REPRESENTATION_DEFAULT ||
	currentRepresentationMode ==
	SoBRLDatabaseSource::REPRESENTATION_WIRE ||
	currentRepresentationMode ==
	SoBRLDatabaseSource::REPRESENTATION_SHADED) {
	const int representationMode =
	    drawMode == SoBRLDatabaseSource::SHADED ?
	    SoBRLDatabaseSource::REPRESENTATION_SHADED :
	    SoBRLDatabaseSource::REPRESENTATION_WIRE;
	const SbString existingRepresentationKey =
	    source->representationKey.getValue();
	const SbString representationKey =
	    existingRepresentationKey.getLength() > 0 ?
	    existingRepresentationKey :
	    database_source_effective_instance_key(source);
	const int representationChanged = source->setRepresentationState(
					      representationKey.getString(), representationMode);
	if (representationChanged > 0)
	    changed = 1;
    }
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setDatabaseSourceInstanceRepresentation(
    const char *sourceInstanceKey,
    const char *sourceRepresentationKey,
    int sourceRepresentationMode)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setRepresentationState(
			    sourceRepresentationKey, sourceRepresentationMode);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setDatabaseSourceMaterialPolicy(const char *sourcePath,
	int materialPolicy)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceMaterialPolicy(
	       sourceInstanceKey.getString(),
	       materialPolicy);
}

int
SoBRLSceneController::setDatabaseSourceInstanceMaterialPolicy(
    const char *sourceInstanceKey,
    int materialPolicy)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setMaterialPolicyState(materialPolicy);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setDatabaseSourcePlacementState(const char *sourcePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstancePlacementState(
	       sourceInstanceKey.getString(),
	       drawMatrixValid, drawMatrix, drawCenterValid, drawCenter,
	       drawSizeValid, drawSize);
}

int
SoBRLSceneController::setDatabaseSourceInstancePlacementState(
    const char *sourceInstanceKey,
    SbBool drawMatrixValid,
    const SbMatrix &drawMatrix,
    SbBool drawCenterValid,
    const SbVec3f &drawCenter,
    SbBool drawSizeValid,
    float drawSize)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setPlacementState(drawMatrixValid,
			drawMatrix, drawCenterValid, drawCenter, drawSizeValid, drawSize);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setDatabaseSourceBoundsState(const char *sourcePath,
	SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceBoundsState(
	       sourceInstanceKey.getString(), boundsValid,
	       boundsMin, boundsMax);
}

int
SoBRLSceneController::setDatabaseSourceInstanceBoundsState(
    const char *sourceInstanceKey,
    SbBool boundsValid,
    const SbVec3f &boundsMin,
    const SbVec3f &boundsMax)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setSourceBoundsState(boundsValid,
			boundsMin, boundsMax);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::markDatabaseSourceStale(const char *sourcePath,
	uint32_t staleReason)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->markDatabaseSourceInstanceStale(sourceInstanceKey.getString(),
	    staleReason);
}

int
SoBRLSceneController::markDatabaseSourceInstanceStale(
    const char *sourceInstanceKey,
    uint32_t staleReason)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    if (!staleReason)
	staleReason = SoBRLDatabaseSource::STALE_SOURCE;

    uint32_t nextSourceRevision = source->sourceRevision.getValue();
    if (staleReason & (SoBRLDatabaseSource::STALE_SOURCE |
		       SoBRLDatabaseSource::STALE_DATABASE))
	nextSourceRevision++;

    const uint32_t nextReason =
	source->staleReason.getValue() | staleReason;
    const int changed =
	source->sourceRevision.getValue() != nextSourceRevision ||
	!source->stale.getValue() ||
	source->staleReason.getValue() != nextReason ||
	source->realizationStatus.getValue() != SoBRLDatabaseSource::UNREALIZED ||
	source->realizationDiagnostic.getValue().getLength() > 0;
    if (!changed)
	return 0;

    source->sourceRevision = nextSourceRevision;
    source->markStale(staleReason);
    this->advanceFrameRevision();
    return 1;
}

int
SoBRLSceneController::setDatabaseSourceRealizationState(const char *sourcePath,
	int realizationStatus,
	uint32_t realizedSourceRevision,
	uint32_t realizedInputsRevision,
	uint32_t staleReason,
	const char *diagnostic)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceRealizationState(
	       sourceInstanceKey.getString(),
	       realizationStatus, realizedSourceRevision, realizedInputsRevision,
	       staleReason, diagnostic);
}

int
SoBRLSceneController::setDatabaseSourceInstanceRealizationState(
    const char *sourceInstanceKey,
    int realizationStatus,
    uint32_t realizedSourceRevision,
    uint32_t realizedInputsRevision,
    uint32_t staleReason,
    const char *diagnostic)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setRealizationState(realizationStatus,
			realizedSourceRevision, realizedInputsRevision, staleReason,
			diagnostic);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setDatabaseSourceRealizationRoleFlags(
    const char *sourcePath,
    int roleFlags)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceRealizationRoleFlags(
	       sourceInstanceKey.getString(),
	       roleFlags);
}

int
SoBRLSceneController::setDatabaseSourceInstanceRealizationRoleFlags(
    const char *sourceInstanceKey,
    int roleFlags)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setRealizationRoleFlags(roleFlags);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setDatabaseSourceRealizationViewPolicy(
    const char *sourcePath,
    SbBool viewDependent,
    float viewScale,
    uint32_t botThreshold,
    float curveScale,
    float pointScale)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceRealizationViewPolicy(
	       sourceInstanceKey.getString(),
	       viewDependent, viewScale, botThreshold, curveScale, pointScale);
}

int
SoBRLSceneController::setDatabaseSourceInstanceRealizationViewPolicy(
    const char *sourceInstanceKey,
    SbBool viewDependent,
    float viewScale,
    uint32_t botThreshold,
    float curveScale,
    float pointScale)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    const int changed = source->setRealizationViewPolicy(viewDependent,
			viewScale, botThreshold, curveScale, pointScale);
    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::moveDatabaseSourceToGroup(const char *sourcePath,
	const char *groupPath)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoBRLDatabaseSource *source = this->findIndexedDatabaseSource(sourcePath);
    SbString sourceInstanceKey = database_source_effective_instance_key(source);
    if (sourceInstanceKey.getLength() == 0)
	return 0;

    return this->moveDatabaseSourceInstanceToGroup(
	       sourceInstanceKey.getString(), groupPath);
}

int
SoBRLSceneController::moveDatabaseSourceInstanceToGroup(
    const char *sourceInstanceKey,
    const char *groupPath)
{
    BRLObolPerformanceTimer timer(BRLOBOL_PERF_SOURCE_MOVE_US);
    if (timer.active())
	brlobol_performance_counter_add(BRLOBOL_PERF_SOURCE_MOVE_CALLS, 1);

    if (!sourceInstanceKey || !sourceInstanceKey[0] || !groupPath)
	return -1;
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *sourceParent = NULL;
    int sourceIndex = -1;
    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (source) {
	sourceParent = this->findIndexedDatabaseSourceInstanceParent(
			   sourceInstanceKey);
	if (sourceParent)
	    sourceIndex = scene_group_find_child_index(sourceParent, source);
    }
    if (!source || !sourceParent || sourceIndex < 0) {
	SoGroup *rootGroup = static_cast<SoGroup *>(this->root);
	source = find_database_source_instance_recursive(rootGroup,
	    sourceInstanceKey, &sourceParent, &sourceIndex);
    }
    if (!source || !sourceParent || sourceIndex < 0)
	return 0;

    SoGroup *targetGroup = this->ensureGroup(groupPath);
    if (!targetGroup)
	return -1;
    if (targetGroup == sourceParent)
	return 0;

    source->ref();
    this->unindexDatabaseSource(source);
    sourceParent->removeChild(sourceIndex);
    targetGroup->addChild(source);
    this->indexDatabaseSource(source, targetGroup);
    source->unref();
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::removeDatabaseSource(const char *sourcePath)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoBRLDatabaseSource *source = this->findIndexedDatabaseSource(sourcePath);
    SbString sourceInstanceKey = database_source_effective_instance_key(source);
    if (sourceInstanceKey.getLength() == 0)
	return 0;

    return this->removeDatabaseSourceInstance(sourceInstanceKey.getString());
}

int
SoBRLSceneController::removeDatabaseSourceInstance(
    const char *sourceInstanceKey)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return 0;
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *sourceParent = NULL;
    int childIndex = -1;
    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (source) {
	sourceParent = this->findIndexedDatabaseSourceInstanceParent(
			   sourceInstanceKey);
	if (sourceParent)
	    childIndex = scene_group_find_child_index(sourceParent, source);
    }
    if (!source || !sourceParent || childIndex < 0) {
	SoGroup *group = static_cast<SoGroup *>(this->root);
	source = find_database_source_instance_recursive(group, sourceInstanceKey,
	    &sourceParent, &childIndex);
    }
    if (!sourceParent || childIndex < 0)
	return 0;

    SoNode *node = sourceParent->getChild(childIndex);
    if (node && node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	this->unindexDatabaseSource(static_cast<SoBRLDatabaseSource *>(node));
    sourceParent->removeChild(childIndex);
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::clearDatabaseSources(void)
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(this->root);
    int removed = clear_database_sources_recursive(group);
    if (removed > 0) {
	this->clearDatabaseSourceIndex();
	this->advanceStructuralRevision();
    }
    return removed;
}

static int
scene_path_component_count(const char *path);

static SbString
database_source_parent_path_from_source_path(const char *sourcePath)
{
    const char *start = skip_leading_slash(sourcePath);
    if (!start || !start[0])
	return "/";

    size_t len = strlen(start);
    while (len > 0 && start[len - 1] == '/')
	len--;
    if (len == 0)
	return "/";

    const char *slash = NULL;
    for (size_t i = 0; i < len; i++) {
	if (start[i] == '/')
	    slash = &start[i];
    }
    if (!slash)
	return "/";

    return SbString(std::string(start, (size_t)(slash - start)).c_str());
}

SbBool
SoBRLSceneController::databaseSourceSummaryForSource(
    SoBRLDatabaseSource *source,
    BRLObolDatabaseSourceSummary &summary) const
{
    summary = BRLObolDatabaseSourceSummary();
    if (!source || !source->getSummary(summary) || !summary.valid)
	return FALSE;

    const SbString instanceKey = database_source_effective_instance_key(source);
    SoGroup *parent = this->findIndexedDatabaseSourceInstanceParent(
			  instanceKey.getString());
    SbString parentGroupPath("");
    if (parent) {
	parentGroupPath = scene_group_summary_path(parent, "");
	if (parentGroupPath.getLength() == 0 &&
	    parent == scene_root_group(this->root))
	    parentGroupPath = "/";
    }
    if (parentGroupPath.getLength() == 0)
	parentGroupPath = database_source_parent_path_from_source_path(
			      source->path.getValue().getString());

    summary.hasParent = TRUE;
    summary.parentGroupPath = parentGroupPath;
    summary.drawTreeDepth =
	scene_path_component_count(parentGroupPath.getString()) + 1;
    return TRUE;
}

SbBool
SoBRLSceneController::getDatabaseSourceSummary(int index,
	BRLObolDatabaseSourceSummary &summary) const
{
    summary = BRLObolDatabaseSourceSummary();
    if (index < 0 || !this->root ||
	!this->root->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    return this->databaseSourceSummaryForSource(
	       this->getDatabaseSource(index), summary);
}

SbBool
SoBRLSceneController::getDatabaseSourceSummaryForPath(
    const char *sourcePath,
    BRLObolDatabaseSourceSummary &summary) const
{
    summary = BRLObolDatabaseSourceSummary();
    if (!sourcePath || !sourcePath[0])
	return FALSE;
    return this->databaseSourceSummaryForSource(
	       this->findIndexedDatabaseSource(sourcePath), summary);
}

SbBool
SoBRLSceneController::getDatabaseSourceSummaryForInstance(
    const char *sourceInstanceKey,
    BRLObolDatabaseSourceSummary &summary) const
{
    summary = BRLObolDatabaseSourceSummary();
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return FALSE;
    return this->databaseSourceSummaryForSource(
	       this->findIndexedDatabaseSourceInstance(sourceInstanceKey),
	       summary);
}

int
SoBRLSceneController::getRealizedShapeSummaryCount(void) const
{
    int count = 0;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (source)
	    count += source->getRealizedShapeSummaryCount();
    }
    return count;
}

SbBool
SoBRLSceneController::getRealizedShapeSummary(int index,
	BRLObolRealizedShapeSummary &summary) const
{
    summary = BRLObolRealizedShapeSummary();
    if (index < 0)
	return FALSE;

    int remaining = index;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (!source)
	    continue;

	const int sourceShapeCount =
	    source->getRealizedShapeSummaryCount();
	if (remaining < sourceShapeCount) {
	    if (!source->getRealizedShapeSummary(remaining, summary))
		return FALSE;
	    summary.ownerSourceIndex = i;
	    if (summary.ownerSourceInstanceKey.getLength() == 0)
		summary.ownerSourceInstanceKey =
		    database_source_effective_instance_key(source);
	    return TRUE;
	}
	remaining -= sourceShapeCount;
    }

    return FALSE;
}

int
SoBRLSceneController::getRealizedMaterialSummaryCount(void) const
{
    int count = 0;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (source)
	    count += source->getRealizedMaterialSummaryCount();
    }
    return count;
}

SbBool
SoBRLSceneController::getRealizedMaterialSummary(int index,
	BRLObolRealizedMaterialSummary &summary) const
{
    summary = BRLObolRealizedMaterialSummary();
    if (index < 0)
	return FALSE;

    int remaining = index;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (!source)
	    continue;

	const int sourceMaterialCount =
	    source->getRealizedMaterialSummaryCount();
	if (remaining < sourceMaterialCount) {
	    if (!source->getRealizedMaterialSummary(remaining, summary))
		return FALSE;
	    summary.ownerSourceIndex = i;
	    if (summary.ownerSourceInstanceKey.getLength() == 0)
		summary.ownerSourceInstanceKey =
		    database_source_effective_instance_key(source);
	    return TRUE;
	}
	remaining -= sourceMaterialCount;
    }

    return FALSE;
}

SbBool
SoBRLSceneController::getRealizedMaterialProperty(int materialIndex,
	int propertyIndex, SbString &groupOut, SbString &nameOut,
	SbString &valueOut) const
{
    if (materialIndex < 0 || propertyIndex < 0)
	return FALSE;

    int remaining = materialIndex;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (!source)
	    continue;

	const int sourceMaterialCount =
	    source->getRealizedMaterialSummaryCount();
	if (remaining < sourceMaterialCount)
	    return source->getRealizedMaterialProperty(remaining,
		    propertyIndex, groupOut, nameOut, valueOut);
	remaining -= sourceMaterialCount;
    }

    return FALSE;
}

static void
scene_tree_summary_fill(const SoNode *node, int depth, SbBool hasParent,
			int ownerSourceIndex, const SbString &ownerSourcePath,
			const SbString &nodePath, BRLObolSceneTreeSummary &summary)
{
    summary = BRLObolSceneTreeSummary();
    if (!node)
	return;

    summary.valid = TRUE;
    summary.hasParent = hasParent;
    summary.drawTreeDepth = depth;
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    if (summary.ownerSourcePath.getLength() == 0) {
	const SbString retainedOwnerPath = scene_shape_owner_source_path(node);
	if (retainedOwnerPath.getLength() > 0)
	    summary.ownerSourcePath = retainedOwnerPath;
    }
    summary.ownerSourceInstanceKey =
	scene_shape_owner_source_instance_key(node);
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
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	summary.nodeKind = BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE;
	summary.path = source->path.getValue();
	summary.displayName = source->displayName.getValue();
	if (summary.ownerSourcePath.getLength() == 0)
	    summary.ownerSourcePath = source->path.getValue();
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		database_source_effective_instance_key(source);
	return;
    }

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	const SoBRLVListShape *shape =
	    static_cast<const SoBRLVListShape *>(node);
	summary.nodeKind = BRLObolSceneTreeSummary::NODE_VLIST_SHAPE;
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
	summary.nodeKind = BRLObolSceneTreeSummary::NODE_MESH_SHAPE;
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
	summary.nodeKind = BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT;
	summary.path = object->sourcePath.getValue();
	summary.sourceName = object->sourceName.getValue();
	summary.sourceType = object->sourceType.getValue();
	summary.sourceId = object->sourceId.getValue();
	summary.displayName = object->materialName.getValue();
	return;
    }

    summary.nodeKind = summary.isGroup ?
		       BRLObolSceneTreeSummary::NODE_GROUP :
		       BRLObolSceneTreeSummary::NODE_OTHER;
    summary.path = scene_group_summary_path(node, nodePath);
}


static const SoNode *
scene_public_realized_shape_node(const SoNode *node)
{
    if (!node)
	return NULL;

    if (node->isOfType(SoBRLVListShape::getClassTypeId()) ||
	node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return node;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *found =
		scene_public_realized_shape_node(group->getChild(i));
	    if (found)
		return found;
	}
    }

    return NULL;
}


static const SoNode *
scene_public_realized_material_node(const SoNode *node)
{
    if (!node)
	return NULL;

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return node;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *found =
		scene_public_realized_material_node(group->getChild(i));
	    if (found)
		return found;
	}
    }

    return NULL;
}


static const SoNode *
scene_public_realized_child_node(const SoNode *node)
{
    const SoNode *shape = scene_public_realized_shape_node(node);
    if (shape)
	return shape;

    const SoNode *material = scene_public_realized_material_node(node);
    if (material)
	return material;

    return node;
}


static int
scene_tree_summary_node_count(const SoNode *node)
{
    if (!node)
	return 0;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	return source->getRealizedTreeSummaryCount();
    }

    int count = 1;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    count += scene_tree_summary_node_count(group->getChild(i));
    }
    return count;
}

static SbBool
find_scene_tree_summary_in_node(const SoNode *node, int &index, int depth,
				SbBool hasParent, int ownerSourceIndex,
				const SbString &ownerSourcePath, const SbString &nodePath,
				BRLObolSceneTreeSummary &summary)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	const int sourceTreeCount = source->getRealizedTreeSummaryCount();
	if (index >= sourceTreeCount) {
	    index -= sourceTreeCount;
	    return FALSE;
	}
	const int sourceTreeIndex = index;
	if (!source->getRealizedTreeSummary(sourceTreeIndex, summary))
	    return FALSE;
	summary.drawTreeDepth += depth;
	summary.hasParent = hasParent ? TRUE : summary.hasParent;
	summary.ownerSourceIndex = ownerSourceIndex;
	if (ownerSourcePath.getLength() > 0)
	    summary.ownerSourcePath = ownerSourcePath;
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		database_source_effective_instance_key(source);
	return TRUE;
    }

    if (index == 0) {
	scene_tree_summary_fill(node, depth, hasParent, ownerSourceIndex,
				ownerSourcePath, nodePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *child = group->getChild(i);
	    const SbString childPath =
		scene_child_summary_path(nodePath, child);
	    if (find_scene_tree_summary_in_node(child, index,
						depth + 1, TRUE, ownerSourceIndex, ownerSourcePath,
						childPath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLSceneController::getSceneTreeSummaryCount(void) const
{
    return scene_tree_summary_node_count(this->root);
}

SbBool
SoBRLSceneController::getSceneTreeSummary(int index,
	BRLObolSceneTreeSummary &summary) const
{
    summary = BRLObolSceneTreeSummary();
    if (index < 0 || !this->root)
	return FALSE;

    if (index == 0) {
	scene_tree_summary_fill(this->root, 0, FALSE, -1, "", "/",
				summary);
	return TRUE;
    }
    index--;

    if (!this->root->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    const SoGroup *group = static_cast<const SoGroup *>(this->root);
    int sourceIndex = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *child = group->getChild(i);
	int childOwnerIndex = -1;
	SbString childOwnerPath("");
	if (child &&
	    child->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    const SoBRLDatabaseSource *source =
		static_cast<const SoBRLDatabaseSource *>(child);
	    childOwnerIndex = sourceIndex;
	    childOwnerPath = source->path.getValue();
	    sourceIndex++;
	}

	const SbString childPath = scene_child_summary_path("/", child);
	if (find_scene_tree_summary_in_node(child, index, 1, TRUE,
					    childOwnerIndex, childOwnerPath, childPath, summary))
	    return TRUE;
    }

    return FALSE;
}

static int
scene_summary_path_equal(const SbString &summaryPath, const char *nodePath)
{
    return scene_path_equal(summaryPath.getString(), nodePath);
}

static int
scene_path_component_count(const char *path)
{
    if (!path || !path[0])
	return 0;

    int count = 0;
    const char *cp = path;
    while (*cp) {
	while (*cp == '/')
	    cp++;
	if (!*cp)
	    break;
	count++;
	while (*cp && *cp != '/')
	    cp++;
    }

    return count;
}

static const SoGroup *
scene_tree_group_find_by_summary_path(const SoNode *node,
				      const char *nodePath,
				      const SbString &fallbackPath)
{
    if (!node || !nodePath || !node->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    const SoGroup *group = static_cast<const SoGroup *>(node);
    const SbString groupPath = scene_group_summary_path(node, fallbackPath);
    if (scene_path_equal(groupPath.getString(), nodePath))
	return group;

    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *child = group->getChild(i);
	if (!child || !child->isOfType(SoGroup::getClassTypeId()))
	    continue;
	const SbString childPath = scene_child_summary_path(groupPath, child);
	const SoGroup *found = scene_tree_group_find_by_summary_path(child,
			       nodePath, childPath);
	if (found)
	    return found;
    }

    return NULL;
}

SbBool
SoBRLSceneController::getSceneTreeSummaryForPath(const char *nodePath,
	BRLObolSceneTreeSummary &summary) const
{
    summary = BRLObolSceneTreeSummary();
    if (!this->root || !nodePath)
	return FALSE;

    const char *normalizedPath = skip_leading_slash(nodePath);
    if (!nodePath[0] || !normalizedPath[0]) {
	scene_tree_summary_fill(this->root, 0, FALSE, -1, "", "/",
	    summary);
	return summary.valid;
    }

    SoBRLDatabaseSource *source = this->findIndexedDatabaseSource(nodePath);
    if (source) {
	const int depth = scene_path_component_count(nodePath);
	scene_tree_summary_fill(source, depth, depth > 0 ? TRUE : FALSE, -1,
	    source->path.getValue(), source->path.getValue(), summary);
	return summary.valid;
    }

    const SoGroup *group = scene_group_find_path_const(this->root, nodePath);
    if (group) {
	const int depth = scene_path_component_count(nodePath);
	scene_tree_summary_fill(group, depth, depth > 0 ? TRUE : FALSE, -1,
	    "", nodePath, summary);
	return summary.valid;
    }

    const SoNode *shape = scene_shape_find_path(this->root, nodePath, NULL);
    if (shape) {
	const int depth = scene_path_component_count(nodePath);
	scene_tree_summary_fill(shape, depth, depth > 0 ? TRUE : FALSE, -1,
	    "", nodePath, summary);
	return summary.valid;
    }

    BRLObolSceneTreeSummary candidate;
    BRLObolSceneTreeSummary fallback;
    const int count = this->getSceneTreeSummaryCount();
    for (int i = 0; i < count; i++) {
	if (!this->getSceneTreeSummary(i, candidate) ||
	    !candidate.valid ||
	    !scene_summary_path_equal(candidate.path, nodePath))
	    continue;
	if (candidate.nodeKind ==
	    BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE) {
	    summary = candidate;
	    return TRUE;
	}
	if (!fallback.valid)
	    fallback = candidate;
    }

    if (fallback.valid) {
	summary = fallback;
	return TRUE;
    }

    return FALSE;
}

SbBool
SoBRLSceneController::getSceneChildTreeSummary(const char *nodePath,
	int childIndex,
	BRLObolSceneTreeSummary &summary) const
{
    summary = BRLObolSceneTreeSummary();
    if (!this->root || !nodePath || childIndex < 0)
	return FALSE;

    BRLObolSceneTreeSummary parentSummary;
    if (!this->getSceneTreeSummaryForPath(nodePath, parentSummary))
	return FALSE;

    SoBRLDatabaseSource *source = this->findDatabaseSource(nodePath);
    if (source) {
	if (childIndex >= source->getNumChildren())
	    return FALSE;

	const SoNode *child = source->getChild(childIndex);
	const SoNode *publicChild =
	    (child && child->isOfType(SoBRLDatabaseSource::getClassTypeId())) ?
	    child : scene_public_realized_child_node(child);
	const SbString childPath =
	    scene_child_summary_path(parentSummary.path, child);
	scene_tree_summary_fill(publicChild,
				parentSummary.drawTreeDepth + 1, TRUE,
				parentSummary.ownerSourceIndex,
				source->path.getValue(), childPath, summary);
	if (!summary.valid)
	    return FALSE;
	if (summary.ownerSourcePath.getLength() == 0)
	    summary.ownerSourcePath = source->path.getValue();
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		database_source_effective_instance_key(source);
	if (summary.path.getLength() == 0)
	    summary.path = source->path.getValue();
	return TRUE;
    }

    const SoGroup *group = this->findGroup(nodePath);
    if (!group)
	group = scene_tree_group_find_by_summary_path(this->root, nodePath,
		SbString("/"));
    if (group) {
	if (childIndex >= group->getNumChildren())
	    return FALSE;

	const SoNode *child = group->getChild(childIndex);
	const SbString childPath = scene_child_summary_path(
				       parentSummary.path, child);
	scene_tree_summary_fill(child, parentSummary.drawTreeDepth + 1,
				TRUE, -1, "", childPath, summary);
	return summary.valid;
    }

    return FALSE;
}

static void
scene_display_summary_fill_common(BRLObolSceneDisplaySummary &summary,
				  int nodeKind, int ownerSourceIndex, const SbString &ownerSourcePath,
				  const SbString &nodePath)
{
    summary = BRLObolSceneDisplaySummary();
    summary.valid = TRUE;
    summary.nodeKind = nodeKind;
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    summary.path = nodePath;
}

template <typename ShapeT>
static void
scene_display_summary_fill_shape(const ShapeT *shape, int nodeKind,
				 int ownerSourceIndex, const SbString &ownerSourcePath,
				 BRLObolSceneDisplaySummary &summary)
{
    SbString effectiveOwnerPath = ownerSourcePath;
    if (shape && effectiveOwnerPath.getLength() == 0 &&
	shape->ownerSourcePath.getValue().getLength() > 0)
	effectiveOwnerPath = shape->ownerSourcePath.getValue();
    scene_display_summary_fill_common(summary, nodeKind, ownerSourceIndex,
				      effectiveOwnerPath,
				      shape ? shape->sourcePath.getValue() : SbString(""));
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
scene_display_summary_fill(const SoNode *node, int ownerSourceIndex,
			   const SbString &ownerSourcePath, const SbString &nodePath,
			   BRLObolSceneDisplaySummary &summary)
{
    summary = BRLObolSceneDisplaySummary();
    if (!node)
	return;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	scene_display_summary_fill_common(summary,
					  BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE,
					  ownerSourceIndex, ownerSourcePath, source->path.getValue());
	summary.ownerSourceInstanceKey =
	    database_source_effective_instance_key(source);
	summary.isDatabaseSource = TRUE;
	summary.hasDrawIntent = source->path.getValue().getLength() > 0;
	summary.intentPath = source->path.getValue();
	summary.intentDrawMode = source->drawMode.getValue();
	summary.visible = source->visible.getValue();
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
	scene_display_summary_fill_shape(
	    static_cast<const SoBRLVListShape *>(node),
	    BRLObolSceneTreeSummary::NODE_VLIST_SHAPE,
	    ownerSourceIndex, ownerSourcePath, summary);
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	scene_display_summary_fill_shape(
	    static_cast<const SoBRLMeshShape *>(node),
	    BRLObolSceneTreeSummary::NODE_MESH_SHAPE,
	    ownerSourceIndex, ownerSourcePath, summary);
	return;
    }

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId())) {
	const SoBRLMaterialObject *object =
	    static_cast<const SoBRLMaterialObject *>(node);
	scene_display_summary_fill_common(summary,
					  BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT,
					  ownerSourceIndex, ownerSourcePath,
					  object->sourcePath.getValue());
	return;
    }

    if (node->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	const SoBRLSceneGroup *group =
	    static_cast<const SoBRLSceneGroup *>(node);
	const SbString retainedPath =
	    scene_group_summary_path(node, nodePath);
	scene_display_summary_fill_common(summary,
					  BRLObolSceneTreeSummary::NODE_GROUP,
					  ownerSourceIndex, ownerSourcePath, retainedPath);
	summary.hasDrawIntent = group->drawIntentValid.getValue();
	if (summary.hasDrawIntent) {
	    summary.intentPath = group->drawIntentPath.getValue();
	    if (summary.intentPath.getLength() == 0)
		summary.intentPath = retainedPath;
	    summary.intentDrawMode = group->drawMode.getValue();
	}
	summary.visible = group->visible.getValue();
	summary.highlighted = group->highlighted.getValue();
	summary.lineStyle = group->lineStyle.getValue();
	summary.lineWidth = group->lineWidth.getValue();
	summary.transparency = group->transparency.getValue();
	summary.drawMode = group->drawMode.getValue();
	summary.materialValid =
	    group->materialColorValid.getValue() ||
	    group->colorOverride.getValue();
	summary.materialRevision = group->materialRevision.getValue();
	if (group->materialColorValid.getValue())
	    summary.materialColor = group->materialColor.getValue();
	else if (group->colorOverride.getValue())
	    summary.materialColor = group->color.getValue();
	return;
    }

    const int nodeKind = node->isOfType(SoGroup::getClassTypeId()) ?
			 BRLObolSceneTreeSummary::NODE_GROUP :
			 BRLObolSceneTreeSummary::NODE_OTHER;
    scene_display_summary_fill_common(summary, nodeKind, ownerSourceIndex,
				      ownerSourcePath, nodePath);
}

static SbBool
find_scene_display_summary_in_node(const SoNode *node, int &index,
				   int ownerSourceIndex, const SbString &ownerSourcePath,
				   const SbString &nodePath, BRLObolSceneDisplaySummary &summary)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	const int sourceDisplayCount =
	    source->getRealizedDisplaySummaryCount();
	if (index >= sourceDisplayCount) {
	    index -= sourceDisplayCount;
	    return FALSE;
	}
	const int sourceDisplayIndex = index;
	if (!source->getRealizedDisplaySummary(sourceDisplayIndex, summary))
	    return FALSE;
	summary.ownerSourceIndex = ownerSourceIndex;
	if (ownerSourcePath.getLength() > 0)
	    summary.ownerSourcePath = ownerSourcePath;
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		database_source_effective_instance_key(source);
	return TRUE;
    }

    if (index == 0) {
	scene_display_summary_fill(node, ownerSourceIndex, ownerSourcePath,
				   nodePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *child = group->getChild(i);
	    const SbString childPath =
		scene_child_summary_path(nodePath, child);
	    if (find_scene_display_summary_in_node(child, index,
						   ownerSourceIndex, ownerSourcePath, childPath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLSceneController::getSceneDisplaySummaryCount(void) const
{
    return this->getSceneTreeSummaryCount();
}

SbBool
SoBRLSceneController::getSceneDisplaySummary(int index,
	BRLObolSceneDisplaySummary &summary) const
{
    summary = BRLObolSceneDisplaySummary();
    if (index < 0 || !this->root)
	return FALSE;

    if (index == 0) {
	scene_display_summary_fill(this->root, -1, "", "/", summary);
	return TRUE;
    }
    index--;

    if (!this->root->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    const SoGroup *group = static_cast<const SoGroup *>(this->root);
    int sourceIndex = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *child = group->getChild(i);
	int childOwnerIndex = -1;
	SbString childOwnerPath("");
	if (child &&
	    child->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    const SoBRLDatabaseSource *source =
		static_cast<const SoBRLDatabaseSource *>(child);
	    childOwnerIndex = sourceIndex;
	    childOwnerPath = source->path.getValue();
	    sourceIndex++;
	}

	const SbString childPath = scene_child_summary_path("/", child);
	if (find_scene_display_summary_in_node(child, index,
					       childOwnerIndex, childOwnerPath, childPath, summary))
	    return TRUE;
    }

    return FALSE;
}

static void
scene_material_summary_from_display(const BRLObolSceneDisplaySummary &display,
				    BRLObolSceneMaterialSummary &summary)
{
    summary = BRLObolSceneMaterialSummary();
    if (!display.valid)
	return;

    summary.valid = TRUE;
    summary.nodeKind = display.nodeKind;
    summary.materialValid =
	(display.nodeKind == BRLObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	 display.nodeKind == BRLObolSceneTreeSummary::NODE_MESH_SHAPE) ?
	display.materialValid : FALSE;
    summary.materialRevision = display.materialRevision;
    summary.materialColor = display.materialColor;
    summary.ownerSourceIndex = display.ownerSourceIndex;
    summary.ownerSourcePath = display.ownerSourcePath;
    summary.ownerSourceInstanceKey = display.ownerSourceInstanceKey;
    summary.path = display.path;
}

int
SoBRLSceneController::getSceneMaterialSummaryCount(void) const
{
    return this->getSceneDisplaySummaryCount();
}

SbBool
SoBRLSceneController::getSceneMaterialSummary(int index,
	BRLObolSceneMaterialSummary &summary) const
{
    summary = BRLObolSceneMaterialSummary();
    BRLObolSceneDisplaySummary display;
    if (!this->getSceneDisplaySummary(index, display))
	return FALSE;

    scene_material_summary_from_display(display, summary);
    return TRUE;
}

static int
scene_bounds_node_kind(const SoNode *node)
{
    if (!node)
	return BRLObolSceneTreeSummary::NODE_UNKNOWN;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return BRLObolSceneTreeSummary::NODE_VLIST_SHAPE;
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return BRLObolSceneTreeSummary::NODE_MESH_SHAPE;
    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT;
    if (node->isOfType(SoGroup::getClassTypeId()))
	return BRLObolSceneTreeSummary::NODE_GROUP;
    return BRLObolSceneTreeSummary::NODE_OTHER;
}

static SbString
scene_bounds_node_path(const SoNode *node)
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
scene_bounds_for_vlist_shape(const SoBRLVListShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    for (int i = 0; i < shape->point.getNum(); i++)
	bounds.extendBy(shape->point[i]);
    return shape->point.getNum() > 0;
}

static SbBool
scene_bounds_for_mesh_shape(const SoBRLMeshShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    for (int i = 0; i < shape->point.getNum(); i++)
	bounds.extendBy(shape->point[i]);
    return shape->point.getNum() > 0;
}

static SbBox3f
scene_bounds_transform_box(const SbBox3f &bounds, const SbMatrix &matrix)
{
    SbBox3f transformed;
    transformed.makeEmpty();
    if (bounds.isEmpty())
	return transformed;

    const SbVec3f bmin = bounds.getMin();
    const SbVec3f bmax = bounds.getMax();
    for (int xi = 0; xi < 2; xi++) {
	for (int yi = 0; yi < 2; yi++) {
	    for (int zi = 0; zi < 2; zi++) {
		const SbVec3f corner(
		    xi ? bmax[0] : bmin[0],
		    yi ? bmax[1] : bmin[1],
		    zi ? bmax[2] : bmin[2]);
		SbVec3f transformedCorner;
		matrix.multVecMatrix(corner, transformedCorner);
		transformed.extendBy(transformedCorner);
	    }
	}
    }

    return transformed;
}

static SbBool
scene_node_is_overlay_intent(const SoNode *node)
{
    if (!node)
	return FALSE;
    if (node->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return static_cast<const SoBRLSceneGroup *>(node)->
	       overlayIntent.getValue();
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->
	       overlayIntent.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->
	       overlayIntent.getValue();
    return FALSE;
}

static SbBool
scene_database_source_uses_realized_placement(
    const SoBRLDatabaseSource *source)
{
    if (!source ||
	source->realizationStatus.getValue() !=
	SoBRLDatabaseSource::REALIZED ||
	(source->realizationRoleFlags.getValue() &
	 SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL))
	return FALSE;

    return (source->getRealizedShapeCount() > 0 ||
	    source->getRealizedMeshCount() > 0 ||
	    source->getRealizedMaterialObjectCount() > 0) ? TRUE : FALSE;
}

static SbBool
scene_bounds_for_node_transformed(const SoNode *node, const SbMatrix &matrix,
				  SbBox3f &bounds, SbBool includeOverlays)
{
    bounds.makeEmpty();
    if (!node)
	return FALSE;

    if (!includeOverlays && scene_node_is_overlay_intent(node))
	return FALSE;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SbBox3f localBounds;
	if (!scene_bounds_for_vlist_shape(
		static_cast<const SoBRLVListShape *>(node), localBounds))
	    return FALSE;
	bounds = scene_bounds_transform_box(localBounds, matrix);
	return bounds.isEmpty() ? FALSE : TRUE;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SbBox3f localBounds;
	if (!scene_bounds_for_mesh_shape(
		static_cast<const SoBRLMeshShape *>(node), localBounds))
	    return FALSE;
	bounds = scene_bounds_transform_box(localBounds, matrix);
	return bounds.isEmpty() ? FALSE : TRUE;
    }

    SbBool valid = FALSE;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	SbMatrix childMatrix = matrix;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    const SoBRLDatabaseSource *source =
		static_cast<const SoBRLDatabaseSource *>(node);
	    SbBool hasSourceTransform = FALSE;
	    for (int i = 0; i < group->getNumChildren(); i++) {
		const SoNode *child = group->getChild(i);
		if (child && child->isOfType(
			SoMatrixTransform::getClassTypeId())) {
		    hasSourceTransform = TRUE;
		    break;
		}
	    }
	    if (!hasSourceTransform &&
		!scene_database_source_uses_realized_placement(source) &&
		source->drawMatrixValid.getValue())
		childMatrix.multRight(source->drawMatrix.getValue());
	}
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *child = group->getChild(i);
	    if (child &&
		child->isOfType(SoMatrixTransform::getClassTypeId())) {
		const SoMatrixTransform *transform =
		    static_cast<const SoMatrixTransform *>(child);
		childMatrix.multRight(transform->matrix.getValue());
		continue;
	    }
	    SbBox3f childBounds;
	    if (scene_bounds_for_node_transformed(child, childMatrix,
						  childBounds, includeOverlays)) {
		bounds.extendBy(childBounds);
		valid = TRUE;
	    }
	}
    }

    if (valid)
	return TRUE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	SbBox3f localBounds;
	if (!source->getSourceBounds(localBounds))
	    return FALSE;

	SbMatrix sourceMatrix = matrix;
	if (source->drawMatrixValid.getValue())
	    sourceMatrix.multRight(source->drawMatrix.getValue());
	bounds = scene_bounds_transform_box(localBounds, sourceMatrix);
	return bounds.isEmpty() ? FALSE : TRUE;
    }

    return valid;
}

static SbBool
scene_bounds_for_node(const SoNode *node, SbBox3f &bounds,
		      SbBool includeOverlays)
{
    return scene_bounds_for_node_transformed(node, SbMatrix::identity(),
	    bounds, includeOverlays);
}

static SbBox3f
scene_autoview_padded_bounds(const SbBox3f &bounds)
{
    SbBox3f padded;
    padded.makeEmpty();
    if (bounds.isEmpty())
	return padded;

    const SbVec3f bmin = bounds.getMin();
    const SbVec3f bmax = bounds.getMax();
    const SbVec3f center = bounds.getCenter();
    float size = bmax[0] - bmin[0];
    if (bmax[1] - bmin[1] > size)
	size = bmax[1] - bmin[1];
    if (bmax[2] - bmin[2] > size)
	size = bmax[2] - bmin[2];

    padded.extendBy(SbVec3f(center[0] - size, center[1] - size,
	    center[2] - size));
    padded.extendBy(SbVec3f(center[0] + size, center[1] + size,
	    center[2] + size));
    return padded;
}

static SbBool
scene_database_source_bounds_recursive(const SoNode *node,
				       SbBox3f &bounds,
				       SbBool padForAutoview)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	SbBox3f sourceBounds;
	if (!scene_bounds_for_node(source, sourceBounds, TRUE))
	    (void)source->getEffectiveSourceBounds(sourceBounds);
	if (sourceBounds.isEmpty())
	    return FALSE;
	if (padForAutoview)
	    sourceBounds = scene_autoview_padded_bounds(sourceBounds);
	bounds.extendBy(sourceBounds);
	return sourceBounds.isEmpty() ? FALSE : TRUE;
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    SbBool valid = FALSE;
    const SoGroup *group = static_cast<const SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (scene_database_source_bounds_recursive(group->getChild(i),
		bounds, padForAutoview))
	    valid = TRUE;
    }

    return valid;
}

SbBool
SoBRLSceneController::getDatabaseSourceBounds(SbBox3f &bounds,
	SbBool padForAutoview) const
{
    bounds.makeEmpty();
    return scene_database_source_bounds_recursive(this->root, bounds,
	    padForAutoview);
}

static void
scene_bounds_summary_fill(const SoNode *node, int ownerSourceIndex,
			  const SbString &ownerSourcePath, const SbString &nodePath,
			  BRLObolSceneBoundsSummary &summary)
{
    summary = BRLObolSceneBoundsSummary();
    if (!node)
	return;

    summary.valid = TRUE;
    summary.nodeKind = scene_bounds_node_kind(node);
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    if (summary.ownerSourcePath.getLength() == 0) {
	const SbString retainedOwnerPath = scene_shape_owner_source_path(node);
	if (retainedOwnerPath.getLength() > 0)
	    summary.ownerSourcePath = retainedOwnerPath;
    }
    summary.ownerSourceInstanceKey =
	scene_shape_owner_source_instance_key(node);
    if (summary.ownerSourceInstanceKey.getLength() == 0 &&
	node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	summary.ownerSourceInstanceKey =
	    database_source_effective_instance_key(source);
    }
    summary.path = summary.nodeKind == BRLObolSceneTreeSummary::NODE_GROUP ?
		   scene_group_summary_path(node, nodePath) : scene_bounds_node_path(node);
    summary.boundsValid = scene_bounds_for_node(node, summary.bounds, TRUE);
}

SbBool
SoBRLSceneController::getSceneSubtreeBounds(const char *nodePath,
	SbBool includeOverlays,
	SbBox3f &bounds) const
{
    bounds.makeEmpty();
    if (!this->root)
	return FALSE;

    const char *path = nodePath ? nodePath : "/";
    const char *normalizedPath = skip_leading_slash(path);
    const SoNode *node = NULL;

    if (!path[0] || !normalizedPath[0])
	node = this->root;
    if (!node)
	node = this->findDatabaseSource(path);
    if (!node)
	node = this->findGroup(path);
    if (!node)
	node = this->findShape(path);

    return scene_bounds_for_node(node, bounds, includeOverlays);
}

static SbBool
find_scene_bounds_summary_in_node(const SoNode *node, int &index,
				  int ownerSourceIndex, const SbString &ownerSourcePath,
				  const SbString &nodePath, BRLObolSceneBoundsSummary &summary)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	const int sourceBoundsCount =
	    source->getRealizedBoundsSummaryCount();
	if (index >= sourceBoundsCount) {
	    index -= sourceBoundsCount;
	    return FALSE;
	}
	const int sourceBoundsIndex = index;
	if (!source->getRealizedBoundsSummary(sourceBoundsIndex, summary))
	    return FALSE;
	summary.ownerSourceIndex = ownerSourceIndex;
	if (ownerSourcePath.getLength() > 0)
	    summary.ownerSourcePath = ownerSourcePath;
	return TRUE;
    }

    if (index == 0) {
	scene_bounds_summary_fill(node, ownerSourceIndex, ownerSourcePath,
				  nodePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *child = group->getChild(i);
	    const SbString childPath =
		scene_child_summary_path(nodePath, child);
	    if (find_scene_bounds_summary_in_node(child, index,
						  ownerSourceIndex, ownerSourcePath, childPath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLSceneController::getSceneBoundsSummaryCount(void) const
{
    return this->getSceneTreeSummaryCount();
}

SbBool
SoBRLSceneController::getSceneBoundsSummary(int index,
	BRLObolSceneBoundsSummary &summary) const
{
    summary = BRLObolSceneBoundsSummary();
    if (index < 0 || !this->root)
	return FALSE;

    if (index == 0) {
	scene_bounds_summary_fill(this->root, -1, "", "/", summary);
	return TRUE;
    }
    index--;

    if (!this->root->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    const SoGroup *group = static_cast<const SoGroup *>(this->root);
    int sourceIndex = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *child = group->getChild(i);
	int childOwnerIndex = -1;
	SbString childOwnerPath("");
	if (child &&
	    child->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    const SoBRLDatabaseSource *source =
		static_cast<const SoBRLDatabaseSource *>(child);
	    childOwnerIndex = sourceIndex;
	    childOwnerPath = source->path.getValue();
	    sourceIndex++;
	}

	const SbString childPath = scene_child_summary_path("/", child);
	if (find_scene_bounds_summary_in_node(child, index,
					      childOwnerIndex, childOwnerPath, childPath, summary))
	    return TRUE;
    }

    return FALSE;
}

unsigned int
SoBRLSceneController::getLastVisitedSourceCount(void) const
{
    return this->lastVisitedSourceCount;
}

unsigned int
SoBRLSceneController::getLastRealizedSourceCount(void) const
{
    return this->lastRealizedSourceCount;
}

unsigned int
SoBRLSceneController::getLastFailedSourceCount(void) const
{
    return this->lastFailedSourceCount;
}

const SbString &
SoBRLSceneController::getLastDiagnostics(void) const
{
    return this->lastDiagnostics;
}

void
SoBRLSceneController::advanceFrameRevision(void)
{
    if (this->mutationBatchDepth > 0) {
	this->mutationBatchFrameRevisionPending = TRUE;
	return;
    }

    this->frameRevision++;
    if (this->frameRevision == 0)
	this->frameRevision++;
}

void
SoBRLSceneController::advanceStructuralRevision(void)
{
    if (this->mutationBatchDepth > 0) {
	this->mutationBatchStructuralRevisionPending = TRUE;
	this->mutationBatchFrameRevisionPending = TRUE;
	return;
    }

    this->structuralRevision++;
    if (this->structuralRevision == 0)
	this->structuralRevision++;
    this->advanceFrameRevision();
}
