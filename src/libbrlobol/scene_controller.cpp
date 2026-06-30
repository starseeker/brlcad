/*                S C E N E _ C O N T R O L L E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/material_object.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/realize_action.h"
#include "brlobol/scene_controller.h"
#include "brlobol/scene_group.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/SbName.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>

#include <string.h>
#include <string>
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
database_source_path_equal(const SoBRLDatabaseSource *source,
	const char *path)
{
    if (!source || !path)
	return 0;
    const char *sourcePath = source->path.getValue().getString();
    if (strcmp(sourcePath, path) == 0)
	return 1;
    return strcmp(skip_leading_slash(sourcePath),
	    skip_leading_slash(path)) == 0;
}

static SoBRLDatabaseSource *
find_database_source_recursive(SoGroup *group, const char *sourcePath,
	SoGroup **parentOut = NULL, int *childIndexOut = NULL)
{
    if (parentOut)
	*parentOut = NULL;
    if (childIndexOut)
	*childIndexOut = -1;
    if (!group || !sourcePath)
	return NULL;

    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    SoBRLDatabaseSource *source =
		static_cast<SoBRLDatabaseSource *>(node);
	    if (database_source_path_equal(source, sourcePath)) {
		if (parentOut)
		    *parentOut = group;
		if (childIndexOut)
		    *childIndexOut = i;
		return source;
	    }
	    continue;
	}
	if (node->isOfType(SoGroup::getClassTypeId())) {
	    SoBRLDatabaseSource *found =
		find_database_source_recursive(
			static_cast<SoGroup *>(node), sourcePath,
			parentOut, childIndexOut);
	    if (found)
		return found;
	}
    }

    return NULL;
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
	    continue;
	}
	if (node->isOfType(SoGroup::getClassTypeId()))
	    count += count_database_sources_recursive(
		    static_cast<const SoGroup *>(node));
    }
    return count;
}

static SoBRLDatabaseSource *
database_source_at_recursive(SoGroup *group, int index, int &seen)
{
    if (!group || index < 0)
	return NULL;

    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    if (seen == index)
		return static_cast<SoBRLDatabaseSource *>(node);
	    seen++;
	    continue;
	}
	if (node->isOfType(SoGroup::getClassTypeId())) {
	    SoBRLDatabaseSource *found = database_source_at_recursive(
		    static_cast<SoGroup *>(node), index, seen);
	    if (found)
		return found;
	}
    }

    return NULL;
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

static SoGroup *
scene_group_find_path(SoNode *sceneRoot, const char *groupPath)
{
    if (!groupPath)
	return NULL;

    SoGroup *current = scene_root_group(sceneRoot);
    if (!current)
	return NULL;

    std::vector<std::string> components;
    scene_group_path_components(groupPath, components);
    for (size_t i = 0; i < components.size(); i++) {
	current = scene_group_find_child(current, components[i].c_str());
	if (!current)
	    return NULL;
    }

    return current;
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

static SoBRLDatabaseSource *
database_source_summary_at_recursive(SoGroup *group, int index, int &seen,
	int groupDepth, const SbString &groupPath, SbString &parentGroupPath,
	int &drawTreeDepth)
{
    if (!group || index < 0)
	return NULL;

    const SbString effectiveGroupPath =
	scene_group_summary_path(group, groupPath);
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node)
	    continue;

	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    if (seen == index) {
		parentGroupPath = effectiveGroupPath;
		drawTreeDepth = groupDepth + 1;
		return static_cast<SoBRLDatabaseSource *>(node);
	    }
	    seen++;
	    continue;
	}

	if (node->isOfType(SoGroup::getClassTypeId())) {
	    const SbString childPath =
		scene_child_summary_path(effectiveGroupPath, node);
	    const SbString childGroupPath =
		scene_group_summary_path(node, childPath);
	    SoBRLDatabaseSource *found =
		database_source_summary_at_recursive(
			static_cast<SoGroup *>(node), index, seen,
			groupDepth + 1, childGroupPath, parentGroupPath,
			drawTreeDepth);
	    if (found)
		return found;
	}
    }

    return NULL;
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

SoBRLSceneController::SoBRLSceneController(void) :
    root(NULL),
    structuralRevision(0),
    frameRevision(0),
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
    action.apply(this->root);
    this->lastVisitedSourceCount = action.getVisitedSourceCount();
    this->lastRealizedSourceCount = action.getRealizedSourceCount();
    this->lastFailedSourceCount = action.getFailedSourceCount();
    this->lastDiagnostics = action.getDiagnostics();
    if (this->lastRealizedSourceCount > 0 || this->lastFailedSourceCount > 0)
	this->advanceFrameRevision();
    return this->lastFailedSourceCount == 0;
}

SoGroup *
SoBRLSceneController::findGroup(const char *groupPath) const
{
    return scene_group_find_path(this->root, groupPath);
}

SoGroup *
SoBRLSceneController::ensureGroup(const char *groupPath)
{
    if (!groupPath)
	return NULL;

    SoGroup *current = scene_root_group(this->root);
    if (!current)
	return NULL;

    std::vector<std::string> components;
    scene_group_path_components(groupPath, components);
    SbBool created = FALSE;
    SbString currentPath("");
    for (size_t i = 0; i < components.size(); i++) {
	const SbString childPath =
	    scene_group_append_path(currentPath, components[i].c_str());
	SoGroup *child = scene_group_find_child(current,
		components[i].c_str());
	if (!child) {
	    SoBRLSceneGroup *newGroup = new SoBRLSceneGroup;
	    newGroup->setName(SbName(components[i].c_str()));
	    newGroup->groupPath = childPath;
	    current->addChild(newGroup);
	    child = newGroup;
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
    SoGroup *group = scene_group_find_path(this->root, groupPath);
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
    SoGroup *group = scene_group_find_path(this->root, groupPath);
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
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::appendChildToGroup(const char *groupPath,
	SoNode *child)
{
    if (!child)
	return -1;

    SoGroup *group = scene_group_find_path(this->root, groupPath);
    if (!group)
	return -1;
    if (scene_group_find_child_index(group, child) >= 0)
	return 0;

    group->addChild(child);
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::removeChildFromGroup(const char *groupPath,
	SoNode *child)
{
    if (!child)
	return -1;

    SoGroup *group = scene_group_find_path(this->root, groupPath);
    if (!group)
	return -1;

    const int childIndex = scene_group_find_child_index(group, child);
    if (childIndex < 0)
	return 0;

    group->removeChild(childIndex);
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::eraseGroupSubpath(const char *parentGroupPath,
	const char *subpath)
{
    SoGroup *parent = scene_group_find_path(this->root, parentGroupPath);
    if (!parent)
	return -1;

    std::vector<std::string> components;
    scene_group_path_components(subpath, components);
    if (components.empty())
	return 0;

    const int erased = scene_group_erase_nested_subpath(parent, components);
    if (erased > 0)
	this->advanceStructuralRevision();
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

    parent->removeChild(target);
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::clearGroup(const char *groupPath)
{
    SoGroup *group = scene_group_find_path(this->root, groupPath);
    if (!group)
	return -1;

    const int removed = group->getNumChildren();
    if (removed <= 0)
	return 0;

    group->removeAllChildren();
    this->advanceStructuralRevision();
    return removed;
}

int
SoBRLSceneController::getGroupChildCount(const char *groupPath) const
{
    const SoGroup *group = scene_group_find_path_const(this->root,
	    groupPath);
    if (!group)
	return -1;
    return group->getNumChildren();
}

int
SoBRLSceneController::getGroupDatabaseSourceCount(
	const char *groupPath) const
{
    const SoGroup *group = scene_group_find_path_const(this->root,
	    groupPath);
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

    SoGroup *targetGroup = scene_group_find_path(this->root, groupPath);
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

    SoGroup *group = static_cast<SoGroup *>(this->root);
    int seen = 0;
    return database_source_at_recursive(group, index, seen);
}

int
SoBRLSceneController::getDatabaseSourceCount(void) const
{
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return 0;

    SoGroup *group = static_cast<SoGroup *>(this->root);
    return count_database_sources_recursive(group);
}

SoBRLDatabaseSource *
SoBRLSceneController::findDatabaseSource(const char *sourcePath) const
{
    if (!sourcePath || !sourcePath[0] || !this->root ||
	    !this->root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    SoGroup *group = static_cast<SoGroup *>(this->root);
    return find_database_source_recursive(group, sourcePath);
}

int
SoBRLSceneController::replaceDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int drawMode,
	uint32_t sourceRevision)
{
    if (!sourcePath || !sourcePath[0])
	return -1;
    if (!database)
	return this->removeDatabaseSource(sourcePath);
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(this->root);
    int childIndex = -1;
    SoBRLDatabaseSource *source =
	find_database_source_recursive(group, sourcePath, NULL,
		&childIndex);
    if (!source)
	source = new SoBRLDatabaseSource;

    if (drawMode != SoBRLDatabaseSource::SHADED)
	drawMode = SoBRLDatabaseSource::WIREFRAME;

    if (sourceRevision == 0)
	sourceRevision = source->sourceRevision.getValue() + 1;

    source->configureDatabaseSource(sourcePath, database, drawMode,
	    sourceRevision);

    if (childIndex < 0)
	group->addChild(source);
    if (childIndex < 0)
	this->advanceStructuralRevision();
    else
	this->advanceFrameRevision();

    return 1;
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
    SoBRLDatabaseSource *source = this->findDatabaseSource(sourcePath);
    if (!source)
	return -1;

    int changed = 0;
    if (sourceRevisionValid &&
	    source->sourceRevision.getValue() != sourceRevision) {
	source->sourceRevision = sourceRevision;
	changed = 1;
    }
    if (source->inputsRevision.getValue() != inputsRevision) {
	source->inputsRevision = inputsRevision;
	changed = 1;
    }
    if (source->visible.getValue() != visible) {
	source->visible = visible;
	changed = 1;
    }
    if (source->highlighted.getValue() != highlighted) {
	source->highlighted = highlighted;
	changed = 1;
    }
    if (source->lineStyle.getValue() != lineStyle) {
	source->lineStyle = lineStyle;
	changed = 1;
    }
    if (source->lineWidth.getValue() != lineWidth) {
	source->lineWidth = lineWidth;
	changed = 1;
    }
    if (scene_group_float_different(source->transparency.getValue(),
	    transparency)) {
	source->transparency = transparency;
	changed = 1;
    }
    if (source->colorOverride.getValue() != colorOverride) {
	source->colorOverride = colorOverride;
	changed = 1;
    }
    if (colorOverride &&
	    !scene_group_color_equal(source->color.getValue(), color)) {
	source->color = color;
	changed = 1;
    }
    if (source->materialColorValid.getValue() != materialColorValid) {
	source->materialColorValid = materialColorValid;
	changed = 1;
    }
    if (materialColorValid &&
	    !scene_group_color_equal(source->materialColor.getValue(),
		materialColor)) {
	source->materialColor = materialColor;
	changed = 1;
    }
    if (source->materialRevision.getValue() != materialRevision) {
	source->materialRevision = materialRevision;
	changed = 1;
    }

    if (changed > 0)
	this->advanceFrameRevision();
    return changed;
}

int
SoBRLSceneController::setDatabaseSourceDrawMode(const char *sourcePath,
	int drawMode)
{
    SoBRLDatabaseSource *source = this->findDatabaseSource(sourcePath);
    if (!source)
	return -1;

    if (drawMode != SoBRLDatabaseSource::SHADED)
	drawMode = SoBRLDatabaseSource::WIREFRAME;
    if (source->drawMode.getValue() == drawMode)
	return 0;

    source->drawMode = drawMode;
    this->advanceFrameRevision();
    return 1;
}

int
SoBRLSceneController::setDatabaseSourceMaterialPolicy(const char *sourcePath,
	int materialPolicy)
{
    SoBRLDatabaseSource *source = this->findDatabaseSource(sourcePath);
    if (!source)
	return -1;

    if (materialPolicy != SoBRLDatabaseSource::MATERIAL_DATABASE)
	materialPolicy = SoBRLDatabaseSource::MATERIAL_INHERIT;
    if (source->materialPolicy.getValue() == materialPolicy)
	return 0;

    source->materialPolicy = materialPolicy;
    this->advanceFrameRevision();
    return 1;
}

int
SoBRLSceneController::markDatabaseSourceStale(const char *sourcePath,
	uint32_t staleReason)
{
    SoBRLDatabaseSource *source = this->findDatabaseSource(sourcePath);
    if (!source)
	return -1;

    if (!staleReason)
	staleReason = SoBRLDatabaseSource::STALE_SOURCE;

    const uint32_t nextReason =
	source->staleReason.getValue() | staleReason;
    const int changed = !source->stale.getValue() ||
	source->staleReason.getValue() != nextReason ||
	source->realizationStatus.getValue() != SoBRLDatabaseSource::UNREALIZED ||
	source->realizationDiagnostic.getValue().getLength() > 0;
    if (!changed)
	return 0;

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
    SoBRLDatabaseSource *source = this->findDatabaseSource(sourcePath);
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
    SoBRLDatabaseSource *source = this->findDatabaseSource(sourcePath);
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
    SoBRLDatabaseSource *source = this->findDatabaseSource(sourcePath);
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
    if (!sourcePath || !sourcePath[0] || !groupPath)
	return -1;
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *rootGroup = static_cast<SoGroup *>(this->root);
    SoGroup *sourceParent = NULL;
    int sourceIndex = -1;
    SoBRLDatabaseSource *source = find_database_source_recursive(rootGroup,
	    sourcePath, &sourceParent, &sourceIndex);
    if (!source || !sourceParent || sourceIndex < 0)
	return 0;

    SoGroup *targetGroup = this->ensureGroup(groupPath);
    if (!targetGroup)
	return -1;
    if (targetGroup == sourceParent)
	return 0;

    source->ref();
    sourceParent->removeChild(sourceIndex);
    targetGroup->addChild(source);
    source->unref();
    this->advanceStructuralRevision();
    return 1;
}

int
SoBRLSceneController::removeDatabaseSource(const char *sourcePath)
{
    if (!sourcePath || !sourcePath[0])
	return 0;
    if (!this->root || !this->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(this->root);
    SoGroup *sourceParent = NULL;
    int childIndex = -1;
    (void)find_database_source_recursive(group, sourcePath, &sourceParent,
	    &childIndex);
    if (!sourceParent || childIndex < 0)
	return 0;

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
    if (removed > 0)
	this->advanceStructuralRevision();
    return removed;
}

SbBool
SoBRLSceneController::getDatabaseSourceSummary(int index,
	BRLObolDatabaseSourceSummary &summary) const
{
    summary = BRLObolDatabaseSourceSummary();
    if (index < 0 || !this->root ||
	    !this->root->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    SoGroup *group = static_cast<SoGroup *>(this->root);
    int seen = 0;
    SbString parentGroupPath("");
    int drawTreeDepth = 0;
    SoBRLDatabaseSource *source = database_source_summary_at_recursive(
	    group, index, seen, 0, "/", parentGroupPath, drawTreeDepth);
    if (!source)
	return FALSE;

    if (!source->getSummary(summary))
	return FALSE;

    summary.hasParent = TRUE;
    summary.parentGroupPath = parentGroupPath;
    summary.drawTreeDepth = drawTreeDepth;
    return TRUE;
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
	if (summary.ownerSourcePath.getLength() == 0)
	    summary.ownerSourcePath = source->path.getValue();
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
    if (!nodePath)
	return 0;

    const char *stored = summaryPath.getString();
    if (strcmp(stored, nodePath) == 0)
	return 1;
    return strcmp(skip_leading_slash(stored),
	    skip_leading_slash(nodePath)) == 0;
}

SbBool
SoBRLSceneController::getSceneTreeSummaryForPath(const char *nodePath,
	BRLObolSceneTreeSummary &summary) const
{
    summary = BRLObolSceneTreeSummary();
    if (!this->root || !nodePath)
	return FALSE;

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
	const SoNode *publicChild = scene_public_realized_child_node(child);
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
	if (summary.path.getLength() == 0)
	    summary.path = source->path.getValue();
	return TRUE;
    }

    const SoGroup *group = this->findGroup(nodePath);
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
scene_bounds_for_node(const SoNode *node, SbBox3f &bounds,
	SbBool includeOverlays)
{
    bounds.makeEmpty();
    if (!node)
	return FALSE;

    if (!includeOverlays && scene_node_is_overlay_intent(node))
	return FALSE;

    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return scene_bounds_for_vlist_shape(
		static_cast<const SoBRLVListShape *>(node), bounds);

    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return scene_bounds_for_mesh_shape(
		static_cast<const SoBRLMeshShape *>(node), bounds);

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	BRLObolSceneBoundsSummary sourceBounds;
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	if (source->getRealizedBoundsSummary(0, sourceBounds) &&
		sourceBounds.boundsValid) {
	    bounds = sourceBounds.bounds;
	    return TRUE;
	}
	return FALSE;
    }

    SbBool valid = FALSE;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SbBox3f childBounds;
	    if (scene_bounds_for_node(group->getChild(i), childBounds,
		    includeOverlays)) {
		bounds.extendBy(childBounds);
		valid = TRUE;
	    }
	}
    }

    return valid;
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
	node = this->findGroup(path);
    if (!node)
	node = this->findDatabaseSource(path);
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
    this->frameRevision++;
    if (this->frameRevision == 0)
	this->frameRevision++;
}

void
SoBRLSceneController::advanceStructuralRevision(void)
{
    this->structuralRevision++;
    if (this->structuralRevision == 0)
	this->structuralRevision++;
    this->advanceFrameRevision();
}
