/*       D R A W _ O B O L _ S C E N E _ R E C O R D S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file draw_obol_scene_records.cpp
 *
 * Database-source and group record enumeration, visibility, selection, and
 * display-state operations for the Obol scene bridge.
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BDrawCache.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSceneController.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BVListShape.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "ged.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "vmath.h"

#include "./draw_obol_bridge_private.hpp"
#include "./draw_obol_overlay_private.hpp"
#include "./ged_bobol_private.hpp"
#include "./ged_draw_private.h"
#include "./ged_private.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <algorithm>
#include <math.h>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

static bool
ged_obol_group_path_is_overlay(BObolSceneController *scene,
			       const char *group_path)
{
    if (!scene || !group_path || !group_path[0] ||
	BU_STR_EQUAL(group_path, "/"))
	return false;

    std::string current = ged_obol_skip_leading_slash(group_path);
    while (!current.empty()) {
	SoGroup *group = scene->findGroup(current.c_str());
	if (group && group->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	    SoBRLSceneGroup *scene_group =
		static_cast<SoBRLSceneGroup *>(group);
	    if (scene_group->overlayIntent.getValue())
		return true;
	}

	const size_t slash = current.rfind('/');
	if (slash == std::string::npos)
	    break;
	current = current.substr(0, slash);
    }

    return false;
}

extern "C" int
ged_draw_obol_database_source_count(
    struct ged *gedp,
    int skip_overlay_groups,
    size_t *out)
{
    if (out)
	*out = 0;
    if (!out || !gedp)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int source_count = scene->getDatabaseSourceCount();
    if (source_count <= 0)
	return 1;

    if (!skip_overlay_groups) {
	*out = (size_t)source_count;
	return 1;
    }

    size_t count = 0;
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;
	if (ged_obol_group_path_is_overlay(scene,
					   summary.parentGroupPath.getString()))
	    continue;
	count++;
    }

    *out = count;
    return 1;
}

extern "C" int
ged_draw_obol_database_source_paths_foreach(
    struct ged *gedp,
    int skip_overlay_groups,
    ged_draw_obol_database_source_path_cb cb,
    void *userdata)
{
    if (!gedp || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    if (source_count <= 0)
	return 1;

    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	if (skip_overlay_groups &&
	    ged_obol_group_path_is_overlay(scene,
					   summary.parentGroupPath.getString()))
	    continue;

	const char *source_path = summary.path.getString();
	if (source_path && source_path[0] && !(*cb)(gedp, source_path,
		userdata))
	    return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_database_source_records_foreach(
    struct ged *gedp,
    int skip_overlay_groups,
    ged_draw_obol_database_source_record_cb cb,
    void *userdata)
{
    if (!gedp || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    if (source_count <= 0)
	return 1;

    std::vector<BObolDatabaseSourceSummary> summaries;
    summaries.reserve(static_cast<size_t>(source_count));
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	if (skip_overlay_groups &&
	    ged_obol_group_path_is_overlay(scene,
					   summary.parentGroupPath.getString()))
	    continue;

	summaries.push_back(summary);
    }

    std::stable_sort(summaries.begin(), summaries.end(),
	[scene](const BObolDatabaseSourceSummary &a,
		const BObolDatabaseSourceSummary &b) {
	const auto active_in_parent = [scene](
	    const BObolDatabaseSourceSummary &summary) -> int {
	    const char *group_path = summary.parentGroupPath.getString();
	    SoGroup *group = (group_path && group_path[0]) ?
			     scene->findGroup(group_path) : NULL;
	    if (!group ||
		!group->isOfType(SoBRLSceneGroup::getClassTypeId()))
		return 0;
	    SoBRLSceneGroup *scene_group =
		static_cast<SoBRLSceneGroup *>(group);
	    const int group_mode =
		ged_obol_lod_draw_mode_to_ged(
		    scene_group->drawMode.getValue());
	    return group_mode ==
		   ged_obol_database_source_summary_ged_mode(summary) ? 1 : 0;
	};
	const int a_active = active_in_parent(a);
	const int b_active = active_in_parent(b);
	if (a_active != b_active)
	    return a_active > b_active;
	return false;
    });

    for (const BObolDatabaseSourceSummary &summary : summaries) {
	struct ged_draw_obol_database_source_record record;
	if (!ged_obol_database_source_record_from_summary(&record, summary))
	    continue;
	if (record.database_path && record.database_path[0] &&
	    !(*cb)(gedp, &record, userdata))
	    return 0;

	SoBRLDatabaseSource *source = scene->findDatabaseSourceInstance(
	    summary.instanceKey.getString());
	if (!source || !source->hasCompactInstanceIndex())
	    continue;
	for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary occurrence;
	    if (!source->getCompactInstanceHandle(i, handle) ||
		!source->getCompactInstanceSummary(handle, occurrence) ||
		occurrence.path.getLength() == 0)
		continue;
	    struct ged_draw_obol_database_source_record occurrence_record =
		record;
	    occurrence_record.database_path = occurrence.path.getString();
	    occurrence_record.visible = record.visible && occurrence.visible;
	    occurrence_record.selected = occurrence.selected ? 1 : 0;
	    occurrence_record.highlighted = occurrence.highlighted ? 1 : 0;
	    if (!(*cb)(gedp, &occurrence_record, userdata))
		return 0;
	}
    }

    return 1;
}

extern "C" int
ged_draw_obol_visible_database_source_records_foreach_fast(
    struct ged *gedp,
    ged_draw_obol_database_source_record_cb cb,
    void *userdata)
{
    if (!gedp || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (!source || !source->visible.getValue())
	    continue;

	const char *path = source->path.getValue().getString();
	const char *instance_key = source->instanceKey.getValue().getString();
	if (!path || !path[0])
	    continue;
	if (!instance_key || !instance_key[0])
	    instance_key = path;

	struct ged_draw_obol_database_source_record record;
	memset(&record, 0, sizeof(record));
	record.valid = 1;
	record.database_path = path;
	record.source_path = path;
	record.instance_key = instance_key;
	record.visible = 1;
	const int representation_mode = source->representationMode.getValue();
	record.draw_mode = representation_mode >= 0 ? representation_mode :
	    ged_obol_database_draw_mode_to_ged(source->drawMode.getValue());
	if (!(*cb)(gedp, &record, userdata))
	    return 0;
    }

    return 1;
}

static SbBool
ged_obol_shape_node_bool_field(SoNode *node,
			       const char *field_name,
			       SbBool fallback)
{
    if (!node || !field_name)
	return fallback;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	if (BU_STR_EQUAL(field_name, "databaseIntent"))
	    return shape->databaseIntent.getValue();
	if (BU_STR_EQUAL(field_name, "overlayIntent"))
	    return shape->overlayIntent.getValue();
	if (BU_STR_EQUAL(field_name, "localSource"))
	    return shape->localSource.getValue();
	if (BU_STR_EQUAL(field_name, "sharedSource"))
	    return shape->sharedSource.getValue();
	if (BU_STR_EQUAL(field_name, "nonDatabaseSource"))
	    return shape->nonDatabaseSource.getValue();
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	if (BU_STR_EQUAL(field_name, "databaseIntent"))
	    return shape->databaseIntent.getValue();
	if (BU_STR_EQUAL(field_name, "overlayIntent"))
	    return shape->overlayIntent.getValue();
	if (BU_STR_EQUAL(field_name, "localSource"))
	    return shape->localSource.getValue();
	if (BU_STR_EQUAL(field_name, "sharedSource"))
	    return shape->sharedSource.getValue();
	if (BU_STR_EQUAL(field_name, "nonDatabaseSource"))
	    return shape->nonDatabaseSource.getValue();
    }

    return fallback;
}

const char *
ged_obol_shape_node_record_role(SoNode *node)
{
    if (!node)
	return "";

    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<SoBRLVListShape *>(node)->
	       recordRole.getValue().getString();

    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<SoBRLMeshShape *>(node)->
	       recordRole.getValue().getString();

    return "";
}

static int
ged_obol_shape_node_is_database_realization(SoNode *node)
{
    if (!node)
	return 0;

    const char *role = ged_obol_shape_node_record_role(node);
    return ged_obol_shape_node_bool_field(node, "databaseIntent", FALSE) &&
	   !ged_obol_shape_node_bool_field(node, "localSource", FALSE) &&
	   !ged_obol_shape_node_bool_field(node, "sharedSource", FALSE) &&
	   !ged_obol_shape_node_bool_field(node, "nonDatabaseSource", FALSE) &&
	   role && BU_STR_EQUAL(role, "database");
}

static int
ged_obol_shape_node_is_auxiliary_record(SoNode *node)
{
    const char *role = ged_obol_shape_node_record_role(node);
    return role && BU_STR_EQUAL(role, "auxiliary");
}

static int
ged_draw_obol_shape_paths_foreach_node(
    struct ged *gedp,
    SoNode *node,
    int skip_overlay_groups,
    ged_draw_obol_database_source_path_cb cb,
    void *userdata)
{
    if (!node)
	return 1;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return 1;

    if (skip_overlay_groups &&
	node->isOfType(SoBRLSceneGroup::getClassTypeId()) &&
	static_cast<SoBRLSceneGroup *>(node)->overlayIntent.getValue())
	return 1;

    if (node->isOfType(SoBRLVListShape::getClassTypeId()) ||
	node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	if (ged_obol_shape_node_is_database_realization(node) ||
	    ged_obol_shape_node_is_auxiliary_record(node))
	    return 1;
	if (skip_overlay_groups &&
	    ged_obol_shape_node_bool_field(node, "overlayIntent", FALSE))
	    return 1;

	const char *shape_path = NULL;
	if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	    shape_path = static_cast<SoBRLVListShape *>(node)->
			 sourcePath.getValue().getString();
	else
	    shape_path = static_cast<SoBRLMeshShape *>(node)->
			 sourcePath.getValue().getString();

	if (shape_path && shape_path[0] && !(*cb)(gedp, shape_path, userdata))
	    return 0;
	return 1;
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return 1;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (!ged_draw_obol_shape_paths_foreach_node(gedp, group->getChild(i),
	    skip_overlay_groups, cb, userdata))
	    return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_shape_paths_foreach(
    struct ged *gedp,
    int skip_overlay_groups,
    ged_draw_obol_database_source_path_cb cb,
    void *userdata)
{
    if (!gedp || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    return ged_draw_obol_shape_paths_foreach_node(gedp, scene->getSceneRoot(),
	    skip_overlay_groups, cb, userdata);
}

extern "C" int
ged_draw_obol_group_database_source_paths_foreach(
    struct ged *gedp,
    const char *group_path,
    ged_draw_obol_database_source_path_cb cb,
    void *userdata)
{
    if (!gedp || !group_path || !group_path[0] || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const std::string target_group =
	ged_obol_group_path_from_record_path(group_path);
    if (target_group.empty())
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	const std::string owner_group_path =
	    ged_obol_database_source_owner_group_path_from_summary(summary);
	if (!ged_obol_path_equal(owner_group_path.c_str(),
				 target_group.c_str()) &&
	    !ged_obol_path_has_prefix(owner_group_path.c_str(),
				      target_group.c_str()))
	    continue;

	const char *source_path = summary.path.getString();
	if (source_path && source_path[0] && !(*cb)(gedp, source_path,
		userdata))
	    return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_group_database_source_records_foreach(
    struct ged *gedp,
    const char *group_path,
    ged_draw_obol_database_source_record_cb cb,
    void *userdata)
{
    if (!gedp || !group_path || !group_path[0] || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const std::string target_group =
	ged_obol_group_path_from_record_path(group_path);
    if (target_group.empty())
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	const std::string owner_group_path =
	    ged_obol_database_source_owner_group_path_from_summary(summary);
	if (!ged_obol_path_equal(owner_group_path.c_str(),
				 target_group.c_str()) &&
	    !ged_obol_path_has_prefix(owner_group_path.c_str(),
				      target_group.c_str()))
	    continue;

	struct ged_draw_obol_database_source_record record;
	if (!ged_obol_database_source_record_from_summary(&record, summary))
	    continue;
	if (record.database_path && record.database_path[0] &&
	    !(*cb)(gedp, &record, userdata))
	    return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_database_source_owner_group_path_for_path(
    struct ged *gedp,
    const char *path,
    struct bu_vls *out)
{
    if (out)
	bu_vls_trunc(out, 0);
    if (!gedp || !path || !path[0] || !out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    BObolDatabaseSourceSummary summary;
    if (!ged_obol_database_source_controller_summary_for_path(scene, path,
	    summary))
	return 0;

    const std::string owner_group_path =
	ged_obol_database_source_owner_group_path_from_summary(summary);
    if (owner_group_path.empty())
	return 0;

    bu_vls_strcpy(out, owner_group_path.c_str());
    return 1;
}

extern "C" int
ged_draw_obol_database_source_remove_for_path(
    struct ged *gedp,
    const char *path)
{
    if (!gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> paths;
    ged_obol_append_unique_path(paths, path);
    return ged_obol_remove_paths(paths, NULL, scene);
}

extern "C" int
ged_draw_obol_database_sources_remove_for_path_prefix(
    struct ged *gedp,
    const char *path_prefix)
{
    if (!gedp || !path_prefix || !path_prefix[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> paths;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (!source)
	    continue;
	const char *source_path = source->path.getValue().getString();
	if (ged_obol_path_has_prefix(source_path, path_prefix))
	    ged_obol_append_unique_path(paths, source_path);
    }
    return ged_obol_remove_paths(paths, NULL, scene);
}

extern "C" int
ged_draw_obol_database_sources_remove_for_path_prefix_in_scope(
    struct ged *gedp,
    const char *path_prefix,
    struct ged_view_context *view_ctx)
{
    return ged_draw_obol_database_sources_remove_for_path_prefix_in_scope_mode(
	       gedp, path_prefix, view_ctx, -1);
}

extern "C" int
ged_draw_obol_database_sources_remove_for_path_prefix_in_scope_mode(
    struct ged *gedp,
    const char *path_prefix,
    struct ged_view_context *view_ctx,
    int mode)
{
    if (!gedp || !path_prefix || !path_prefix[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	const char *source_path = summary.path.getString();
	if (ged_obol_path_has_prefix(source_path, path_prefix))
	    ged_obol_append_database_source_instance_key(instance_keys,
		    summary);
    }
    return ged_obol_remove_instance_keys(instance_keys, scene);
}

extern "C" int
ged_draw_obol_active_database_sources_remove_for_path_prefix(
    struct ged *gedp,
    const char *path_prefix)
{
    return ged_draw_obol_active_database_sources_remove_for_path_prefix_mode(
	       gedp, path_prefix, -1);
}

extern "C" int
ged_draw_obol_active_database_sources_remove_for_path_prefix_mode(
    struct ged *gedp,
    const char *path_prefix,
    int mode)
{
    if (!gedp || !path_prefix || !path_prefix[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	const char *source_path = summary.path.getString();
	if (ged_obol_path_has_prefix(source_path, path_prefix))
	    ged_obol_append_database_source_instance_key(instance_keys,
		    summary);
    }
    return ged_obol_remove_instance_keys(instance_keys, scene);
}

extern "C" int
ged_draw_obol_database_sources_remove_for_component_name(
    struct ged *gedp,
    const char *name,
    int nonroot_only,
    int mode)
{
    if (!gedp || !name || !name[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> targets;
    ged_obol_append_unique_path(targets, name);

    std::vector<std::string> paths;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;
	if (!ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	for (const std::string &target : targets) {
	    if (ged_obol_path_has_component_name(summary.path.getString(),
						 target.c_str(), nonroot_only ? 1 : 0)) {
		ged_obol_append_unique_path(paths, summary.path.getString());
		break;
	    }
	}
    }

    return ged_obol_remove_paths(paths, NULL, scene, mode);
}

extern "C" int
ged_draw_obol_database_sources_clear(struct ged *gedp)
{
    if (!gedp)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    int changed = (scene->clearDatabaseSources() > 0) ? 1 : 0;
    if (ged_obol_prune_empty_groups(scene))
	changed = 1;
    if (changed)
	scene->realizePending();
    return changed;
}

extern "C" int
ged_draw_obol_database_sources_clear_in_scope(struct ged *gedp,
	struct ged_view_context *view_ctx)
{
    if (!gedp)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    int changed = ged_obol_clear_database_sources_in_scope(scene, view_ctx);
    if (changed)
	scene->realizePending();
    return changed;
}

int
ged_obol_scene_clear_controller(BObolSceneController *scene)
{
    if (!scene)
	return 0;

    int changed = (scene->clearDatabaseSources() > 0) ? 1 : 0;

    std::vector<std::string> group_paths;
    const int summary_count = scene->getSceneTreeSummaryCount();
    for (int i = 0; i < summary_count; i++) {
	BObolSceneTreeSummary tree_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !tree_summary.valid ||
	    tree_summary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	    tree_summary.path.getLength() == 0 ||
	    BU_STR_EQUAL(tree_summary.path.getString(), "/"))
	    continue;
	group_paths.push_back(tree_summary.path.getString());
    }

    std::sort(group_paths.begin(), group_paths.end(),
    [](const std::string &a, const std::string &b) {
	return a.size() > b.size();
    });
    for (const std::string &path : group_paths) {
	if (scene->removeGroup(path.c_str()) > 0)
	    changed = 1;
    }

    if (changed)
	scene->realizePending();
    return changed;
}

extern "C" int
ged_draw_obol_scene_clear(struct ged *gedp)
{
    if (!gedp)
	return 0;

    return ged_obol_scene_clear_controller(
	       ged_draw_obol_scene_controller(gedp));
}

extern "C" int
ged_draw_obol_active_scene_clear(struct ged *gedp)
{
    if (!gedp)
	return 0;

    return ged_obol_scene_clear_controller(
	       ged_draw_obol_scene_controller(gedp));
}

extern "C" int
ged_draw_obol_groups_remove_for_component_name(
    struct ged *gedp,
    const char *name)
{
    if (!gedp || !name || !name[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> paths;
    const int summary_count = scene->getSceneTreeSummaryCount();
    for (int i = summary_count - 1; i >= 0; i--) {
	BObolSceneTreeSummary tree_summary;
	BObolSceneDisplaySummary display_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !scene->getSceneDisplaySummary(i, display_summary) ||
	    !tree_summary.valid ||
	    !display_summary.valid ||
	    tree_summary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	    tree_summary.path.getLength() == 0 ||
	    BU_STR_EQUAL(tree_summary.path.getString(), "/") ||
	    !display_summary.hasDrawIntent ||
	    !ged_obol_intent_is_ged_draw_group(
		display_summary.intentPath) ||
	    !ged_obol_path_has_component_name(
		tree_summary.path.getString(), name, 0))
	    continue;

	SoGroup *group = scene->findGroup(tree_summary.path.getString());
	if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    continue;
	SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
	if (scene_group->overlayIntent.getValue())
	    continue;

	ged_obol_append_unique_path(paths, tree_summary.path.getString());
    }

    int changed = 0;
    for (const std::string &path : paths) {
	if (scene->removeGroup(path.c_str()) > 0)
	    changed = 1;
    }
    if (changed)
	scene->realizePending();
    return changed;
}

extern "C" int
ged_draw_obol_group_remove_for_path(
    struct ged *gedp,
    const char *path)
{
    if (!gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    if (group_path.empty())
	return 0;

    const int removed = scene->removeGroup(group_path.c_str());
    if (removed > 0)
	scene->realizePending();
    return removed > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_clear_for_path(
    struct ged *gedp,
    const char *path)
{
    if (!gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    if (group_path.empty())
	return 0;

    const int cleared = scene->clearGroup(group_path.c_str());
    if (cleared > 0)
	scene->realizePending();
    return cleared > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_rename_for_path(
    struct ged *gedp,
    const char *path,
    const char *new_path)
{
    if (!gedp || !path || !path[0] || !new_path || !new_path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    const std::string target_group_path =
	ged_obol_group_path_from_record_path(new_path);
    if (group_path.empty() || target_group_path.empty())
	return 0;

    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    if (group_path == target_group_path)
	return 1;

    std::string old_parent, old_leaf, new_parent, new_leaf;
    if (!ged_obol_group_parent_leaf(group_path, old_parent, old_leaf) ||
	!ged_obol_group_parent_leaf(target_group_path, new_parent,
				    new_leaf) ||
	old_parent != new_parent)
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    int group_draw_mode = scene_group->drawMode.getValue();
    int fallback_draw_mode = scene_group->fallbackDrawMode.getValue();
    if (fallback_draw_mode == BOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BOBOL_LOD_DRAW_WIRE;
    SbBool overlay = scene_group->overlayIntent.getValue();
    uint32_t revalidation_revision =
	scene_group->revalidationRevision.getValue();

    if (scene->renameGroup(group_path.c_str(), new_leaf.c_str()) <= 0)
	return 0;

    group = scene->findGroup(target_group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    const std::string intent_path =
	ged_obol_group_intent_path(target_group_path.c_str());
    int changed = scene->setGroupDrawIntent(target_group_path.c_str(),
					    intent_path.c_str(), group_draw_mode, fallback_draw_mode, overlay,
					    revalidation_revision);
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_group_erase_subpath_for_path(
    struct ged *gedp,
    const char *parent_path,
    const char *subpath)
{
    if (!gedp || !parent_path || !parent_path[0] ||
	!subpath || !subpath[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path =
	ged_obol_group_path_from_record_path(parent_path);
    const std::string relative_path =
	ged_obol_group_path_from_record_path(subpath);
    if (group_path.empty() || relative_path.empty())
	return 0;

    const int erased = scene->eraseGroupSubpath(group_path.c_str(),
		       relative_path.c_str());
    if (erased > 0)
	scene->realizePending();
    return erased > 0 ? 1 : 0;
}


static int ged_obol_collect_view_database_sources(
    struct ged_view_context *view_ctx,
    BObolViewController *controller,
    void *userdata);

static std::set<SoBRLDatabaseSource *>
ged_obol_attached_database_sources(struct ged *gedp)
{
    std::set<SoBRLDatabaseSource *> sources;
    BObolSceneController *owned = ged_draw_obol_scene_controller(gedp);
    if (owned) {
	for (int i = 0; i < owned->getDatabaseSourceCount(); i++)
	    sources.insert(owned->getDatabaseSource(i));
    }
    ged_bobol_view_controllers_foreach(gedp,
	ged_obol_collect_view_database_sources, &sources);
    return sources;
}


extern "C" int
ged_draw_obol_source_visibility_frontier_set(
    struct ged *gedp,
    const char *root_path,
    struct ged_view_context *view_ctx,
    int mode,
    const char *const *paths,
    size_t path_count)
{
    if (!gedp || !root_path || !root_path[0])
	return 0;

    std::vector<SbString> frontier;
    frontier.reserve(path_count);
    for (size_t i = 0; paths && i < path_count; i++) {
	if (paths[i] && paths[i][0])
	    frontier.push_back(SbString(paths[i]));
    }

    int changed = 0;
    const std::set<SoBRLDatabaseSource *> sources =
	ged_obol_attached_database_sources(gedp);
    for (SoBRLDatabaseSource *source : sources) {
	BObolDatabaseSourceSummary summary;
	if (!source || !source->getSummary(summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	const char *source_path = summary.path.getString();
	if (!ged_obol_path_equal(source_path, root_path) &&
	    ged_obol_semantic_path_string(source_path) !=
		ged_obol_semantic_path_string(root_path))
	    continue;
	if (source->setCompactInstanceVisibilityFrontier(frontier) > 0)
	    changed++;
    }
    return changed;
}


extern "C" int
ged_draw_obol_source_visibility_overrides_set(
    struct ged *gedp,
    const char *root_path,
    struct ged_view_context *view_ctx,
    int mode,
    const char *const *paths,
    const int *visible,
    size_t rule_count)
{
    if (!gedp || !root_path || !root_path[0] ||
	(rule_count && (!paths || !visible)))
	return 0;

    std::vector<SbString> rule_paths;
    std::vector<SbBool> rule_states;
    rule_paths.reserve(rule_count);
    rule_states.reserve(rule_count);
    for (size_t i = 0; i < rule_count; i++) {
	if (!paths[i] || !paths[i][0])
	    continue;
	rule_paths.push_back(SbString(paths[i]));
	rule_states.push_back(visible[i] ? TRUE : FALSE);
    }

    int changed = 0;
    const std::set<SoBRLDatabaseSource *> sources =
	ged_obol_attached_database_sources(gedp);
    for (SoBRLDatabaseSource *source : sources) {
	BObolDatabaseSourceSummary summary;
	if (!source || !source->getSummary(summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	const char *source_path = summary.path.getString();
	if (!ged_obol_path_equal(source_path, root_path) &&
	    ged_obol_semantic_path_string(source_path) !=
		ged_obol_semantic_path_string(root_path))
	    continue;
	if (source->setCompactInstanceVisibilityOverrides(
		rule_paths, rule_states) > 0)
	    changed++;
    }
    return changed;
}


extern "C" int
ged_draw_obol_source_visibility_frontier_clear(
    struct ged *gedp,
    const char *root_path,
    struct ged_view_context *view_ctx,
    int mode)
{
    if (!gedp || !root_path || !root_path[0])
	return 0;

    int changed = 0;
    const std::set<SoBRLDatabaseSource *> sources =
	ged_obol_attached_database_sources(gedp);
    for (SoBRLDatabaseSource *source : sources) {
	BObolDatabaseSourceSummary summary;
	if (!source || !source->getSummary(summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	const char *source_path = summary.path.getString();
	if (!ged_obol_path_equal(source_path, root_path) &&
	    ged_obol_semantic_path_string(source_path) !=
		ged_obol_semantic_path_string(root_path))
	    continue;
	if (source->clearCompactInstanceVisibilityFrontier() > 0)
	    changed++;
    }
    return changed;
}


static int
ged_obol_frontier_visibility_apply_cb(
    const struct ged_draw_frontier_visibility_record *record,
    void *userdata)
{
    struct ged *gedp = static_cast<struct ged *>(userdata);
    if (!gedp || !record || !record->root_path || !record->root_path[0])
	return 1;
    if (record->clear) {
	(void)ged_draw_obol_source_visibility_frontier_clear(gedp,
	    record->root_path, record->view, record->mode);
	return 1;
    }
    (void)ged_draw_obol_source_visibility_overrides_set(gedp,
	record->root_path, record->view, record->mode, record->paths,
	record->visible, record->rule_count);
    return 1;
}


void
ged_obol_frontier_visibility_changes_apply(struct ged *gedp)
{
    if (!gedp || !ged_draw_obol_scene_controller_attached(gedp))
	return;
    (void)ged_draw_frontier_visibility_changes_foreach(gedp,
	ged_obol_frontier_visibility_apply_cb, gedp);
}


void
ged_obol_frontier_visibility_snapshot_apply(struct ged *gedp)
{
    if (!gedp || !ged_draw_obol_scene_controller_attached(gedp))
	return;
    (void)ged_draw_frontier_visibility_snapshot_foreach(gedp,
	ged_obol_frontier_visibility_apply_cb, gedp);
}

static int
ged_bobol_database_source_update_display_for_path_impl(
    struct ged_bobol_publication_context *publication,
    struct ged *gedp,
    const char *path,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision)
{
    if ((color_valid && !color) ||
	(material_color_valid && !material_color))
	return 0;

    BObolSceneController *scene =
	(publication && publication->active) ? publication->scene :
	ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0,
	    publication);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;

	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    return 0;

	SbBool nextVisible = visible_valid ?
			     (visible ? TRUE : FALSE) : summary.visible;
	SbBool nextSelected = selected_valid ?
			      (selected ? TRUE : FALSE) : summary.selected;
	SbBool nextHighlighted = highlighted_valid ?
				 (highlighted ? TRUE : FALSE) :
				 summary.highlighted;
	int nextDrawMode = draw_mode_valid ?
			   ged_obol_database_draw_mode_from_ged(draw_mode) :
			   summary.drawMode;
	int nextLineStyle = line_style_valid ?
			    line_style : summary.lineStyle;
	int nextLineWidth = line_width_valid ?
			    line_width : summary.lineWidth;
	float nextTransparency = transparency_valid ?
				 static_cast<float>(transparency) :
				 summary.transparency;

	SbBool nextColorOverride = summary.colorOverride;
	SbColor nextColor = summary.color;
	if (color_valid) {
	    nextColorOverride = TRUE;
	    nextColor = ged_obol_color_from_rgb(color);
	}

	SbBool nextMaterialColorValid = summary.materialColorValid;
	SbColor nextMaterialColor = summary.materialColor;
	uint32_t nextMaterialRevision = summary.materialRevision;
	if (material_color_valid) {
	    nextMaterialColorValid = TRUE;
	    nextMaterialColor = ged_obol_color_from_rgb(material_color);
	    if (material_revision_valid) {
		nextMaterialRevision =
		    ged_obol_fold_revision(material_revision);
	    } else {
		nextMaterialRevision++;
		if (!nextMaterialRevision)
		    nextMaterialRevision = 1;
	    }
	}

	int draw_mode_changed = scene->setDatabaseSourceInstanceDrawMode(
				    source_instance_key.c_str(), nextDrawMode);
	if (draw_mode_changed < 0)
	    return 0;

	if (draw_mode_valid) {
	    const char *representation_key =
		summary.representationKey.getString();
	    if (!representation_key || !representation_key[0])
		representation_key = source_instance_key.c_str();
	    int representation_changed =
		scene->setDatabaseSourceInstanceRepresentation(
		    source_instance_key.c_str(),
		    representation_key,
		    ged_obol_database_representation_mode_from_ged(
			draw_mode));
	    if (representation_changed < 0)
		return 0;
	}

	int changed = scene->setDatabaseSourceInstanceState(
			  source_instance_key.c_str(),
			  FALSE,
			  summary.sourceRevision,
			  summary.inputsRevision,
			  nextVisible,
			  nextSelected,
			  nextHighlighted,
			  nextLineStyle,
			  nextLineWidth,
			  nextTransparency,
			  nextColorOverride,
			  nextColor,
			  nextMaterialColorValid,
			  nextMaterialColor,
			  nextMaterialRevision);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_update_display_for_path(
    struct ged *gedp,
    const char *path,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision)
{
    return ged_bobol_database_source_update_display_for_path_impl(NULL, gedp,
	path, visible_valid, visible, selected_valid, selected,
	highlighted_valid, highlighted, draw_mode_valid, draw_mode,
	line_style_valid, line_style, line_width_valid, line_width,
	transparency_valid, transparency, color_valid, color,
	material_color_valid, material_color, material_revision_valid,
	material_revision);
}

int
ged_bobol_database_source_update_display_for_path(
    struct ged_bobol_publication_context *publication,
    const char *path,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision)
{
    if (!publication || !publication->active || !publication->gedp)
	return 0;
    return ged_bobol_database_source_update_display_for_path_impl(publication,
	publication->gedp, path, visible_valid, visible, selected_valid,
	selected, highlighted_valid, highlighted, draw_mode_valid, draw_mode,
	line_style_valid, line_style, line_width_valid, line_width,
	transparency_valid, transparency, color_valid, color,
	material_color_valid, material_color, material_revision_valid,
	material_revision);
}

extern "C" int
ged_draw_obol_database_source_set_selected_for_instance_key(
    struct ged *gedp,
    const char *instance_key,
    int selected)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !instance_key || !instance_key[0])
	return 0;

    SoBRLDatabaseSource *source =
	scene->findDatabaseSourceInstance(instance_key);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    const int changed = scene->setDatabaseSourceInstanceState(
			    instance_key,
			    FALSE,
			    summary.sourceRevision,
			    summary.inputsRevision,
			    summary.visible,
			    selected ? TRUE : FALSE,
			    summary.highlighted,
			    summary.lineStyle,
			    summary.lineWidth,
			    summary.transparency,
			    summary.colorOverride,
			    summary.color,
			    summary.materialColorValid,
			    summary.materialColor,
			    summary.materialRevision);
    return changed > 0 ? 1 : 0;
}

static void
ged_obol_collect_database_sources(SoNode *node,
	std::set<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	sources.insert(static_cast<SoBRLDatabaseSource *>(node));
	return;
    }
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	ged_obol_collect_database_sources(group->getChild(i), sources);
}

static int
ged_obol_collect_view_database_sources(
    struct ged_view_context *UNUSED(view_ctx),
    BObolViewController *controller, void *userdata)
{
    std::set<SoBRLDatabaseSource *> *sources =
	static_cast<std::set<SoBRLDatabaseSource *> *>(userdata);
    if (controller && sources)
	ged_obol_collect_database_sources(controller->getRenderSceneRoot(),
	    *sources);
    return 1;
}

static bool
ged_obol_path_in_selected_set(const char *candidate,
	const char *const *paths, size_t path_count)
{
    if (!candidate || !candidate[0])
	return false;
    for (size_t i = 0; i < path_count; i++) {
	if (paths[i] && paths[i][0] &&
	    (ged_obol_path_has_prefix(candidate, paths[i]) ||
	     ged_obol_path_has_semantic_prefix(candidate, paths[i])))
	    return true;
    }
    return false;
}

extern "C" int
ged_draw_obol_database_sources_sync_selected_paths(
    struct ged *gedp,
    const char *const *paths,
    size_t path_count)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp)
	return 0;

    std::set<SoBRLDatabaseSource *> sources;
    BObolSceneController *owned = ged_draw_obol_scene_controller(gedp);
    if (owned) {
	for (int i = 0; i < owned->getDatabaseSourceCount(); i++)
	    sources.insert(owned->getDatabaseSource(i));
    }

    ged_bobol_view_controllers_foreach(gedp,
	ged_obol_collect_view_database_sources, &sources);
    int applied = 0;
    std::vector<SbString> compact_selected_paths;
    compact_selected_paths.reserve(path_count);
    for (size_t i = 0; i < path_count; i++) {
	if (paths[i] && paths[i][0])
	    compact_selected_paths.push_back(SbString(paths[i]));
    }
    for (SoBRLDatabaseSource *source : sources) {
	BObolDatabaseSourceSummary sourceSummary;
	if (!source || !source->getSummary(sourceSummary) ||
	    !sourceSummary.valid)
	    continue;
	/* Retain the selection rule before deferred compact geometry arrives.
	 * Index installation and streamed additions will then apply it directly,
	 * rather than performing a delayed full-scene catch-up during the next
	 * unrelated draw or erase operation. */
	applied += source->syncCompactInstanceSelectedPaths(
	    compact_selected_paths);

	/* Compact occurrences own their individual selection presentation.  Once
	 * they exist, also selecting their aggregate source styles the retired
	 * overview/root proxy a second time and can make a white box reappear over
	 * an otherwise correct selected model. */
	const SbBool sourceSelected = source->getCompactInstanceCount() > 0 ?
	    FALSE : (ged_obol_path_in_selected_set(
		sourceSummary.path.getString(), paths, path_count) ? TRUE : FALSE);
	if (sourceSummary.selected == sourceSelected)
	    continue;
	const int changed = source->setDisplayState(FALSE,
	    sourceSummary.sourceRevision, sourceSummary.inputsRevision,
	    sourceSummary.visible, sourceSelected, sourceSummary.highlighted,
	    sourceSummary.lineStyle, sourceSummary.lineWidth,
	    sourceSummary.transparency, sourceSummary.colorOverride,
	    sourceSummary.color, sourceSummary.materialColorValid,
	    sourceSummary.materialColor, sourceSummary.materialRevision);
	if (changed > 0)
	    applied++;
    }
    return applied;
}

struct ged_obol_selection_scene_collect_ctx {
    std::set<BObolSceneController *> *scenes;
};

static int
ged_obol_selection_scene_collect(struct ged_view_context *UNUSED(view_ctx),
    BObolViewController *controller, void *userdata)
{
    ged_obol_selection_scene_collect_ctx *ctx =
	static_cast<ged_obol_selection_scene_collect_ctx *>(userdata);
    if (ctx && ctx->scenes && controller && controller->getSceneController())
	ctx->scenes->insert(controller->getSceneController());
    return 1;
}

static void
ged_obol_selection_sources_for_path(BObolSceneController *scene,
    const char *path, std::set<SoBRLDatabaseSource *> &sources)
{
    if (!scene || !path || !path[0])
	return;

    std::string candidate(ged_obol_skip_leading_slash(path));
    while (!candidate.empty()) {
	SoBRLDatabaseSource *source =
	    scene->findDatabaseSource(candidate.c_str());
	if (source)
	    sources.insert(source);

	const int count =
	    scene->getDatabaseSourceInstanceCountForPath(candidate.c_str());
	for (int i = 0; i < count; i++) {
	    BObolDatabaseSourceSummary summary;
	    if (!scene->getDatabaseSourceInstanceSummaryForPath(
		    candidate.c_str(), i, summary) || !summary.valid ||
		summary.instanceKey.getLength() <= 0)
		continue;
	    source = scene->findDatabaseSourceInstance(
		summary.instanceKey.getString());
	    if (source)
		sources.insert(source);
	}

	const size_t slash = candidate.rfind('/');
	if (slash == std::string::npos)
	    break;
	candidate.erase(slash);
    }
}

extern "C" int
ged_draw_obol_database_sources_apply_selection_delta(
    struct ged *gedp,
    const char *const *added_paths,
    size_t added_count,
    const char *const *removed_paths,
    size_t removed_count,
    const char *const *selected_paths,
    size_t selected_count)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp)
	return 0;

    std::set<BObolSceneController *> scenes;
    BObolSceneController *owned = ged_draw_obol_scene_controller(gedp);
    if (owned)
	scenes.insert(owned);
    ged_obol_selection_scene_collect_ctx scene_ctx = {&scenes};
    ged_bobol_view_controllers_foreach(gedp,
	ged_obol_selection_scene_collect, &scene_ctx);

    std::vector<SbString> added;
    std::vector<SbString> removed;
    added.reserve(added_count);
    removed.reserve(removed_count);
    for (size_t i = 0; added_paths && i < added_count; i++) {
	if (added_paths[i] && added_paths[i][0])
	    added.push_back(SbString(added_paths[i]));
    }
    for (size_t i = 0; removed_paths && i < removed_count; i++) {
	if (removed_paths[i] && removed_paths[i][0])
	    removed.push_back(SbString(removed_paths[i]));
    }

    int applied = 0;
    for (BObolSceneController *scene : scenes) {
	std::set<SoBRLDatabaseSource *> sources;
	for (size_t i = 0; added_paths && i < added_count; i++)
	    ged_obol_selection_sources_for_path(scene, added_paths[i], sources);
	for (size_t i = 0; removed_paths && i < removed_count; i++)
	    ged_obol_selection_sources_for_path(scene, removed_paths[i], sources);

	for (SoBRLDatabaseSource *source : sources) {
	    BObolDatabaseSourceSummary summary;
	    if (!source || !source->getSummary(summary) || !summary.valid)
		continue;
	    applied += source->applyCompactInstanceSelectionDelta(added,
		removed);
	    const SbBool selected = source->getCompactInstanceCount() > 0 ?
		FALSE : (ged_obol_path_in_selected_set(summary.path.getString(),
		    selected_paths, selected_count) ? TRUE : FALSE);
	    if (summary.selected != selected) {
		const int changed = source->setDisplayState(FALSE,
		    summary.sourceRevision, summary.inputsRevision,
		    summary.visible, selected, summary.highlighted,
		    summary.lineStyle, summary.lineWidth,
		    summary.transparency, summary.colorOverride,
		    summary.color, summary.materialColorValid,
		    summary.materialColor, summary.materialRevision);
		if (changed > 0)
		    applied++;
	    }
	}
    }
    return applied;
}

template <typename ShapeT>
static int
ged_obol_shape_update_display_typed(
    BObolSceneController *scene,
    const char *path,
    ShapeT *shape,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision)
{
    if (!scene || !path || !path[0] || !shape)
	return 0;

    if (draw_mode_valid) {
	const int draw_changed = scene->setShapeDrawState(path,
				 ged_obol_lod_draw_mode_from_ged(draw_mode),
				 shape->databaseIntent.getValue(),
				 shape->overlayIntent.getValue(),
				 shape->hudIntent.getValue());
	if (draw_changed < 0)
	    return 0;
    }

    const SbBool next_visible = visible_valid ?
				(visible ? TRUE : FALSE) : shape->visible.getValue();
    const SbBool next_selected = selected_valid ?
				 (selected ? TRUE : FALSE) : shape->selected.getValue();
    const SbBool next_highlighted = highlighted_valid ?
				    (highlighted ? TRUE : FALSE) : shape->highlighted.getValue();
    const int next_line_style = line_style_valid ?
				line_style : shape->lineStyle.getValue();
    const int next_line_width = line_width_valid ?
				line_width : shape->lineWidth.getValue();
    const float next_transparency = transparency_valid ?
				    static_cast<float>(transparency) : shape->transparency.getValue();

    SbBool next_color_override = shape->colorOverride.getValue();
    SbColor next_color = shape->color.getValue();
    if (color_valid) {
	next_color_override = TRUE;
	next_color = ged_obol_color_from_rgb(color);
    }

    SbBool next_material_color_valid = shape->materialColorValid.getValue();
    SbColor next_material_color = shape->materialColor.getValue();
    uint32_t next_material_revision = shape->materialRevision.getValue();
    if (material_color_valid) {
	next_material_color_valid = TRUE;
	next_material_color = ged_obol_color_from_rgb(material_color);
	if (material_revision_valid) {
	    next_material_revision =
		ged_obol_fold_revision(material_revision);
	} else {
	    next_material_revision++;
	    if (!next_material_revision)
		next_material_revision = 1;
	}
    }

    const int changed = scene->setShapeDisplayState(path,
			next_visible,
			next_selected,
			next_highlighted,
			next_line_style,
			next_line_width,
			next_transparency,
			next_color_override,
			next_color,
			next_material_color_valid,
			next_material_color,
			next_material_revision);
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_shape_update_display_for_path(
    struct ged *gedp,
    const char *path,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision)
{
    if ((color_valid && !color) ||
	(material_color_valid && !material_color))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node)
	return 0;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	return ged_obol_shape_update_display_typed(scene, path,
		static_cast<SoBRLVListShape *>(node),
		visible_valid, visible,
		selected_valid, selected,
		highlighted_valid, highlighted,
		draw_mode_valid, draw_mode,
		line_style_valid, line_style,
		line_width_valid, line_width,
		transparency_valid, transparency,
		color_valid, color,
		material_color_valid, material_color,
		material_revision_valid, material_revision);
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	return ged_obol_shape_update_display_typed(scene, path,
		static_cast<SoBRLMeshShape *>(node),
		visible_valid, visible,
		selected_valid, selected,
		highlighted_valid, highlighted,
		draw_mode_valid, draw_mode,
		line_style_valid, line_style,
		line_width_valid, line_width,
		transparency_valid, transparency,
		color_valid, color,
		material_color_valid, material_color,
		material_revision_valid, material_revision);
    }

    return 0;
}

extern "C" int
ged_draw_obol_database_source_update_appearance_for_path(
    struct ged *gedp,
    const char *path,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_override_valid,
    int color_override,
    int color_valid,
    const unsigned char color[3])
{
    if (color_valid && !color)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    BObolDatabaseSourceDisplayPatch patch;
    patch.lineWidthValid = line_width_valid ? TRUE : FALSE;
    patch.lineWidth = line_width;
    patch.transparencyValid = transparency_valid ? TRUE : FALSE;
    patch.transparency = static_cast<float>(transparency);
    patch.colorOverrideValid = color_override_valid ? TRUE : FALSE;
    patch.colorOverride = color_override ? TRUE : FALSE;
    patch.colorValid = color_valid ? TRUE : FALSE;
    if (color_valid)
	patch.color = ged_obol_color_from_rgb(color);

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	int changed = scene->setDatabaseSourceInstanceDisplayPatch(
			  source_instance_key.c_str(), patch);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_group_update_display_for_path(
    struct ged *gedp,
    const char *path,
    int visible_valid,
    int visible)
{
    if (!gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    int changed = scene->setGroupDisplayState(group_path.c_str(),
		  visible_valid ? (visible ? TRUE : FALSE) :
		  scene_group->visible.getValue(),
		  scene_group->selected.getValue(),
		  scene_group->highlighted.getValue(),
		  scene_group->lineStyle.getValue(),
		  scene_group->lineWidth.getValue(),
		  scene_group->transparency.getValue(),
		  scene_group->colorOverride.getValue(),
		  scene_group->color.getValue(),
		  scene_group->materialColorValid.getValue(),
		  scene_group->materialColor.getValue(),
		  scene_group->materialRevision.getValue());
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_ensure_for_path(
    struct ged *gedp,
    const char *path,
    const char *intent_path,
    int mode,
    int overlay)
{
    if (!gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    if (group_path.empty())
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoGroup *group = scene->ensureGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    const int next_draw_mode = mode >= 0 ?
			       ged_obol_lod_draw_mode_from_ged(mode) :
			       scene_group->drawMode.getValue();
    int fallback_draw_mode = scene_group->fallbackDrawMode.getValue();
    if (fallback_draw_mode == BOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BOBOL_LOD_DRAW_WIRE;
    const SbBool next_overlay = overlay >= 0 ?
				(overlay ? TRUE : FALSE) :
				scene_group->overlayIntent.getValue();

    const std::string target_path =
	(intent_path && intent_path[0]) ?
	ged_obol_group_path_from_record_path(intent_path) : group_path;
    const std::string draw_intent_path =
	ged_obol_group_intent_path(
	    target_path.empty() ? group_path.c_str() :
	    target_path.c_str());

    const int changed = scene->setGroupDrawIntent(group_path.c_str(),
			draw_intent_path.c_str(),
			next_draw_mode,
			fallback_draw_mode,
			next_overlay,
			scene_group->revalidationRevision.getValue());
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_record_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_group_record_summary *out)
{
    if (!out || !gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    const char *record_path = ged_obol_group_record_path(scene_group);
    if (record_path && record_path[0])
	out->path = record_path;
    out->draw_mode = ged_obol_lod_draw_mode_to_ged(
			 scene_group->drawMode.getValue());
    out->transparency = scene_group->transparency.getValue();
    out->visible = scene_group->visible.getValue() ? 1 : 0;
    out->is_overlay = scene_group->overlayIntent.getValue() ? 1 : 0;
    return 1;
}


static int
ged_draw_obol_group_paths_foreach_node(
    struct ged *gedp,
    BObolSceneController *scene,
    SoNode *node,
    int skip_overlay_groups,
    ged_draw_obol_group_path_cb cb,
    void *userdata)
{
    if (!node)
	return 1;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return 1;

    if (node->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	SoBRLSceneGroup *group = static_cast<SoBRLSceneGroup *>(node);
	const char *group_path = group->groupPath.getValue().getString();
	if ((!skip_overlay_groups ||
	     !ged_obol_group_path_is_overlay(scene, group_path)) &&
	    group_path && group_path[0] &&
	    !BU_STR_EQUAL(group_path, "/")) {
	    const char *record_path = ged_obol_group_record_path(group);
	    if (record_path && record_path[0] &&
		!(*cb)(gedp, record_path, userdata))
		return 0;
	}
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return 1;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (!ged_draw_obol_group_paths_foreach_node(gedp, scene,
		group->getChild(i), skip_overlay_groups, cb, userdata))
	    return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_group_paths_foreach(
    struct ged *gedp,
    int skip_overlay_groups,
    ged_draw_obol_group_path_cb cb,
    void *userdata)
{
    if (!gedp || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    return ged_draw_obol_group_paths_foreach_node(gedp, scene,
	scene->getSceneRoot(), skip_overlay_groups, cb, userdata);
}


extern "C" int
ged_draw_obol_group_display_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_scene_display_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    int source_like_local_group = 0;
    for (int i = 0; i < scene_group->getNumChildren(); i++) {
	SoNode *child = scene_group->getChild(i);
	if (child &&
	    (child->isOfType(SoBRLVListShape::getClassTypeId()) ||
	     child->isOfType(SoBRLMeshShape::getClassTypeId())) &&
	    ged_obol_shape_node_bool_field(child, "databaseIntent", FALSE) &&
	    ged_obol_shape_node_bool_field(child, "localSource", FALSE)) {
	    source_like_local_group = 1;
	    break;
	}
    }

    out->valid = 1;
    out->is_database_source = source_like_local_group;
    out->has_draw_intent =
	scene_group->drawIntentValid.getValue() ? 1 : 0;
    out->intent_path = scene_group->drawIntentPath.getValue().getString();
    out->intent_draw_mode = ged_obol_lod_draw_mode_to_ged(
				scene_group->drawMode.getValue());
    out->visible = scene_group->visible.getValue() ? 1 : 0;
    out->selected = scene_group->selected.getValue() ? 1 : 0;
    out->highlighted = scene_group->highlighted.getValue() ? 1 : 0;
    out->line_style = scene_group->lineStyle.getValue();
    out->line_width = scene_group->lineWidth.getValue();
    out->transparency = scene_group->transparency.getValue();
    out->draw_mode = ged_obol_lod_draw_mode_to_ged(
			 scene_group->drawMode.getValue());
    out->material_valid =
	(scene_group->materialColorValid.getValue() ||
	 scene_group->colorOverride.getValue()) ? 1 : 0;
    if (scene_group->materialColorValid.getValue())
	ged_obol_rgb_from_color(scene_group->materialColor.getValue(),
				out->material_color);
    else if (scene_group->colorOverride.getValue())
	ged_obol_rgb_from_color(scene_group->color.getValue(),
				out->material_color);
    return 1;
}

extern "C" int
ged_draw_obol_group_shape_count_for_path(
    struct ged *gedp,
    const char *path,
    int *out)
{
    if (out)
	*out = 0;
    if (!out || !gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    int count = scene->getGroupDatabaseSourceCount(group_path.c_str());
    if (count < 0)
	count = 0;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	const std::string owner_group_path =
	    ged_obol_database_source_owner_group_path_from_summary(summary);
	if (ged_obol_path_equal(owner_group_path.c_str(),
				group_path.c_str()))
	    continue;
	if (ged_obol_path_has_prefix(owner_group_path.c_str(),
				     group_path.c_str()))
	    count++;
    }

    *out = count;
    return 1;
}


extern "C" int
ged_draw_obol_group_descendant_group_count_for_path(
    struct ged *gedp,
    const char *path,
    size_t *out)
{
    if (out)
	*out = 0;
    if (!out || !gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    const int count = scene->getGroupDescendantGroupCount(
			  group_path.c_str());
    if (count < 0)
	return 0;

    *out = static_cast<size_t>(count);
    return 1;
}


extern "C" int
ged_draw_obol_group_child_count_for_path(
    struct ged *gedp,
    const char *path,
    size_t *out)
{
    if (out)
	*out = 0;
    if (!out || !gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    const int count = scene->getGroupChildCount(group_path.c_str());
    if (count < 0)
	return 0;

    *out = static_cast<size_t>(count);
    return 1;
}


extern "C" int
ged_draw_obol_group_update_appearance_for_path(
    struct ged *gedp,
    const char *path,
    const struct ged_draw_appearance_settings *settings)
{
    if (!gedp || !path || !path[0] || !settings)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    const char *retained_intent =
	scene_group->drawIntentPath.getValue().getString();
    const std::string next_intent_path =
	(retained_intent && retained_intent[0]) ?
	std::string(retained_intent) :
	ged_obol_group_intent_path(group_path.c_str());
    int fallback_draw_mode = scene_group->fallbackDrawMode.getValue();
    if (fallback_draw_mode == BOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BOBOL_LOD_DRAW_WIRE;
    int intent_changed = scene->setGroupDrawIntent(group_path.c_str(),
			 next_intent_path.c_str(),
			 ged_obol_lod_draw_mode_from_ged(settings->draw_mode),
			 fallback_draw_mode,
			 scene_group->overlayIntent.getValue(),
			 scene_group->revalidationRevision.getValue());
    if (intent_changed < 0)
	return 0;

    const SbColor next_color = ged_obol_color_from_rgb(settings->color);
    int changed = scene->setGroupDisplayState(group_path.c_str(),
		  scene_group->visible.getValue(),
		  scene_group->selected.getValue(),
		  scene_group->highlighted.getValue(),
	  scene_group->lineStyle.getValue(),
	  settings->s_line_width,
	  ged_obol_transparency_from_appearance_opacity(settings->transparency),
		  settings->color_override ? TRUE : FALSE,
		  next_color,
		  scene_group->materialColorValid.getValue(),
		  scene_group->materialColor.getValue(),
		  scene_group->materialRevision.getValue());
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_appearance_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_appearance_settings *settings)
{
    if (!gedp || !path || !path[0] || !settings)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    struct ged_draw_appearance_settings next =
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
    next.draw_mode = ged_obol_lod_draw_mode_to_ged(
			 scene_group->drawMode.getValue());
    next.transparency = ged_obol_appearance_opacity_from_transparency(
	scene_group->transparency.getValue());
    next.color_override = scene_group->colorOverride.getValue() ? 1 : 0;
    ged_obol_rgb_from_color(scene_group->color.getValue(), next.color);
    next.s_line_width = scene_group->lineWidth.getValue();
    *settings = next;
    return 1;
}

extern "C" int
ged_draw_obol_group_update_draw_intent_for_path(
    struct ged *gedp,
    const char *path,
    const char *intent_path,
    int mode_valid,
    int mode,
    int overlay_valid,
    int overlay)
{
    if (!gedp)
	return 0;

    std::string group_path = ged_obol_group_path_from_record_path(path);
    std::string target_group_path =
	(intent_path && intent_path[0]) ?
	ged_obol_group_path_from_record_path(intent_path) : group_path;
    if (group_path.empty())
	group_path = target_group_path;
    if (target_group_path.empty())
	target_group_path = group_path;
    if (group_path.empty())
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group && target_group_path != group_path) {
	group = scene->findGroup(target_group_path.c_str());
	if (group)
	    group_path = target_group_path;
    }
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    if (target_group_path != group_path) {
	std::string old_parent, old_leaf, new_parent, new_leaf;
	if (ged_obol_group_parent_leaf(group_path, old_parent, old_leaf) &&
	    ged_obol_group_parent_leaf(target_group_path, new_parent,
				       new_leaf) &&
	    old_parent == new_parent &&
	    scene->renameGroup(group_path.c_str(), new_leaf.c_str()) > 0) {
	    group_path = target_group_path;
	    group = scene->findGroup(group_path.c_str());
	}
	if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    return 0;
    }

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    int next_draw_mode = mode_valid ?
			 ged_obol_lod_draw_mode_from_ged(mode) :
			 scene_group->drawMode.getValue();
    int fallback_draw_mode = scene_group->fallbackDrawMode.getValue();
    if (fallback_draw_mode == BOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BOBOL_LOD_DRAW_WIRE;
    SbBool next_overlay = overlay_valid ?
			  (overlay ? TRUE : FALSE) : scene_group->overlayIntent.getValue();
    uint32_t revalidation_revision =
	scene_group->revalidationRevision.getValue();

    std::string next_intent_path;
    if (intent_path && intent_path[0]) {
	next_intent_path =
	    ged_obol_group_intent_path(target_group_path.c_str());
    } else {
	const char *retained_intent =
	    scene_group->drawIntentPath.getValue().getString();
	if (retained_intent && retained_intent[0])
	    next_intent_path = retained_intent;
	else
	    next_intent_path = ged_obol_group_intent_path(group_path.c_str());
    }

    int changed = scene->setGroupDrawIntent(group_path.c_str(),
					    next_intent_path.c_str(),
					    next_draw_mode,
					    fallback_draw_mode,
					    next_overlay,
					    revalidation_revision);
    return changed >= 0 ? 1 : 0;
}

int
ged_obol_database_source_exact_draw_mode_to_ged(
    struct ged *gedp,
    const BObolDatabaseSourceSummary &summary,
    SoBRLDatabaseSource *source)
{
    if (summary.representationMode >= 0)
	return summary.representationMode;

    const int source_ged_mode =
	ged_obol_database_draw_mode_to_ged(summary.drawMode);
    int exact_ged_mode = source_ged_mode;
    if (source_ged_mode == GED_DRAW_MODE_SHADED && source &&
	source->hasRealizedMeshGeometry())
	exact_ged_mode = GED_DRAW_MODE_SHADED_BOTS;

    const std::string owner_group_path =
	ged_obol_database_source_owner_group_path_from_summary(summary);
    if (owner_group_path.empty())
	return exact_ged_mode;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    SoGroup *group = scene ? scene->findGroup(owner_group_path.c_str()) :
		     NULL;
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return exact_ged_mode;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    const int group_ged_mode =
	ged_obol_lod_draw_mode_to_ged(scene_group->drawMode.getValue());
    if (group_ged_mode != GED_DRAW_MODE_WIRE ||
	source_ged_mode == GED_DRAW_MODE_WIRE)
	return group_ged_mode;
    return exact_ged_mode;
}

extern "C" int
ged_draw_obol_database_source_display_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_scene_display_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!ged_obol_database_source_controller_summary_for_source(scene,
	    source, summary))
	return 0;
    if (!summary.valid)
	return 0;

    int exact_draw_mode =
	ged_obol_database_source_exact_draw_mode_to_ged(gedp, summary,
	    source);

    out->valid = 1;
    out->is_database_source = 1;
    out->has_draw_intent = 1;
    out->intent_path = source->path.getValue().getString();
    out->intent_draw_mode = exact_draw_mode;
    out->visible = summary.visible ? 1 : 0;
    out->selected = summary.selected ? 1 : 0;
    out->highlighted = summary.highlighted ? 1 : 0;
    out->line_style = summary.lineStyle;
    out->line_width = summary.lineWidth;
    out->transparency = ged_obol_reported_transparency(summary.transparency);
    out->draw_mode = exact_draw_mode;
    out->material_valid = (summary.materialColorValid ||
			   summary.databaseMaterialColorValid ||
			   summary.colorOverride) ? 1 : 0;
    if (out->material_valid)
	ged_obol_rgb_from_color(ged_obol_summary_material_color(summary),
				out->material_color);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_draw_state_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_obol_draw_state_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    SoBRLDatabaseSource *source = NULL;
    BObolDatabaseSourceSummary source_summary;
    int have_source = 0;
    int best_score = -1;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *candidate =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	BObolDatabaseSourceSummary candidate_summary;
	if (!ged_obol_database_source_controller_summary_for_source(scene,
		candidate, candidate_summary) || !candidate_summary.valid)
	    continue;

	int score = 0;
	if (candidate_summary.drawMatrixValid)
	    score += 4;
	if (candidate_summary.lineStyle)
	    score += 2;
	if (candidate_summary.lineWidth)
	    score += 1;
	if (!have_source || score > best_score) {
	    source = candidate;
	    source_summary = candidate_summary;
	    have_source = 1;
	    best_score = score;
	}
    }
    if (!have_source || !source)
	return 0;

    out->valid = 1;
    out->draw_mode_valid = 1;
    out->draw_mode = ged_obol_database_source_exact_draw_mode_to_ged(gedp,
		     source_summary, source);
    out->line_style = source_summary.lineStyle;
    if (source_summary.drawMatrixValid) {
	out->draw_mat_valid = 1;
	ged_obol_mat_from_sbmatrix(source_summary.drawMatrix, out->draw_mat);
    }

    if (!out->draw_mat_valid) {
	const int count = source->getRealizedDisplaySummaryCount();
	for (int i = 0; i < count; i++) {
	    BObolSceneDisplaySummary display_summary;
	    if (!source->getRealizedDisplaySummary(i, display_summary) ||
		!display_summary.valid ||
		!display_summary.drawMatrixValid)
		continue;
	    if (display_summary.nodeKind !=
		BObolSceneTreeSummary::NODE_VLIST_SHAPE &&
		display_summary.nodeKind !=
		BObolSceneTreeSummary::NODE_MESH_SHAPE)
		continue;
	    out->draw_mat_valid = 1;
	    ged_obol_mat_from_sbmatrix(display_summary.drawMatrix,
				       out->draw_mat);
	    break;
	}
    }

    const auto matrix_delta = [](const mat_t mat) -> fastf_t {
	fastf_t delta = 0.0;
	for (int i = 0; i < 16; i++) {
	    const fastf_t expected =
		(i == 0 || i == 5 || i == 10 || i == 15) ? 1.0 : 0.0;
	    delta += fabs(mat[i] - expected);
	}
	return delta;
    };

    fastf_t best_matrix_delta =
	out->draw_mat_valid ? matrix_delta(out->draw_mat) : -1.0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *candidate =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	BObolDatabaseSourceSummary candidate_summary;
	if (!ged_obol_database_source_controller_summary_for_source(scene,
		candidate, candidate_summary) || !candidate_summary.valid)
	    continue;
	if (candidate_summary.lineStyle > out->line_style)
	    out->line_style = candidate_summary.lineStyle;
	if (candidate_summary.drawMatrixValid) {
	    mat_t candidate_mat;
	    ged_obol_mat_from_sbmatrix(candidate_summary.drawMatrix,
				       candidate_mat);
	    const fastf_t candidate_delta = matrix_delta(candidate_mat);
	    if (!out->draw_mat_valid ||
		candidate_delta > best_matrix_delta + VUNITIZE_TOL) {
		out->draw_mat_valid = 1;
		MAT_COPY(out->draw_mat, candidate_mat);
		best_matrix_delta = candidate_delta;
	    }
	}
    }

    return 1;
}
