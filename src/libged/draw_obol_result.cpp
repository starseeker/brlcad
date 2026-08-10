/*                D R A W _ O B O L _ R E S U L T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file draw_obol_result.cpp
 *
 * Single-responsibility GED command-result scene adapter.  A result scene owns
 * only its command transaction record; geometry and metadata are published
 * directly into the endpoint-selected libBObol feature store.
 */

#include "common.h"

#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "bg/line_layer.h"
#include "ged.h"
#include "ged/result.h"
#include "ged/view_feature.h"

#include "./draw_obol_bridge_private.hpp"

#include <Inventor/nodes/SoNode.h>

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#define GED_RESULT_SCENE_MAGIC 0x67646373u

struct ged_result_scene {
    uint32_t magic;
    struct ged_view_context *view_ctx;
    BObolViewController *controller;
    BObolFeatureOwner owner;
    BObolFeatureScope scope;
    ged_result_callback result_cb;
    void *result_cb_data;
    int changed;
    int failed;
};

static int
ged_obol_command_scene_valid(const struct ged_result_scene *scene)
{
    return scene && scene->magic == GED_RESULT_SCENE_MAGIC &&
	   scene->controller;
}

static void
ged_obol_command_scene_notify(
    struct ged_result_scene *scene,
    int status,
    const char *name,
    const char *command,
    const char *diagnostic,
    BObolFeatureHandle handle = BObolFeatureHandle())
{
    if (!ged_obol_command_scene_valid(scene) || !scene->result_cb)
	return;

    struct ged_result_event result =
	    GED_RESULT_INIT;
    result.status = status;
    result.feature_name = name;
    result.command = command;
    result.diagnostic = diagnostic;
    if (handle.isValid()) {
	result.feature_id = handle.id;
	result.feature_revision = handle.revision;
    }
    scene->result_cb(&result, scene->result_cb_data);
}

static int
ged_obol_command_scene_writable(
    struct ged_result_scene *scene,
    const char *name,
    const char *command)
{
    if (!ged_obol_command_scene_valid(scene))
	return 0;

    if (!scene->controller->features().commandOwnerGenerationCurrent(
	    scene->owner)) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name, command,
				      "stale command-scene generation");
	return 0;
    }

    return 1;
}

static BObolFeatureOwner
ged_obol_command_scene_owner(
    struct ged_view_context *view_ctx,
    const struct ged_result_desc *desc)
{
    BObolFeatureOwner owner;
    const char *owner_id = (desc && desc->owner_id && desc->owner_id[0]) ?
			   desc->owner_id : "ged-command";
    const char *owner_role =
	(desc && desc->owner_role && desc->owner_role[0]) ?
	desc->owner_role : "command-scene";
    const char *run_id = (desc && desc->run_id && desc->run_id[0]) ?
			 desc->run_id : NULL;

    std::string id(owner_id);
    if (run_id) {
	id += "#";
	id += run_id;
    }
    owner.ownerToken = NULL;
    owner.ownerId = id.c_str();
    owner.ownerRole = owner_role;
    owner.generation = desc ? desc->generation : 0;
    if (owner.ownerId.getLength() == 0)
	owner.ownerId = ged_obol_view_scope_name(view_ctx).c_str();
    return owner;
}

static BObolOverlayInfo
ged_obol_command_scene_overlay_info(
    const struct ged_result_scene *scene,
    const char *source_path)
{
    BObolOverlayInfo info =
	ged_obol_model_overlay_info(scene ? scene->view_ctx : NULL,
				    BObolOverlayClass::CommandResult,
				    BObolOverlayLifecycle::PerCommand,
				    BObolOverlayOrder::PostTransparent,
				    source_path);
    return info;
}

static int
ged_obol_command_scene_remove_feature(
    struct ged_result_scene *scene,
    const char *name,
    const char *command)
{
    if (!ged_obol_command_scene_writable(scene, name, command) || !name)
	return 0;

    const unsigned int scope_mask =
	scene->scope == BObolFeatureScope::Local ?
	BOBOL_FEATURE_SCOPE_LOCAL : BOBOL_FEATURE_SCOPE_SHARED;
    BObolFeatureHandle handle =
	scene->controller->features().findOwned(name, scope_mask,
	    &scene->owner);
    const int removed = scene->controller->features().removeOwned(name,
			scope_mask, &scene->owner) ? 1 : 0;
    if (removed)
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_REMOVED, name,
				      command ? command : "remove", NULL, handle);
    return removed;
}

extern "C" struct ged_result_scene *
ged_result_begin(
    struct ged_view_context *view_ctx,
    const struct ged_result_desc *desc)
{
    const int local = desc && desc->local;
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller)
	return NULL;

    struct ged_result_scene *scene = new ged_result_scene;
    scene->magic = GED_RESULT_SCENE_MAGIC;
    scene->view_ctx = view_ctx;
    scene->controller = controller;
    scene->owner = ged_obol_command_scene_owner(view_ctx, desc);
    scene->controller->features().markCommandOwnerGeneration(scene->owner);
    scene->scope = local ?
		   BObolFeatureScope::Local : BObolFeatureScope::Shared;
    scene->result_cb = desc ? desc->result_cb : NULL;
    scene->result_cb_data = desc ? desc->result_cb_data : NULL;
    scene->changed = 0;
    scene->failed = 0;
    ged_obol_command_scene_notify(scene, GED_RESULT_ACCEPTED,
				  NULL, "begin", NULL);
    return scene;
}

extern "C" size_t
ged_result_features_remove_prefix(
    struct ged_result_scene *scene,
    const char *prefix)
{
    if (!prefix) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, NULL, "removePrefix",
					  "missing feature prefix");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, prefix, "removePrefix"))
	return 0;

    size_t removed = 0;
    if (scene->scope == BObolFeatureScope::Local) {
	removed = scene->controller->features().removePrefix(prefix,
		  BOBOL_FEATURE_SCOPE_LOCAL, &scene->owner);
    } else {
	removed = scene->controller->features().removePrefix(prefix,
		  BOBOL_FEATURE_SCOPE_SHARED, &scene->owner);
    }
    scene->changed += static_cast<int>(removed);
    if (removed) {
	char diagnostic[64];
	snprintf(diagnostic, sizeof(diagnostic), "%zu feature(s) removed",
		 removed);
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_REMOVED, prefix, "removePrefix",
				      diagnostic);
    }
    return removed;
}

extern "C" int
ged_result_line_layer_builder_replace(
    struct ged_result_scene *scene,
    const char *name,
    const struct bg_line_layer_builder *builder,
    const struct ged_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, NULL,
					  "lineLayerBuilderReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name,
					 "lineLayerBuilderReplace"))
	return 0;

    if (!builder || !bg_line_layer_builder_point_count(builder)) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"lineLayerBuilderReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle =
	scene->controller->features().publishLineLayerBuilder(name,
	    scene->scope, builder, style ? &obol_style : NULL,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name,
				      "lineLayerBuilderReplace", "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_RESULT_UPDATED, name,
				  "lineLayerBuilderReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_result_line_layers_replace(
    struct ged_result_scene *scene,
    const char *name,
    const struct ged_annotation_line_layer *layers,
    size_t layer_count,
    const struct ged_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, NULL,
					  "lineLayersReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name,
					 "lineLayersReplace"))
	return 0;

    size_t live_layers = 0;
    for (size_t i = 0; layers && i < layer_count; i++) {
	if (layers[i].points && layers[i].point_count)
	    live_layers++;
    }
    if (!layers || !live_layers) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"lineLayersReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    std::vector<BObolLineLayer> obol_layers;
    obol_layers.reserve(live_layers);
    for (size_t i = 0; i < layer_count; i++) {
	if (!layers[i].points || !layers[i].point_count)
	    continue;

	BObolLineLayer layer;
	layer.name = layers[i].name ? layers[i].name : name;
	layer.style = ged_obol_feature_style_from_ged(&layers[i].style);
	layer.points.reserve(layers[i].point_count);
	layer.commands.reserve(layers[i].point_count);
	for (size_t j = 0; j < layers[i].point_count; j++) {
	    layer.points.push_back(SbVec3f(
				       static_cast<float>(layers[i].points[j][0]),
				       static_cast<float>(layers[i].points[j][1]),
				       static_cast<float>(layers[i].points[j][2])));
	    const int command = layers[i].commands ?
				layers[i].commands[j] : -1;
	    layer.commands.push_back(ged_obol_line_command_from_ged(command,
				     j));
	}
	obol_layers.push_back(layer);
    }

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle =
	scene->controller->features().publishLineLayers(name,
	    scene->scope, obol_layers, style ? &obol_style : NULL,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name,
				      "lineLayersReplace", "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_RESULT_UPDATED, name, "lineLayersReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_result_line_set_replace(
    struct ged_result_scene *scene,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, NULL,
					  "lineSetReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name, "lineSetReplace"))
	return 0;

    if (point_count && !points) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name, "lineSetReplace",
				      "missing line payload");
	return 0;
    }

    if (!points || !point_count) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"lineSetReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle =
	scene->controller->features().publishLineSet(name, scene->scope,
	    ged_obol_points_from_ged(points, point_count),
	    ged_obol_commands_from_ged(cmds, point_count),
	    style ? &obol_style : NULL, &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name, "lineSetReplace",
				      "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_RESULT_UPDATED, name, "lineSetReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_result_point_set_replace(
    struct ged_result_scene *scene,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, NULL,
					  "pointSetReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name, "pointSetReplace"))
	return 0;

    if (point_count && !points) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name, "pointSetReplace",
				      "missing point payload");
	return 0;
    }

    if (!points || !point_count) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"pointSetReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle =
	scene->controller->features().publishPointSet(name, scene->scope,
	    ged_obol_points_from_ged(points, point_count),
	    style ? &obol_style : NULL, &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name, "pointSetReplace",
				      "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_RESULT_UPDATED, name, "pointSetReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_result_indexed_face_set_replace(
    struct ged_result_scene *scene,
    const char *name,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count,
    const struct ged_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, NULL,
					  "indexedFaceSetReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name,
					 "indexedFaceSetReplace"))
	return 0;

    if ((point_count && !points) || (index_count && !indices)) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name,
				      "indexedFaceSetReplace", "missing indexed-face payload");
	return 0;
    }

    if (!points || !point_count || !indices || !index_count) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"indexedFaceSetReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle =
	scene->controller->features().publishIndexedFaceSet(name,
	    scene->scope,
	    ged_obol_points_from_ged(points, point_count),
	    ged_obol_vectors_from_ged(normals, normal_count),
	    ged_obol_indices_from_ged(indices, index_count),
	    style ? &obol_style : NULL, &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name,
				      "indexedFaceSetReplace", "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_RESULT_UPDATED, name, "indexedFaceSetReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_result_hud_label_replace(
    struct ged_result_scene *scene,
    const char *name,
    const struct ged_diagnostic_hud_label *label)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, NULL, "hudLabelReplace",
					  "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name, "hudLabelReplace"))
	return 0;

    if (!label || !label->text || !label->text[0]) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"hudLabelReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    BObolFeatureStyle obol_style;
    obol_style.hasVisible = TRUE;
    obol_style.visible = TRUE;
    obol_style.hasSelectable = TRUE;
    obol_style.selectable = TRUE;

    std::vector<BObolLabel> labels;
    labels.push_back(ged_obol_label_from_hud(*label));
    BObolFeatureHandle handle =
	scene->controller->features().publishHudLabels(name,
	    scene->scope, labels, &obol_style, &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name, "hudLabelReplace",
				      "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_RESULT_UPDATED, name, "hudLabelReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_result_custom_node_replace(
    struct ged_result_scene *scene,
    const char *name,
    ged_result_custom_node_cb node_cb,
    void *node_cb_data,
    const struct ged_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, NULL,
					  "customNodeReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name,
					 "customNodeReplace"))
	return 0;

    if (!node_cb) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"customNodeReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    struct ged_result_custom_node_request request =
	    GED_RESULT_SCENE_CUSTOM_NODE_REQUEST_INIT;
    request.feature_name = name;
    request.owner_id = scene->owner.ownerId.getString();
    request.owner_role = scene->owner.ownerRole.getString();
    request.generation = scene->owner.generation;
    request.local = scene->scope == BObolFeatureScope::Local ? 1 : 0;

    SoNode *node = static_cast<SoNode *>(node_cb(&request,
					 node_cb_data));
    if (!node) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name,
				      "customNodeReplace", "custom Coin node provider returned null");
	return 0;
    }

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle =
	scene->controller->features().publishCustomNode(name,
	    scene->scope, node, style ? &obol_style : NULL,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name,
				      "customNodeReplace", "custom Coin node publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_RESULT_UPDATED, name,
				  "customNodeReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_result_feature_metadata_replace(
    struct ged_result_scene *scene,
    const char *name,
    const struct ged_result_metadata *metadata,
    size_t metadata_count)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, NULL, "metadataReplace",
					  "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name, "metadataReplace"))
	return 0;

    BObolFeatureHandle handle =
	scene->controller->features().findOwned(name,
	    scene->scope == BObolFeatureScope::Local ?
	    BOBOL_FEATURE_SCOPE_LOCAL : BOBOL_FEATURE_SCOPE_SHARED,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name, "metadataReplace",
				      "owned feature not found");
	return 0;
    }

    std::vector<BObolFeatureMetadata> obol_metadata;
    obol_metadata.reserve(metadata_count);
    for (size_t i = 0; metadata && i < metadata_count; i++) {
	if (!metadata[i].key || !metadata[i].key[0])
	    continue;
	BObolFeatureMetadata item;
	item.key = metadata[i].key;
	item.value = metadata[i].value ? metadata[i].value : "";
	obol_metadata.push_back(item);
    }

    if (!scene->controller->features().replaceMetadata(handle,
	    obol_metadata)) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name, "metadataReplace",
				      "metadata replace failed", handle);
	return 0;
    }

    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_RESULT_UPDATED, name, "metadataReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_result_feature_primitive_metadata_replace(
    struct ged_result_scene *scene,
    const char *name,
    int primitive,
    const struct ged_result_metadata *metadata,
    size_t metadata_count)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, NULL,
					  "primitiveMetadataReplace", "missing feature name");
	}
	return 0;
    }
    if (primitive < 0) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_RESULT_FAILED, name,
					  "primitiveMetadataReplace",
					  "invalid primitive index");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name,
					 "primitiveMetadataReplace"))
	return 0;

    BObolFeatureHandle handle =
	scene->controller->features().findOwned(name,
	    scene->scope == BObolFeatureScope::Local ?
	    BOBOL_FEATURE_SCOPE_LOCAL : BOBOL_FEATURE_SCOPE_SHARED,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name,
				      "primitiveMetadataReplace", "owned feature not found");
	return 0;
    }

    std::vector<BObolFeatureMetadata> obol_metadata;
    obol_metadata.reserve(metadata_count);
    for (size_t i = 0; metadata && i < metadata_count; i++) {
	if (!metadata[i].key || !metadata[i].key[0])
	    continue;
	BObolFeatureMetadata item;
	item.key = metadata[i].key;
	item.value = metadata[i].value ? metadata[i].value : "";
	obol_metadata.push_back(item);
    }

    if (!scene->controller->features().replacePrimitiveMetadata(handle,
	    static_cast<int32_t>(primitive), obol_metadata)) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_RESULT_FAILED, name,
				      "primitiveMetadataReplace",
				      "primitive metadata replace failed", handle);
	return 0;
    }

    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_RESULT_UPDATED, name,
				  "primitiveMetadataReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_result_commit(struct ged_result_scene *scene)
{
    if (!ged_obol_command_scene_valid(scene))
	return 0;

    const int ret = (scene->failed ||
		     !scene->controller->features().commandOwnerGenerationCurrent(
			 scene->owner)) ? 0 : 1;
    ged_obol_command_scene_notify(scene,
				  ret ? GED_RESULT_ACCEPTED :
				  GED_RESULT_FAILED,
				  NULL, "commit",
				  ret ? NULL : "command-scene commit rejected");
    scene->magic = 0;
    delete scene;
    return ret;
}

extern "C" void
ged_result_abort(struct ged_result_scene *scene)
{
    if (!ged_obol_command_scene_valid(scene))
	return;

    scene->magic = 0;
    delete scene;
}

