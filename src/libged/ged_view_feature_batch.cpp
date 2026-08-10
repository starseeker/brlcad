/*              G E D _ V I E W _ F E A T U R E _ B A T C H . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_view_feature_batch.cpp
 *
 * Single-responsibility GED view-feature batch adapter.  A batch owns deep
 * copies of its staged geometry and metadata until commit, then publishes the
 * complete operation sequence into the endpoint-selected feature store.
 */

#include "common.h"

#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "bg/line_layer.h"
#include "ged.h"
#include "ged/view_feature_batch.h"
#include "ged/view_feature.h"

#include "./draw_obol_bridge_private.hpp"

#include <cstdio>
#include <cstring>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#define GED_VIEW_FEATURE_BATCH_MAGIC 0x67646373u

struct ged_view_feature_batch {
    uint32_t magic;
    struct ged_view_context *view_ctx;
    BObolViewController *controller;
    BObolFeatureOwner owner;
    BObolFeatureScope scope;
    ged_view_feature_batch_callback event_cb;
    void *event_cb_data;
    std::vector<std::function<int(struct ged_view_feature_batch *)>> operations;
    int committing;
    int changed;
    int failed;
};

struct ged_view_feature_staged_style {
    struct ged_view_feature_style value = GED_VIEW_FEATURE_STYLE_INIT;
    bool valid = false;
};

struct ged_view_feature_staged_line_layer {
    std::string name;
    std::vector<fastf_t> points;
    std::vector<int> commands;
    struct ged_view_feature_style style = GED_VIEW_FEATURE_STYLE_INIT;
};

static ged_view_feature_staged_style
ged_view_feature_stage_style(const struct ged_view_feature_style *style)
{
    ged_view_feature_staged_style staged;
    if (style) {
	staged.value = *style;
	staged.valid = true;
    }
    return staged;
}

static std::vector<fastf_t>
ged_view_feature_stage_vectors(const fastf_t *values, size_t vector_count)
{
    if (!values || !vector_count)
	return {};
    return std::vector<fastf_t>(values, values + vector_count * 3);
}

static int
ged_obol_feature_batch_valid(const struct ged_view_feature_batch *scene)
{
    return scene && scene->magic == GED_VIEW_FEATURE_BATCH_MAGIC &&
	   scene->controller;
}

static void
ged_obol_feature_batch_notify(
    struct ged_view_feature_batch *scene,
    enum ged_view_feature_batch_status status,
    const char *name,
    const char *command,
    const char *diagnostic,
    BObolFeatureHandle handle = BObolFeatureHandle())
{
    if (!ged_obol_feature_batch_valid(scene) || !scene->event_cb)
	return;

    struct ged_view_feature_batch_event result =
	    GED_VIEW_FEATURE_BATCH_EVENT_INIT;
    result.status = status;
    result.feature_name = name;
    result.command = command;
    result.diagnostic = diagnostic;
    if (handle.isValid()) {
	result.feature_id = handle.id;
	result.feature_revision = handle.revision;
    }
    scene->event_cb(&result, scene->event_cb_data);
}

static int
ged_obol_feature_batch_writable(
    struct ged_view_feature_batch *scene,
    const char *name,
    const char *command)
{
    if (!ged_obol_feature_batch_valid(scene))
	return 0;

    if (!scene->controller->features().commandOwnerGenerationCurrent(
	    scene->owner)) {
	scene->failed = 1;
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name, command,
				      "stale feature-batch generation");
	return 0;
    }

    return 1;
}

static BObolFeatureOwner
ged_obol_feature_batch_owner(
    struct ged_view_context *view_ctx,
    const struct ged_view_feature_batch_desc *desc)
{
    BObolFeatureOwner owner;
    const char *owner_id = (desc && desc->owner_id && desc->owner_id[0]) ?
			   desc->owner_id : "ged-command";
    const char *owner_role =
	(desc && desc->owner_role && desc->owner_role[0]) ?
	desc->owner_role : "feature-batch";
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
ged_obol_feature_batch_overlay_info(
    const struct ged_view_feature_batch *scene,
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
ged_obol_feature_batch_remove_feature(
    struct ged_view_feature_batch *scene,
    const char *name,
    const char *command)
{
    if (!ged_obol_feature_batch_writable(scene, name, command) || !name)
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
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_REMOVED, name,
				      command ? command : "remove", NULL, handle);
    return removed;
}

extern "C" struct ged_view_feature_batch *
ged_view_feature_batch_begin(
    struct ged_view_context *view_ctx,
    const struct ged_view_feature_batch_desc *desc)
{
    const int local = desc && desc->local;
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller)
	return NULL;

    struct ged_view_feature_batch *scene = new ged_view_feature_batch;
    scene->magic = GED_VIEW_FEATURE_BATCH_MAGIC;
    scene->view_ctx = view_ctx;
    scene->controller = controller;
    scene->owner = ged_obol_feature_batch_owner(view_ctx, desc);
    scene->controller->features().markCommandOwnerGeneration(scene->owner);
    scene->scope = local ?
		   BObolFeatureScope::Local : BObolFeatureScope::Shared;
    scene->event_cb = desc ? desc->event_cb : NULL;
    scene->event_cb_data = desc ? desc->event_cb_data : NULL;
    scene->committing = 0;
    scene->changed = 0;
    scene->failed = 0;
    ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_ACCEPTED,
				  NULL, "begin", NULL);
    return scene;
}

extern "C" size_t
ged_view_feature_batch_remove_prefix(
    struct ged_view_feature_batch *scene,
    const char *prefix)
{
    if (!prefix) {
	if (ged_obol_feature_batch_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_feature_batch_notify(scene,
					  GED_VIEW_FEATURE_BATCH_FAILED, NULL, "removePrefix",
					  "missing feature prefix");
	}
	return 0;
    }
    if (!ged_obol_feature_batch_valid(scene))
	return 0;

    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, prefix, "removePrefix"))
	    return 0;
	const std::string staged_prefix(prefix);
	const unsigned int scope_mask =
	    scene->scope == BObolFeatureScope::Local ?
	    BOBOL_FEATURE_SCOPE_LOCAL : BOBOL_FEATURE_SCOPE_SHARED;
	const size_t planned = scene->controller->features().countPrefix(
	    prefix, scope_mask, &scene->owner);
	scene->operations.emplace_back(
	    [staged_prefix](struct ged_view_feature_batch *batch) {
		(void)ged_view_feature_batch_remove_prefix(batch,
		    staged_prefix.c_str());
		return batch->failed ? 0 : 1;
	    });
	return planned;
    }

    if (!ged_obol_feature_batch_writable(scene, prefix, "removePrefix"))
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
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_REMOVED, prefix, "removePrefix",
				      diagnostic);
    }
    return removed;
}

extern "C" int
ged_view_feature_batch_line_layer_builder_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const struct bg_line_layer_builder *builder,
    const struct ged_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_feature_batch_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_feature_batch_notify(scene,
					  GED_VIEW_FEATURE_BATCH_FAILED, NULL,
					  "lineLayerBuilderReplace", "missing feature name");
	}
	return 0;
    }
    if (!ged_obol_feature_batch_valid(scene))
	return 0;

    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name,
		"lineLayerBuilderReplace"))
	    return 0;
	const std::string staged_name(name);
	const ged_view_feature_staged_style staged_style =
	    ged_view_feature_stage_style(style);
	std::shared_ptr<struct bg_line_layer_builder> staged_builder;
	if (builder && bg_line_layer_builder_point_count(builder)) {
	    staged_builder.reset(bg_line_layer_builder_create(),
		[](struct bg_line_layer_builder *value) {
		    bg_line_layer_builder_free(value);
		});
	    const size_t layer_count =
		bg_line_layer_builder_layer_count(builder);
	    for (size_t i = 0; i < layer_count; i++) {
		const struct bg_line_layer *layer =
		    bg_line_layer_builder_layer_at(builder, i);
		if (!layer)
		    continue;
		unsigned char r = 255, g = 255, b = 255;
		(void)bg_line_layer_color(layer, &r, &g, &b);
		const size_t point_count = bg_line_layer_point_count(layer);
		for (size_t j = 0; j < point_count; j++) {
		    point_t point;
		    int command = BG_GEOMETRY_LINE_DRAW;
		    if (!bg_line_layer_point_at(layer, j, point) ||
			!bg_line_layer_command_at(layer, j, &command) ||
			!bg_line_layer_builder_add(staged_builder.get(), r, g,
			    b, point, command)) {
			scene->failed = 1;
			return 0;
		    }
		}
	    }
	}
	scene->operations.emplace_back(
	    [staged_name, staged_builder, staged_style](
		struct ged_view_feature_batch *batch) {
		return ged_view_feature_batch_line_layer_builder_replace(
		    batch, staged_name.c_str(), staged_builder.get(),
		    staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }

    if (!ged_obol_feature_batch_writable(scene, name,
					 "lineLayerBuilderReplace"))
	return 0;

    if (!builder || !bg_line_layer_builder_point_count(builder)) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
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
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name,
				      "lineLayerBuilderReplace", "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_feature_batch_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene,
				  GED_VIEW_FEATURE_BATCH_UPDATED, name,
				  "lineLayerBuilderReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_line_layers_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const struct ged_annotation_line_layer *layers,
    size_t layer_count,
    const struct ged_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_feature_batch_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_feature_batch_notify(scene,
					  GED_VIEW_FEATURE_BATCH_FAILED, NULL,
					  "lineLayersReplace", "missing feature name");
	}
	return 0;
    }
    if (!ged_obol_feature_batch_valid(scene))
	return 0;

    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name,
		"lineLayersReplace"))
	    return 0;
	const std::string staged_name(name);
	const ged_view_feature_staged_style staged_style =
	    ged_view_feature_stage_style(style);
	std::vector<ged_view_feature_staged_line_layer> staged_layers;
	staged_layers.reserve(layer_count);
	for (size_t i = 0; layers && i < layer_count; i++) {
	    ged_view_feature_staged_line_layer staged;
	    staged.name = layers[i].name ? layers[i].name : "";
	    staged.style = layers[i].style;
	    staged.points = ged_view_feature_stage_vectors(
		reinterpret_cast<const fastf_t *>(layers[i].points),
		layers[i].point_count);
	    if (layers[i].commands && layers[i].point_count)
		staged.commands.assign(layers[i].commands,
		    layers[i].commands + layers[i].point_count);
	    staged_layers.push_back(std::move(staged));
	}
	scene->operations.emplace_back(
	    [staged_name, staged_style, staged_layers](
		struct ged_view_feature_batch *batch) {
		std::vector<struct ged_annotation_line_layer> layers_copy;
		layers_copy.reserve(staged_layers.size());
		for (const ged_view_feature_staged_line_layer &staged :
		    staged_layers) {
		    struct ged_annotation_line_layer layer =
			GED_ANNOTATION_LINE_LAYER_INIT;
		    layer.name = staged.name.empty() ? NULL :
			staged.name.c_str();
		    layer.points = staged.points.empty() ? NULL :
			reinterpret_cast<const point_t *>(staged.points.data());
		    layer.commands = staged.commands.empty() ? NULL :
			staged.commands.data();
		    layer.point_count = staged.points.size() / 3;
		    layer.style = staged.style;
		    layers_copy.push_back(layer);
		}
		return ged_view_feature_batch_line_layers_replace(batch,
		    staged_name.c_str(), layers_copy.empty() ? NULL :
		    layers_copy.data(), layers_copy.size(),
		    staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }

    if (!ged_obol_feature_batch_writable(scene, name,
					 "lineLayersReplace"))
	return 0;

    size_t live_layers = 0;
    for (size_t i = 0; layers && i < layer_count; i++) {
	if (layers[i].points && layers[i].point_count)
	    live_layers++;
    }
    if (!layers || !live_layers) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
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
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name,
				      "lineLayersReplace", "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_feature_batch_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene,
				  GED_VIEW_FEATURE_BATCH_UPDATED, name, "lineLayersReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_line_set_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_feature_batch_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_feature_batch_notify(scene,
					  GED_VIEW_FEATURE_BATCH_FAILED, NULL,
					  "lineSetReplace", "missing feature name");
	}
	return 0;
    }
    if (!ged_obol_feature_batch_valid(scene))
	return 0;

    if (point_count && !points) {
	if (ged_obol_feature_batch_valid(scene))
	    scene->failed = 1;
	return 0;
    }
    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name, "lineSetReplace"))
	    return 0;
	const std::string staged_name(name);
	const ged_view_feature_staged_style staged_style =
	    ged_view_feature_stage_style(style);
	const std::vector<fastf_t> staged_points =
	    ged_view_feature_stage_vectors(
		reinterpret_cast<const fastf_t *>(points), point_count);
	std::vector<int> staged_commands;
	if (cmds && point_count)
	    staged_commands.assign(cmds, cmds + point_count);
	scene->operations.emplace_back(
	    [staged_name, staged_style, staged_points, staged_commands](
		struct ged_view_feature_batch *batch) {
		return ged_view_feature_batch_line_set_replace(batch,
		    staged_name.c_str(), staged_points.empty() ? NULL :
		    reinterpret_cast<const point_t *>(staged_points.data()),
		    staged_commands.empty() ? NULL : staged_commands.data(),
		    staged_points.size() / 3,
		    staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }

    if (!ged_obol_feature_batch_writable(scene, name, "lineSetReplace"))
	return 0;

    if (point_count && !points) {
	scene->failed = 1;
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name, "lineSetReplace",
				      "missing line payload");
	return 0;
    }

    if (!points || !point_count) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
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
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name, "lineSetReplace",
				      "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_feature_batch_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene,
				  GED_VIEW_FEATURE_BATCH_UPDATED, name, "lineSetReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_point_set_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_feature_batch_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_feature_batch_notify(scene,
					  GED_VIEW_FEATURE_BATCH_FAILED, NULL,
					  "pointSetReplace", "missing feature name");
	}
	return 0;
    }
    if (!ged_obol_feature_batch_valid(scene))
	return 0;

    if (point_count && !points) {
	if (ged_obol_feature_batch_valid(scene))
	    scene->failed = 1;
	return 0;
    }
    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name, "pointSetReplace"))
	    return 0;
	const std::string staged_name(name);
	const ged_view_feature_staged_style staged_style =
	    ged_view_feature_stage_style(style);
	const std::vector<fastf_t> staged_points =
	    ged_view_feature_stage_vectors(
		reinterpret_cast<const fastf_t *>(points), point_count);
	scene->operations.emplace_back(
	    [staged_name, staged_style, staged_points](
		struct ged_view_feature_batch *batch) {
		return ged_view_feature_batch_point_set_replace(batch,
		    staged_name.c_str(), staged_points.empty() ? NULL :
		    reinterpret_cast<const point_t *>(staged_points.data()),
		    staged_points.size() / 3,
		    staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }

    if (!ged_obol_feature_batch_writable(scene, name, "pointSetReplace"))
	return 0;

    if (point_count && !points) {
	scene->failed = 1;
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name, "pointSetReplace",
				      "missing point payload");
	return 0;
    }

    if (!points || !point_count) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
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
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name, "pointSetReplace",
				      "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_feature_batch_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene,
				  GED_VIEW_FEATURE_BATCH_UPDATED, name, "pointSetReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_indexed_face_set_replace(
    struct ged_view_feature_batch *scene,
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
	if (ged_obol_feature_batch_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_feature_batch_notify(scene,
					  GED_VIEW_FEATURE_BATCH_FAILED, NULL,
					  "indexedFaceSetReplace", "missing feature name");
	}
	return 0;
    }
    if (!ged_obol_feature_batch_valid(scene))
	return 0;

    if ((point_count && !points) || (normal_count && !normals) ||
	(index_count && !indices)) {
	if (ged_obol_feature_batch_valid(scene))
	    scene->failed = 1;
	return 0;
    }
    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name,
		"indexedFaceSetReplace"))
	    return 0;
	const std::string staged_name(name);
	const ged_view_feature_staged_style staged_style =
	    ged_view_feature_stage_style(style);
	const std::vector<fastf_t> staged_points =
	    ged_view_feature_stage_vectors(
		reinterpret_cast<const fastf_t *>(points), point_count);
	const std::vector<fastf_t> staged_normals =
	    ged_view_feature_stage_vectors(
		reinterpret_cast<const fastf_t *>(normals), normal_count);
	std::vector<int> staged_indices;
	if (indices && index_count)
	    staged_indices.assign(indices, indices + index_count);
	scene->operations.emplace_back(
	    [staged_name, staged_style, staged_points, staged_normals,
		staged_indices](struct ged_view_feature_batch *batch) {
		return ged_view_feature_batch_indexed_face_set_replace(batch,
		    staged_name.c_str(), staged_points.empty() ? NULL :
		    reinterpret_cast<const point_t *>(staged_points.data()),
		    staged_points.size() / 3,
		    staged_normals.empty() ? NULL :
		    reinterpret_cast<const vect_t *>(staged_normals.data()),
		    staged_normals.size() / 3,
		    staged_indices.empty() ? NULL : staged_indices.data(),
		    staged_indices.size(),
		    staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }

    if (!ged_obol_feature_batch_writable(scene, name,
					 "indexedFaceSetReplace"))
	return 0;

    if ((point_count && !points) || (index_count && !indices)) {
	scene->failed = 1;
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name,
				      "indexedFaceSetReplace", "missing indexed-face payload");
	return 0;
    }

    if (!points || !point_count || !indices || !index_count) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
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
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name,
				      "indexedFaceSetReplace", "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_feature_batch_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene,
				  GED_VIEW_FEATURE_BATCH_UPDATED, name, "indexedFaceSetReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_hud_label_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const struct ged_diagnostic_hud_label *label)
{
    if (!name) {
	if (ged_obol_feature_batch_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_feature_batch_notify(scene,
					  GED_VIEW_FEATURE_BATCH_FAILED, NULL, "hudLabelReplace",
					  "missing feature name");
	}
	return 0;
    }
    if (!ged_obol_feature_batch_valid(scene))
	return 0;

    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name, "hudLabelReplace"))
	    return 0;
	const std::string staged_name(name);
	const bool has_label = label && label->text && label->text[0];
	const std::string staged_label_id =
	    label && label->label_id ? label->label_id : "";
	const std::string staged_text = has_label ? label->text : "";
	struct ged_diagnostic_hud_label staged_label =
	    GED_DIAGNOSTIC_HUD_LABEL_INIT;
	if (label)
	    staged_label = *label;
	scene->operations.emplace_back(
	    [staged_name, has_label, staged_label_id, staged_text,
		staged_label](struct ged_view_feature_batch *batch) {
		if (!has_label)
		    return ged_view_feature_batch_hud_label_replace(batch,
			staged_name.c_str(), NULL);
		struct ged_diagnostic_hud_label label_copy = staged_label;
		label_copy.label_id = staged_label_id.empty() ? NULL :
		    staged_label_id.c_str();
		label_copy.text = staged_text.c_str();
		return ged_view_feature_batch_hud_label_replace(batch,
		    staged_name.c_str(), &label_copy);
	    });
	return 1;
    }

    if (!ged_obol_feature_batch_writable(scene, name, "hudLabelReplace"))
	return 0;

    if (!label || !label->text || !label->text[0]) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
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
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name, "hudLabelReplace",
				      "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_feature_batch_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene,
				  GED_VIEW_FEATURE_BATCH_UPDATED, name, "hudLabelReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_metadata_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const struct ged_view_feature_metadata *metadata,
    size_t metadata_count)
{
    if (!name) {
	if (ged_obol_feature_batch_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_feature_batch_notify(scene,
					  GED_VIEW_FEATURE_BATCH_FAILED, NULL, "metadataReplace",
					  "missing feature name");
	}
	return 0;
    }
    if (!ged_obol_feature_batch_valid(scene))
	return 0;

    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name, "metadataReplace"))
	    return 0;
	const std::string staged_name(name);
	std::vector<std::pair<std::string, std::string>> staged_metadata;
	staged_metadata.reserve(metadata_count);
	for (size_t i = 0; metadata && i < metadata_count; i++) {
	    if (metadata[i].key && metadata[i].key[0])
		staged_metadata.emplace_back(metadata[i].key,
		    metadata[i].value ? metadata[i].value : "");
	}
	scene->operations.emplace_back(
	    [staged_name, staged_metadata](
		struct ged_view_feature_batch *batch) {
		std::vector<struct ged_view_feature_metadata> metadata_copy;
		metadata_copy.reserve(staged_metadata.size());
		for (const auto &item : staged_metadata) {
		    struct ged_view_feature_metadata metadata_item =
			GED_VIEW_FEATURE_METADATA_INIT;
		    metadata_item.key = item.first.c_str();
		    metadata_item.value = item.second.c_str();
		    metadata_copy.push_back(metadata_item);
		}
		return ged_view_feature_batch_metadata_replace(batch,
		    staged_name.c_str(), metadata_copy.empty() ? NULL :
		    metadata_copy.data(), metadata_copy.size());
	    });
	return 1;
    }

    if (!ged_obol_feature_batch_writable(scene, name, "metadataReplace"))
	return 0;

    BObolFeatureHandle handle =
	scene->controller->features().findOwned(name,
	    scene->scope == BObolFeatureScope::Local ?
	    BOBOL_FEATURE_SCOPE_LOCAL : BOBOL_FEATURE_SCOPE_SHARED,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name, "metadataReplace",
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
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name, "metadataReplace",
				      "metadata replace failed", handle);
	return 0;
    }

    scene->changed++;
    ged_obol_feature_batch_notify(scene,
				  GED_VIEW_FEATURE_BATCH_UPDATED, name, "metadataReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_primitive_metadata_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    int primitive,
    const struct ged_view_feature_metadata *metadata,
    size_t metadata_count)
{
    if (!name) {
	if (ged_obol_feature_batch_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_feature_batch_notify(scene,
					  GED_VIEW_FEATURE_BATCH_FAILED, NULL,
					  "primitiveMetadataReplace", "missing feature name");
	}
	return 0;
    }
    if (!ged_obol_feature_batch_valid(scene))
	return 0;
    if (primitive < 0) {
	if (ged_obol_feature_batch_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_feature_batch_notify(scene,
					  GED_VIEW_FEATURE_BATCH_FAILED, name,
					  "primitiveMetadataReplace",
					  "invalid primitive index");
	}
	return 0;
    }

    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name,
		"primitiveMetadataReplace"))
	    return 0;
	const std::string staged_name(name);
	std::vector<std::pair<std::string, std::string>> staged_metadata;
	staged_metadata.reserve(metadata_count);
	for (size_t i = 0; metadata && i < metadata_count; i++) {
	    if (metadata[i].key && metadata[i].key[0])
		staged_metadata.emplace_back(metadata[i].key,
		    metadata[i].value ? metadata[i].value : "");
	}
	scene->operations.emplace_back(
	    [staged_name, primitive, staged_metadata](
		struct ged_view_feature_batch *batch) {
		std::vector<struct ged_view_feature_metadata> metadata_copy;
		metadata_copy.reserve(staged_metadata.size());
		for (const auto &item : staged_metadata) {
		    struct ged_view_feature_metadata metadata_item =
			GED_VIEW_FEATURE_METADATA_INIT;
		    metadata_item.key = item.first.c_str();
		    metadata_item.value = item.second.c_str();
		    metadata_copy.push_back(metadata_item);
		}
		return ged_view_feature_batch_primitive_metadata_replace(
		    batch, staged_name.c_str(), primitive,
		    metadata_copy.empty() ? NULL : metadata_copy.data(),
		    metadata_copy.size());
	    });
	return 1;
    }

    if (!ged_obol_feature_batch_writable(scene, name,
					 "primitiveMetadataReplace"))
	return 0;

    BObolFeatureHandle handle =
	scene->controller->features().findOwned(name,
	    scene->scope == BObolFeatureScope::Local ?
	    BOBOL_FEATURE_SCOPE_LOCAL : BOBOL_FEATURE_SCOPE_SHARED,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name,
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
	ged_obol_feature_batch_notify(scene,
				      GED_VIEW_FEATURE_BATCH_FAILED, name,
				      "primitiveMetadataReplace",
				      "primitive metadata replace failed", handle);
	return 0;
    }

    scene->changed++;
    ged_obol_feature_batch_notify(scene,
				  GED_VIEW_FEATURE_BATCH_UPDATED, name,
				  "primitiveMetadataReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_commit(struct ged_view_feature_batch *scene)
{
    if (!ged_obol_feature_batch_valid(scene))
	return 0;

    int ret = (scene->failed ||
	!scene->controller->features().commandOwnerGenerationCurrent(
	    scene->owner)) ? 0 : 1;
    if (ret) {
	scene->committing = 1;
	for (const auto &operation : scene->operations) {
	    if (!operation || !operation(scene)) {
		scene->failed = 1;
		ret = 0;
		break;
	    }
	}
	scene->committing = 0;
    }
    ged_obol_feature_batch_notify(scene,
				  ret ? GED_VIEW_FEATURE_BATCH_ACCEPTED :
				  GED_VIEW_FEATURE_BATCH_FAILED,
				  NULL, "commit",
				  ret ? NULL : "feature-batch commit rejected");
    scene->magic = 0;
    delete scene;
    return ret;
}

extern "C" void
ged_view_feature_batch_abort(struct ged_view_feature_batch *scene)
{
    if (!ged_obol_feature_batch_valid(scene))
	return;

    scene->magic = 0;
    delete scene;
}
