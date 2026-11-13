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
#include "bv.h"
#include "ged.h"
#include "ged/view_feature_batch.h"
#include "ged/view_feature.h"

#include "./draw_obol_bridge_private.hpp"
#include "./ged_draw_private.h"

#include <algorithm>
#include <array>
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
    BObolOverlayClass overlay_class;
    BObolOverlayLifecycle lifecycle;
    BObolOverlayOrder overlay_order;
    ged_view_feature_batch_callback event_cb;
    void *event_cb_data;
    std::vector<std::function<int(struct ged_view_feature_batch *)>> operations;
    int committing;
    int changed;
    int failed;
};


void
ged_view_feature_batch_desc_init(struct ged_view_feature_batch_desc *desc)
{
    if (!desc)
	return;
    *desc = {};
    desc->overlay_class = GED_VIEW_FEATURE_OVERLAY_CLASS_COMMAND_RESULT;
    desc->lifecycle = GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND;
    desc->overlay_order = GED_VIEW_FEATURE_OVERLAY_ORDER_POST_TRANSPARENT;
}


void
ged_view_feature_batch_event_init(struct ged_view_feature_batch_event *event)
{
    if (!event)
	return;
    *event = {};
    event->status = GED_VIEW_FEATURE_BATCH_NONE;
}


void
ged_view_feature_line_layer_init(struct ged_view_feature_line_layer *layer)
{
    if (!layer)
	return;
    *layer = {};
    ged_view_feature_style_init(&layer->style);
}


void
ged_view_feature_metadata_init(struct ged_view_feature_metadata *metadata)
{
    if (!metadata)
	return;
    *metadata = {};
}

struct ged_view_feature_staged_style {
    struct ged_view_feature_style value = ged_view_feature_style_default();
    bool valid = false;
};

struct ged_view_feature_staged_line_layer {
    std::string name;
    std::vector<fastf_t> points;
    std::vector<int> commands;
    struct ged_view_feature_style style = ged_view_feature_style_default();
};

struct ged_view_feature_staged_label {
    std::string text;
    point_t point = VINIT_ZERO;
    int color_valid = 0;
    unsigned char color[3] = {0, 0, 0};
    int line_flag = 0;
    point_t target = VINIT_ZERO;
    int anchor = 0;
    int arrow = 0;
    fastf_t font_size = 0.0;
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

static BObolLabel
ged_view_feature_batch_label_to_obol(const struct ged_view_feature_label &src)
{
    BObolLabel label;
    label.text = src.text ? src.text : "";
    label.point = SbVec3f(static_cast<float>(src.point[X]),
	static_cast<float>(src.point[Y]),
	static_cast<float>(src.point[Z]));
    if (src.color_valid) {
	label.hasColor = TRUE;
	label.color = SbColor(src.color[0] / 255.0f,
	    src.color[1] / 255.0f, src.color[2] / 255.0f);
    }
    if (src.font_size > 0.0)
	label.fontSize = static_cast<float>(src.font_size);
    label.hasLeader = src.line_flag ? TRUE : FALSE;
    label.target = SbVec3f(static_cast<float>(src.target[X]),
	static_cast<float>(src.target[Y]),
	static_cast<float>(src.target[Z]));
    label.anchor = src.anchor;
    label.arrow = src.arrow ? TRUE : FALSE;
    return label;
}

static SbVec3f
ged_view_feature_batch_hud_point(struct ged_view_context *view_ctx,
	const point_t point)
{
    const struct bv *view = bv_context_view_const(
	(const struct bv_context *)view_ctx);
    const int width = view ? bv_width_get(view) : 0;
    const int height = view ? bv_height_get(view) : 0;
    return SbVec3f(static_cast<float>((point[X] + 1.0) * 0.5 * width),
	static_cast<float>((point[Y] + 1.0) * 0.5 * height), 0.0f);
}

static BObolOverlayInfo
ged_view_feature_batch_hud_overlay_info(
	const struct ged_view_feature_batch *scene, const char *name)
{
    return ged_obol_model_overlay_info(scene ? scene->view_ctx : NULL,
	BObolOverlayClass::Faceplate, BObolOverlayLifecycle::PerView,
	BObolOverlayOrder::PostTransparent, name);
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
	    ged_view_feature_batch_event_default();
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
				    scene ? scene->overlay_class :
				    BObolOverlayClass::CommandResult,
				    scene ? scene->lifecycle :
				    BObolOverlayLifecycle::PerCommand,
				    scene ? scene->overlay_order :
				    BObolOverlayOrder::PostTransparent,
				    source_path);
    return info;
}

static BObolOverlayClass
ged_obol_feature_batch_overlay_class(int value)
{
    switch (value) {
	case GED_VIEW_FEATURE_OVERLAY_CLASS_NONE: return BObolOverlayClass::None;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_FACEPLATE: return BObolOverlayClass::Faceplate;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_EDIT_HANDLE: return BObolOverlayClass::EditHandle;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_MEASURE: return BObolOverlayClass::Measure;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_SELECTION_RUBBER_BAND: return BObolOverlayClass::SelectionRubberBand;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_SNAP_GUIDE: return BObolOverlayClass::SnapGuide;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_DIAGNOSTIC: return BObolOverlayClass::Diagnostic;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_TCL_OVERLAY: return BObolOverlayClass::TclOverlay;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_POLYGON_EDIT: return BObolOverlayClass::PolygonEdit;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_SKETCH_EDIT: return BObolOverlayClass::SketchEdit;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION: return BObolOverlayClass::UserAnnotation;
	case GED_VIEW_FEATURE_OVERLAY_CLASS_COMMAND_RESULT:
	case GED_VIEW_FEATURE_OVERLAY_CLASS_UNKNOWN:
	default: return BObolOverlayClass::CommandResult;
    }
}

static BObolOverlayLifecycle
ged_obol_feature_batch_lifecycle(int value)
{
    switch (value) {
	case GED_VIEW_FEATURE_LIFECYCLE_NONE: return BObolOverlayLifecycle::None;
	case GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT: return BObolOverlayLifecycle::Persistent;
	case GED_VIEW_FEATURE_LIFECYCLE_PER_FRAME: return BObolOverlayLifecycle::PerFrame;
	case GED_VIEW_FEATURE_LIFECYCLE_PER_TOOL: return BObolOverlayLifecycle::PerTool;
	case GED_VIEW_FEATURE_LIFECYCLE_PER_VIEW: return BObolOverlayLifecycle::PerView;
	case GED_VIEW_FEATURE_LIFECYCLE_SHARED_VIEW_SET: return BObolOverlayLifecycle::SharedViewSet;
	case GED_VIEW_FEATURE_LIFECYCLE_AUTO_REMOVE_ON_SOURCE: return BObolOverlayLifecycle::AutoRemoveOnSource;
	case GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND:
	case GED_VIEW_FEATURE_LIFECYCLE_UNKNOWN:
	default: return BObolOverlayLifecycle::PerCommand;
    }
}

static BObolOverlayOrder
ged_obol_feature_batch_overlay_order(int value)
{
    switch (value) {
	case GED_VIEW_FEATURE_OVERLAY_ORDER_MODEL: return BObolOverlayOrder::Model;
	case GED_VIEW_FEATURE_OVERLAY_ORDER_SCREEN: return BObolOverlayOrder::Screen;
	case GED_VIEW_FEATURE_OVERLAY_ORDER_XRAY: return BObolOverlayOrder::XRay;
	case GED_VIEW_FEATURE_OVERLAY_ORDER_POST_TRANSPARENT:
	case GED_VIEW_FEATURE_OVERLAY_ORDER_UNKNOWN:
	default: return BObolOverlayOrder::PostTransparent;
    }
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
    scene->overlay_class = ged_obol_feature_batch_overlay_class(
	desc ? desc->overlay_class : GED_VIEW_FEATURE_OVERLAY_CLASS_COMMAND_RESULT);
    scene->lifecycle = ged_obol_feature_batch_lifecycle(
	desc ? desc->lifecycle : GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND);
    scene->overlay_order = ged_obol_feature_batch_overlay_order(
	desc ? desc->overlay_order :
	GED_VIEW_FEATURE_OVERLAY_ORDER_POST_TRANSPARENT);
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
    const struct ged_view_feature_line_layer *layers,
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
		std::vector<struct ged_view_feature_line_layer> layers_copy;
		layers_copy.reserve(staged_layers.size());
		for (const ged_view_feature_staged_line_layer &staged :
		    staged_layers) {
		    struct ged_view_feature_line_layer layer =
			ged_view_feature_line_layer_default();
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
ged_view_feature_batch_labels_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const struct ged_view_feature_label *labels,
    size_t label_count,
    const struct ged_view_feature_style *style)
{
    if (!ged_obol_feature_batch_valid(scene) || !name || !name[0] ||
	(label_count && !labels)) {
	if (ged_obol_feature_batch_valid(scene))
	    scene->failed = 1;
	return 0;
    }

    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name, "labelsReplace"))
	    return 0;
	const std::string staged_name(name);
	const ged_view_feature_staged_style staged_style =
	    ged_view_feature_stage_style(style);
	std::vector<ged_view_feature_staged_label> staged_labels;
	staged_labels.reserve(label_count);
	for (size_t i = 0; i < label_count; i++) {
	    ged_view_feature_staged_label staged;
	    staged.text = labels[i].text ? labels[i].text : "";
	    VMOVE(staged.point, labels[i].point);
	    staged.color_valid = labels[i].color_valid;
	    staged.color[0] = labels[i].color[0];
	    staged.color[1] = labels[i].color[1];
	    staged.color[2] = labels[i].color[2];
	    staged.line_flag = labels[i].line_flag;
	    VMOVE(staged.target, labels[i].target);
	    staged.anchor = labels[i].anchor;
	    staged.arrow = labels[i].arrow;
	    staged.font_size = labels[i].font_size;
	    staged_labels.push_back(staged);
	}
	scene->operations.emplace_back(
	    [staged_name, staged_style, staged_labels](
		struct ged_view_feature_batch *batch) {
		std::vector<struct ged_view_feature_label> copied;
		copied.reserve(staged_labels.size());
		for (const ged_view_feature_staged_label &staged : staged_labels) {
		    struct ged_view_feature_label label = {};
		    label.text = staged.text.c_str();
		    VMOVE(label.point, staged.point);
		    label.color_valid = staged.color_valid;
		    label.color[0] = staged.color[0];
		    label.color[1] = staged.color[1];
		    label.color[2] = staged.color[2];
		    label.line_flag = staged.line_flag;
		    VMOVE(label.target, staged.target);
		    label.anchor = staged.anchor;
		    label.arrow = staged.arrow;
		    label.font_size = staged.font_size;
		    copied.push_back(label);
		}
		return ged_view_feature_batch_labels_replace(batch,
		    staged_name.c_str(), copied.empty() ? NULL : copied.data(),
		    copied.size(), staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }

    if (!ged_obol_feature_batch_writable(scene, name, "labelsReplace"))
	return 0;
    if (!label_count) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
	    "labelsReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    std::vector<BObolLabel> obol_labels;
    obol_labels.reserve(label_count);
    for (size_t i = 0; i < label_count; i++)
	obol_labels.push_back(ged_view_feature_batch_label_to_obol(labels[i]));
    BObolFeatureStyle obol_style = ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle = scene->controller->features().publishLabels(
	name, scene->scope, obol_labels, style ? &obol_style : NULL,
	&scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_FAILED,
	    name, "labelsReplace", "feature publish failed");
	return 0;
    }
    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
	ged_obol_feature_batch_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_UPDATED,
	name, "labelsReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_arrow_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_view_feature_style *style)
{
    if (!ged_obol_feature_batch_valid(scene) || !name || !name[0] ||
	(point_count && !points)) {
	if (ged_obol_feature_batch_valid(scene))
	    scene->failed = 1;
	return 0;
    }
    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name, "arrowReplace"))
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
		return ged_view_feature_batch_arrow_replace(batch,
		    staged_name.c_str(), staged_points.empty() ? NULL :
		    reinterpret_cast<const point_t *>(staged_points.data()),
		    staged_points.size() / 3,
		    staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }
    if (!ged_obol_feature_batch_writable(scene, name, "arrowReplace"))
	return 0;
    if (!point_count) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
	    "arrowReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }
    BObolFeatureStyle obol_style = ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle = scene->controller->features().publishArrow(
	name, scene->scope, ged_obol_points_from_ged(points, point_count),
	style ? &obol_style : NULL, &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_FAILED,
	    name, "arrowReplace", "feature publish failed");
	return 0;
    }
    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
	ged_obol_feature_batch_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_UPDATED,
	name, "arrowReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_axes_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const point_t *centers,
    size_t center_count,
    fastf_t half_size,
    const struct ged_view_feature_style *style)
{
    if (!ged_obol_feature_batch_valid(scene) || !name || !name[0] ||
	(center_count && !centers) || half_size < 0.0) {
	if (ged_obol_feature_batch_valid(scene))
	    scene->failed = 1;
	return 0;
    }
    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name, "axesReplace"))
	    return 0;
	const std::string staged_name(name);
	const ged_view_feature_staged_style staged_style =
	    ged_view_feature_stage_style(style);
	const std::vector<fastf_t> staged_centers =
	    ged_view_feature_stage_vectors(
		reinterpret_cast<const fastf_t *>(centers), center_count);
	scene->operations.emplace_back(
	    [staged_name, staged_style, staged_centers, half_size](
		struct ged_view_feature_batch *batch) {
		return ged_view_feature_batch_axes_replace(batch,
		    staged_name.c_str(), staged_centers.empty() ? NULL :
		    reinterpret_cast<const point_t *>(staged_centers.data()),
		    staged_centers.size() / 3, half_size,
		    staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }
    if (!ged_obol_feature_batch_writable(scene, name, "axesReplace"))
	return 0;
    if (!center_count) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
	    "axesReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }
    BObolFeatureStyle obol_style = ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle = scene->controller->features().publishAxes(
	name, scene->scope, ged_obol_points_from_ged(centers, center_count),
	static_cast<float>(half_size), style ? &obol_style : NULL,
	&scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_FAILED,
	    name, "axesReplace", "feature publish failed");
	return 0;
    }
    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
	ged_obol_feature_batch_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_UPDATED,
	name, "axesReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_hud_labels_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const struct ged_view_feature_label *labels,
    size_t label_count,
    const struct ged_view_feature_style *style)
{
    if (!ged_obol_feature_batch_valid(scene) || !name || !name[0] ||
	(label_count && !labels)) {
	if (ged_obol_feature_batch_valid(scene))
	    scene->failed = 1;
	return 0;
    }
    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name, "hudLabelsReplace"))
	    return 0;
	const std::string staged_name(name);
	const ged_view_feature_staged_style staged_style =
	    ged_view_feature_stage_style(style);
	std::vector<ged_view_feature_staged_label> staged_labels;
	staged_labels.reserve(label_count);
	for (size_t i = 0; i < label_count; i++) {
	    ged_view_feature_staged_label staged;
	    staged.text = labels[i].text ? labels[i].text : "";
	    VMOVE(staged.point, labels[i].point);
	    staged.color_valid = labels[i].color_valid;
	    VMOVE(staged.color, labels[i].color);
	    staged.line_flag = labels[i].line_flag;
	    VMOVE(staged.target, labels[i].target);
	    staged.anchor = labels[i].anchor;
	    staged.arrow = labels[i].arrow;
	    staged.font_size = labels[i].font_size;
	    staged_labels.push_back(staged);
	}
	scene->operations.emplace_back(
	    [staged_name, staged_style, staged_labels](
		struct ged_view_feature_batch *batch) {
		std::vector<struct ged_view_feature_label> copied;
		copied.reserve(staged_labels.size());
		for (const ged_view_feature_staged_label &staged : staged_labels) {
		    struct ged_view_feature_label label = {};
		    label.text = staged.text.c_str();
		    VMOVE(label.point, staged.point);
		    label.color_valid = staged.color_valid;
		    VMOVE(label.color, staged.color);
		    label.line_flag = staged.line_flag;
		    VMOVE(label.target, staged.target);
		    label.anchor = staged.anchor;
		    label.arrow = staged.arrow;
		    label.font_size = staged.font_size;
		    copied.push_back(label);
		}
		return ged_view_feature_batch_hud_labels_replace(batch,
		    staged_name.c_str(), copied.empty() ? NULL : copied.data(),
		    copied.size(), staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }
    if (!ged_obol_feature_batch_writable(scene, name, "hudLabelsReplace"))
	return 0;
    if (!label_count) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
	    "hudLabelsReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }
    std::vector<BObolLabel> obol_labels;
    obol_labels.reserve(label_count);
    for (size_t i = 0; i < label_count; i++) {
	BObolLabel label = ged_view_feature_batch_label_to_obol(labels[i]);
	label.point = ged_view_feature_batch_hud_point(scene->view_ctx,
	    labels[i].point);
	if (labels[i].line_flag)
	    label.target = ged_view_feature_batch_hud_point(scene->view_ctx,
		labels[i].target);
	obol_labels.push_back(label);
    }
    BObolFeatureStyle obol_style = ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle = scene->controller->features().publishHudLabels(
	name, scene->scope, obol_labels, style ? &obol_style : NULL,
	&scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	return 0;
    }
    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
	ged_view_feature_batch_hud_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_UPDATED,
	name, "hudLabelsReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_hud_line_set_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style)
{
    if (!ged_obol_feature_batch_valid(scene) || !name || !name[0] ||
	(point_count && !points)) {
	if (ged_obol_feature_batch_valid(scene))
	    scene->failed = 1;
	return 0;
    }
    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name, "hudLineSetReplace"))
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
		return ged_view_feature_batch_hud_line_set_replace(batch,
		    staged_name.c_str(), staged_points.empty() ? NULL :
		    reinterpret_cast<const point_t *>(staged_points.data()),
		    staged_commands.empty() ? NULL : staged_commands.data(),
		    staged_points.size() / 3,
		    staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }
    if (!ged_obol_feature_batch_writable(scene, name, "hudLineSetReplace"))
	return 0;
    if (!point_count) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
	    "hudLineSetReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }
    std::vector<SbVec3f> screen_points;
    screen_points.reserve(point_count);
    for (size_t i = 0; i < point_count; i++)
	screen_points.push_back(ged_view_feature_batch_hud_point(
	    scene->view_ctx, points[i]));
    BObolFeatureStyle obol_style = ged_obol_feature_style_from_ged(style);
    obol_style.hud = TRUE;
    BObolFeatureHandle handle = scene->controller->features().publishLineSet(
	name, scene->scope, screen_points,
	ged_obol_commands_from_ged(cmds, point_count), &obol_style,
	&scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	return 0;
    }
    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
	ged_view_feature_batch_hud_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_UPDATED,
	name, "hudLineSetReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_hud_line_layers_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const struct ged_view_feature_line_layer *layers,
    size_t layer_count,
    const struct ged_view_feature_style *style)
{
    if (!ged_obol_feature_batch_valid(scene) || !name || !name[0]) {
	if (ged_obol_feature_batch_valid(scene))
	    scene->failed = 1;
	return 0;
    }
    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name,
		"hudLineLayersReplace"))
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
		std::vector<struct ged_view_feature_line_layer> copied;
		copied.reserve(staged_layers.size());
		for (const ged_view_feature_staged_line_layer &staged :
		    staged_layers) {
		    struct ged_view_feature_line_layer layer =
			ged_view_feature_line_layer_default();
		    layer.name = staged.name.empty() ? NULL : staged.name.c_str();
		    layer.points = staged.points.empty() ? NULL :
			reinterpret_cast<const point_t *>(staged.points.data());
		    layer.commands = staged.commands.empty() ? NULL :
			staged.commands.data();
		    layer.point_count = staged.points.size() / 3;
		    layer.style = staged.style;
		    copied.push_back(layer);
		}
		return ged_view_feature_batch_hud_line_layers_replace(batch,
		    staged_name.c_str(), copied.empty() ? NULL : copied.data(),
		    copied.size(), staged_style.valid ? &staged_style.value : NULL);
	    });
	return 1;
    }
    if (!ged_obol_feature_batch_writable(scene, name,
	    "hudLineLayersReplace"))
	return 0;
    std::vector<BObolLineLayer> obol_layers;
    for (size_t i = 0; layers && i < layer_count; i++) {
	if (!layers[i].points || !layers[i].point_count)
	    continue;
	BObolLineLayer layer;
	layer.name = layers[i].name ? layers[i].name : name;
	layer.style = ged_obol_feature_style_from_ged(&layers[i].style);
	layer.style.hud = TRUE;
	for (size_t j = 0; j < layers[i].point_count; j++) {
	    layer.points.push_back(ged_view_feature_batch_hud_point(
		scene->view_ctx, layers[i].points[j]));
	    layer.commands.push_back(ged_obol_line_command_from_ged(
		layers[i].commands ? layers[i].commands[j] : -1, j));
	}
	obol_layers.push_back(std::move(layer));
    }
    if (obol_layers.empty()) {
	(void)ged_obol_feature_batch_remove_feature(scene, name,
	    "hudLineLayersReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }
    BObolFeatureStyle obol_style = ged_obol_feature_style_from_ged(style);
    obol_style.hud = TRUE;
    BObolFeatureHandle handle = scene->controller->features().publishLineLayers(
	name, scene->scope, obol_layers, &obol_style, &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	return 0;
    }
    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
	ged_view_feature_batch_hud_overlay_info(scene, name));
    scene->changed++;
    ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_UPDATED,
	name, "hudLineLayersReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_view_feature_batch_hud_axes_replace(
    struct ged_view_feature_batch *scene,
    const char *name,
    const struct bv_axes_state *axes,
    const mat_t rotation)
{
    if (!ged_obol_feature_batch_valid(scene) || !name || !name[0]) {
	if (ged_obol_feature_batch_valid(scene))
	    scene->failed = 1;
	return 0;
    }
    if (!scene->committing) {
	if (!ged_obol_feature_batch_writable(scene, name, "hudAxesReplace"))
	    return 0;
	const std::string staged_name(name);
	const bool has_axes = axes != NULL;
	struct bv_axes_state staged_axes = BV_AXES_STATE_INIT;
	if (axes)
	    staged_axes = *axes;
	std::array<fastf_t, 16> staged_rotation = {};
	const bool has_rotation = rotation != NULL;
	if (rotation)
	    std::copy(rotation, rotation + 16, staged_rotation.begin());
	scene->operations.emplace_back(
	    [staged_name, has_axes, staged_axes, has_rotation, staged_rotation](
		struct ged_view_feature_batch *batch) {
		mat_t rotation_copy;
		if (has_rotation)
		    std::copy(staged_rotation.begin(), staged_rotation.end(),
			rotation_copy);
		return ged_view_feature_batch_hud_axes_replace(batch,
		    staged_name.c_str(), has_axes ? &staged_axes : NULL,
		    has_rotation ? rotation_copy : NULL);
	    });
	return 1;
    }
    if (!ged_draw_obol_view_context_hud_axes_replace(scene->view_ctx, name,
	    axes, rotation)) {
	scene->failed = 1;
	return 0;
    }
    scene->changed++;
    ged_obol_feature_batch_notify(scene, GED_VIEW_FEATURE_BATCH_UPDATED,
	name, "hudAxesReplace", NULL);
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
			ged_view_feature_metadata_default();
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
			ged_view_feature_metadata_default();
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
