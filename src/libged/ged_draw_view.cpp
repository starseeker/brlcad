/*                  G E D _ D R A W _ V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file ged_draw_view.cpp
 *
 * C++ GED view-feature integration backed exclusively by Obol/libBObol.
 *
 * Absent Obol view state is treated as a no-op/failure condition, not as
 * permission to synthesize alternate retained records.
 */

#include "common.h"

#include <cfloat>
#include <string.h>

#include "bu/color.h"
#include "bu/malloc.h"
#include "bu/ptbl.h"
#include "bu/vls.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "ged/event.h"
#include "rt/db_io.h"
#include "./ged_bobol_private.hpp"
#include "./ged_draw_private.h"
#include "./ged_draw_view_private.h"

struct ged_pick_record {
    struct bu_vls path;
    fastf_t hit_dist;
    struct ged_pick_detail detail;
};

struct ged_pick_result {
    struct ged_pick_record *records;
    size_t count;
    size_t capacity;
};

struct ged_bobol_feature_lookup {
    BObolViewController *controller = nullptr;
    BObolFeatureHandle handle;
};


void
ged_view_feature_style_init(struct ged_view_feature_style *style)
{
    if (!style)
	return;
    memset(style, 0, sizeof(*style));
    style->visible = -1;
    style->selectable = -1;
    style->line_width = -1;
    style->line_style = -1;
    style->arrow = -1;
    style->arrow_tip_length = -1.0;
    style->arrow_tip_width = -1.0;
}


void
ged_view_feature_summary_init(struct ged_view_feature_summary *summary)
{
    if (!summary)
	return;
    memset(summary, 0, sizeof(*summary));
    summary->kind = GED_VIEW_FEATURE_KIND_UNKNOWN;
    summary->scope = GED_VIEW_FEATURE_SCOPE_UNKNOWN;
    summary->overlay_class = GED_VIEW_FEATURE_OVERLAY_CLASS_UNKNOWN;
    summary->lifecycle = GED_VIEW_FEATURE_LIFECYCLE_UNKNOWN;
}


void
ged_pick_detail_init(struct ged_pick_detail *detail)
{
    if (!detail)
	return;
    memset(detail, 0, sizeof(*detail));
    detail->primitive_index = -1;
    detail->face_vertex_index[0] = -1;
    detail->face_vertex_index[1] = -1;
    detail->face_vertex_index[2] = -1;
    detail->nearest_face_vertex_index = -1;
}


void
ged_view_edit_transaction_init(struct ged_view_edit_transaction *transaction)
{
    if (!transaction)
	return;
    memset(transaction, 0, sizeof(*transaction));
    transaction->event = GED_VIEW_EDIT_PREVIEW_UPDATE;
    transaction->feature = GED_VIEW_EDIT_REF_NULL;
    transaction->color[0] = 255;
    transaction->color[1] = 255;
    transaction->color[2] = 255;
}

static BObolFeatureOwner
ged_bobol_feature_owner(struct ged_view_context *view_ctx)
{
    BObolFeatureOwner owner;
    owner.ownerToken = view_ctx;
    char owner_id[64] = {0};
    snprintf(owner_id, sizeof(owner_id), "%p", static_cast<void *>(view_ctx));
    owner.ownerId = owner_id;
    owner.ownerRole = "view";
    return owner;
}

static BObolViewController *
ged_bobol_shared_feature_controller(struct ged_view_context *view_ctx)
{
    struct ged *gedp = ged_view_context_owner(view_ctx);
    return gedp ? ged_draw_obol_controller(gedp) : nullptr;
}

static struct ged_bobol_feature_lookup
ged_bobol_feature_find(struct ged_view_context *view_ctx, const char *name)
{
    struct ged_bobol_feature_lookup lookup;
    if (!name || !name[0])
	return lookup;

    BObolViewController *view_controller =
	ged_bobol_view_controller(view_ctx);
    if (view_controller) {
	const BObolFeatureOwner owner = ged_bobol_feature_owner(view_ctx);
	lookup.handle = view_controller->features().findOwned(name,
	    BOBOL_FEATURE_SCOPE_LOCAL, &owner);
	/* Public named feature queries address the retained object in this view,
	 * independent of which command producer owns its lifecycle.  Requiring
	 * the synthetic view owner here made batch-published local features
	 * readable by the export adapter but invisible to the public style and
	 * existence APIs. */
	if (!lookup.handle.isValid())
	    lookup.handle = view_controller->features().find(name,
		BOBOL_FEATURE_SCOPE_LOCAL);
	if (!lookup.handle.isValid())
	    lookup.handle = view_controller->features().find(name,
		BOBOL_FEATURE_SCOPE_SHARED);
	if (lookup.handle.isValid()) {
	    lookup.controller = view_controller;
	    return lookup;
	}
    }

    BObolViewController *shared_controller =
	ged_bobol_shared_feature_controller(view_ctx);
    if (!shared_controller || shared_controller == view_controller)
	return lookup;
    lookup.handle = shared_controller->features().find(name,
	BOBOL_FEATURE_SCOPE_SHARED);
    if (lookup.handle.isValid())
	lookup.controller = shared_controller;
    return lookup;
}

static unsigned char
ged_bobol_rgb_channel(float value)
{
    if (value <= 0.0f)
	return 0;
    if (value >= 1.0f)
	return 255;
    return static_cast<unsigned char>(value * 255.0f + 0.5f);
}

static void
ged_bobol_rgb_from_color(const SbColor &color, unsigned char rgb[3])
{
    rgb[0] = ged_bobol_rgb_channel(color[0]);
    rgb[1] = ged_bobol_rgb_channel(color[1]);
    rgb[2] = ged_bobol_rgb_channel(color[2]);
}

static SbColor
ged_bobol_color_from_rgb(const unsigned char rgb[3])
{
    return SbColor(static_cast<float>(rgb[0]) / 255.0f,
	static_cast<float>(rgb[1]) / 255.0f,
	static_cast<float>(rgb[2]) / 255.0f);
}

static BObolFeatureStyle
ged_bobol_feature_style_from_ged(const struct ged_view_feature_style *style)
{
    BObolFeatureStyle out;
    if (!style)
	return out;
    if (style->visible != -1) {
	out.hasVisible = TRUE;
	out.visible = style->visible ? TRUE : FALSE;
    }
    if (style->selectable != -1) {
	out.hasSelectable = TRUE;
	out.selectable = style->selectable ? TRUE : FALSE;
    }
    if (style->color_valid) {
	out.hasColor = TRUE;
	out.color = ged_bobol_color_from_rgb(style->color);
    }
    if (style->line_width >= 0) {
	out.hasLineWidth = TRUE;
	out.lineWidth = style->line_width;
    }
    if (style->line_style >= 0) {
	out.hasLineStyle = TRUE;
	out.lineStyle = style->line_style;
    }
    if (style->arrow != -1) {
	out.hasArrow = TRUE;
	out.arrow = style->arrow ? TRUE : FALSE;
    }
    if (style->arrow_tip_length >= 0.0 || style->arrow_tip_width >= 0.0) {
	out.hasArrowTip = TRUE;
	out.arrowTipLength = style->arrow_tip_length >= 0.0 ?
	    static_cast<float>(style->arrow_tip_length) : 0.0f;
	out.arrowTipWidth = style->arrow_tip_width >= 0.0 ?
	    static_cast<float>(style->arrow_tip_width) : 0.0f;
    }
    return out;
}

static void
ged_bobol_feature_style_to_ged(struct ged_view_feature_style *dst,
    const BObolFeatureStyle &src)
{
    struct ged_view_feature_style init = ged_view_feature_style_default();
    *dst = init;
    if (src.hasVisible)
	dst->visible = src.visible ? 1 : 0;
    if (src.hasSelectable)
	dst->selectable = src.selectable ? 1 : 0;
    if (src.hasColor) {
	dst->color_valid = 1;
	ged_bobol_rgb_from_color(src.color, dst->color);
    }
    if (src.hasLineWidth)
	dst->line_width = src.lineWidth;
    if (src.hasLineStyle)
	dst->line_style = src.lineStyle;
    if (src.hasArrow)
	dst->arrow = src.arrow ? 1 : 0;
    if (src.hasArrowTip) {
	dst->arrow_tip_length = src.arrowTipLength;
	dst->arrow_tip_width = src.arrowTipWidth;
    }
}

static int
ged_bobol_feature_kind_to_ged(BObolFeatureKind kind)
{
    switch (kind) {
	case BObolFeatureKind::Lines:
	    return GED_VIEW_FEATURE_KIND_LINES;
	case BObolFeatureKind::IndexedLines:
	    return GED_VIEW_FEATURE_KIND_INDEXED_LINES;
	case BObolFeatureKind::Points:
	    return GED_VIEW_FEATURE_KIND_POINTS;
	case BObolFeatureKind::Labels:
	    return GED_VIEW_FEATURE_KIND_LABELS;
	case BObolFeatureKind::Arrow:
	    return GED_VIEW_FEATURE_KIND_ARROW;
	case BObolFeatureKind::Axes:
	    return GED_VIEW_FEATURE_KIND_AXES;
	case BObolFeatureKind::LineLayer:
	    return GED_VIEW_FEATURE_KIND_LINE_LAYER;
	case BObolFeatureKind::EditPreview:
	    return GED_VIEW_FEATURE_KIND_EDIT_PREVIEW;
	case BObolFeatureKind::IndexedFaceSet:
	    return GED_VIEW_FEATURE_KIND_INDEXED_FACE_SET;
	case BObolFeatureKind::PolygonOverlay:
	    return GED_VIEW_FEATURE_KIND_POLYGON_OVERLAY;
	case BObolFeatureKind::HudLabel:
	    return GED_VIEW_FEATURE_KIND_HUD_LABEL;
	case BObolFeatureKind::CustomNode:
	    return GED_VIEW_FEATURE_KIND_CUSTOM_NODE;
	case BObolFeatureKind::Unknown:
	default:
	    return GED_VIEW_FEATURE_KIND_UNKNOWN;
    }
}

static int
ged_bobol_feature_scope_to_ged(BObolFeatureScope scope)
{
    return scope == BObolFeatureScope::Local ?
	GED_VIEW_FEATURE_SCOPE_LOCAL : GED_VIEW_FEATURE_SCOPE_SHARED;
}

static int
ged_bobol_overlay_class_to_ged(BObolOverlayClass overlay_class)
{
    switch (overlay_class) {
	case BObolOverlayClass::None:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_NONE;
	case BObolOverlayClass::Faceplate:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_FACEPLATE;
	case BObolOverlayClass::EditHandle:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_EDIT_HANDLE;
	case BObolOverlayClass::Measure:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_MEASURE;
	case BObolOverlayClass::SelectionRubberBand:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_SELECTION_RUBBER_BAND;
	case BObolOverlayClass::SnapGuide:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_SNAP_GUIDE;
	case BObolOverlayClass::CommandResult:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_COMMAND_RESULT;
	case BObolOverlayClass::Diagnostic:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_DIAGNOSTIC;
	case BObolOverlayClass::TclOverlay:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_TCL_OVERLAY;
	case BObolOverlayClass::PolygonEdit:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_POLYGON_EDIT;
	case BObolOverlayClass::SketchEdit:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_SKETCH_EDIT;
	case BObolOverlayClass::UserAnnotation:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION;
	default:
	    return GED_VIEW_FEATURE_OVERLAY_CLASS_UNKNOWN;
    }
}

static int
ged_bobol_overlay_lifecycle_to_ged(BObolOverlayLifecycle lifecycle)
{
    switch (lifecycle) {
	case BObolOverlayLifecycle::None:
	    return GED_VIEW_FEATURE_LIFECYCLE_NONE;
	case BObolOverlayLifecycle::Persistent:
	    return GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT;
	case BObolOverlayLifecycle::PerFrame:
	    return GED_VIEW_FEATURE_LIFECYCLE_PER_FRAME;
	case BObolOverlayLifecycle::PerCommand:
	    return GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND;
	case BObolOverlayLifecycle::PerTool:
	    return GED_VIEW_FEATURE_LIFECYCLE_PER_TOOL;
	case BObolOverlayLifecycle::PerView:
	    return GED_VIEW_FEATURE_LIFECYCLE_PER_VIEW;
	case BObolOverlayLifecycle::SharedViewSet:
	    return GED_VIEW_FEATURE_LIFECYCLE_SHARED_VIEW_SET;
	case BObolOverlayLifecycle::AutoRemoveOnSource:
	    return GED_VIEW_FEATURE_LIFECYCLE_AUTO_REMOVE_ON_SOURCE;
	default:
	    return GED_VIEW_FEATURE_LIFECYCLE_UNKNOWN;
    }
}

static int
ged_pick_result_reserve(struct ged_pick_result *result,
			     size_t capacity)
{
    if (!result)
	return 0;
    if (capacity <= result->capacity)
	return 1;

    size_t new_capacity = result->capacity ? result->capacity : 4;
    while (new_capacity < capacity)
	new_capacity *= 2;

    result->records = (struct ged_pick_record *)bu_realloc(
	    result->records, new_capacity * sizeof(struct ged_pick_record),
	    "ged draw pick result records");
    for (size_t i = result->capacity; i < new_capacity; i++) {
	BU_VLS_INIT(&result->records[i].path);
	result->records[i].hit_dist = 0.0;
	result->records[i].detail = ged_pick_detail_default();
    }
    result->capacity = new_capacity;
    return 1;
}

struct ged_pick_result *
ged_pick_result_create(void)
{
    struct ged_pick_result *result;
    BU_ALLOC(result, struct ged_pick_result);
    result->records = NULL;
    result->count = 0;
    result->capacity = 0;
    return result;
}

void
ged_pick_result_free(struct ged_pick_result *result)
{
    if (!result)
	return;
    for (size_t i = 0; i < result->capacity; i++)
	bu_vls_free(&result->records[i].path);
    bu_free(result->records, "ged draw pick result records");
    bu_free(result, "ged draw pick result");
}

size_t
ged_pick_result_count(const struct ged_pick_result *result)
{
    return result ? result->count : 0;
}

int
ged_pick_result_path(const struct ged_pick_result *result,
			  size_t index,
			  struct bu_vls *path_out)
{
    if (!result || !path_out || index >= result->count)
	return 0;
    bu_vls_strcpy(path_out, bu_vls_cstr(&result->records[index].path));
    return 1;
}

fastf_t
ged_pick_result_hit_dist(const struct ged_pick_result *result,
			      size_t index)
{
    if (!result || index >= result->count)
	return 0.0;
    return result->records[index].hit_dist;
}

int
ged_pick_result_detail(const struct ged_pick_result *result,
	size_t index, struct ged_pick_detail *detail)
{
    if (!result || !detail || index >= result->count)
	return 0;
    *detail = result->records[index].detail;
    return 1;
}

int
ged_pick_result_append_detail(struct ged_pick_result *result,
	const char *path, fastf_t hit_dist,
	const struct ged_pick_detail *detail)
{
    if (!result || !path || !path[0])
	return 0;
    if (!ged_pick_result_reserve(result, result->count + 1))
	return 0;

    struct ged_pick_record *record = &result->records[result->count++];
    bu_vls_strcpy(&record->path, path);
    record->hit_dist = hit_dist;
    record->detail = detail ? *detail : ged_pick_detail_default();
    return 1;
}

int
ged_pick_result_append_path(struct ged_pick_result *result,
				 const char *path,
				 fastf_t hit_dist)
{
    return ged_pick_result_append_detail(result, path, hit_dist, NULL);
}

int
ged_pick_result_append_copy(struct ged_pick_result *dest,
				 const struct ged_pick_result *src,
				 size_t index,
				 fastf_t hit_dist)
{
    if (!dest || !src || index >= src->count)
	return 0;
    return ged_pick_result_append_detail(dest,
	    bu_vls_cstr(&src->records[index].path), hit_dist,
	    &src->records[index].detail);
}

struct ged_pick_result *
ged_pick_result_filter_first(const struct ged_pick_result *src)
{
    struct ged_pick_result *result = ged_pick_result_create();
    if (!result)
	return NULL;
    if (src && src->count > 0)
	(void)ged_pick_result_append_copy(result, src, 0,
		src->records[0].hit_dist);
    return result;
}

struct ged_pick_result *
ged_pick_point(struct ged_view_context *view_ctx, int x, int y, int first_only)
{
    return ged_draw_obol_view_context_pick_point(view_ctx, x, y, first_only);
}

struct ged_pick_result *
ged_pick_nearest(struct ged_view_context *view_ctx, int x, int y)
{
    return ged_draw_obol_view_context_pick_nearest(view_ctx, x, y);
}

struct ged_pick_result *
ged_pick_rect(struct ged_view_context *view_ctx, int x0, int y0, int x1, int y1)
{
    return ged_draw_obol_view_context_pick_rect(view_ctx, x0, y0, x1, y1);
}

int
ged_selection_available(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_view_context_selection_available(view_ctx);
}


size_t
ged_view_selection_count(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_view_context_selection_count(view_ctx);
}


int
ged_view_selection_visit(
	struct ged_view_context *view_ctx,
	ged_view_selection_visit_cb cb,
	void *data)
{
    return ged_draw_obol_view_context_selection_path_foreach(view_ctx, cb,
	    data);
}


int
ged_view_selection_clear(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_view_context_selection_clear(view_ctx);
}

int
ged_view_selection_set_pick(
	struct ged_view_context *view_ctx,
	const struct ged_pick_result *result,
	ged_view_selection_path_cb cb,
	void *data)
{
    if (!view_ctx)
	return 0;
    if (!ged_view_selection_clear(view_ctx))
	return 0;
    if (!result)
	return 1;

    int ok = 1;
    for (size_t i = 0; i < result->count; i++) {
	const char *path = bu_vls_cstr(&result->records[i].path);
	if (!path || !path[0])
	    continue;
	if (!ged_selection_add_path(view_ctx,
		GED_SELECTION_SELECTED_PATH, path))
	    ok = 0;
	if (cb)
	    cb(path, data);
    }
    return ok;
}


int
ged_selection_contains_path(
	struct ged_view_context *view_ctx,
	enum ged_selection_kind kind,
	const char *path)
{
    return ged_draw_obol_view_context_selection_contains_path(view_ctx,
	    (int)kind, path);
}


int
ged_selection_add_path(
	struct ged_view_context *view_ctx,
	enum ged_selection_kind kind,
	const char *path)
{
    return ged_draw_obol_view_context_selection_add_path(view_ctx,
	    (int)kind, path);
}


int
ged_selection_remove_path(
	struct ged_view_context *view_ctx,
	enum ged_selection_kind kind,
	const char *path)
{
    return ged_draw_obol_view_context_selection_remove_path(view_ctx,
	    (int)kind, path);
}


int
ged_selection_set_path(
	struct ged_view_context *view_ctx,
	enum ged_selection_kind kind,
	const char *path)
{
    return ged_draw_obol_view_context_selection_set_path(view_ctx,
	    (int)kind, path);
}


int
ged_view_feature_exists(struct ged_view_context *view_ctx, const char *name)
{
    return ged_bobol_feature_find(view_ctx, name).handle.isValid() ? 1 : 0;
}


int
ged_view_feature_remove(struct ged_view_context *view_ctx, const char *name)
{
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    return lookup.controller && lookup.handle.isValid() &&
	lookup.controller->features().remove(lookup.handle) ? 1 : 0;
}


int
ged_view_edit_remove_ref(struct ged_view_context *view_ctx,
	ged_view_edit_ref ref)
{
    return ged_draw_obol_view_feature_remove_ref(view_ctx, ref);
}


int
ged_view_feature_get_summary(
	struct ged_view_context *view_ctx,
	const char *name,
	struct ged_view_feature_summary *summary)
{
    if (!summary)
	return 0;
    memset(summary, 0, sizeof(*summary));
    if (!name || !name[0])
	return 0;

    BObolFeatureSummary obol_summary;
    BObolViewController *view_controller =
	ged_bobol_view_controller(view_ctx);
    if (view_controller) {
	const BObolFeatureOwner owner = ged_bobol_feature_owner(view_ctx);
	if (!view_controller->features().summaryOwned(name, obol_summary,
		BOBOL_FEATURE_SCOPE_LOCAL, &owner))
	    return 0;
    }
    if (!obol_summary.exists) {
	BObolViewController *shared_controller =
	    ged_bobol_shared_feature_controller(view_ctx);
	if (!shared_controller || !shared_controller->features().summary(name,
		obol_summary, BOBOL_FEATURE_SCOPE_SHARED))
	    return 0;
    }
    if (!obol_summary.exists)
	return 1;

    summary->exists = 1;
    summary->visible = obol_summary.visible ? 1 : 0;
    summary->kind = ged_bobol_feature_kind_to_ged(obol_summary.kind);
    summary->scope = ged_bobol_feature_scope_to_ged(obol_summary.scope);
    summary->overlay_class =
	ged_bobol_overlay_class_to_ged(obol_summary.overlay.overlayClass);
    summary->lifecycle =
	ged_bobol_overlay_lifecycle_to_ged(obol_summary.overlay.lifecycle);
    summary->child_count = obol_summary.childCount;
    summary->geometry_command_count = obol_summary.commandCount;
    summary->metadata_count = obol_summary.metadataCount;
    summary->primitive_metadata_count = obol_summary.primitiveMetadataCount;
    summary->selected_primitive_count = obol_summary.selectedPrimitiveCount;
    summary->highlighted_primitive_count =
	obol_summary.highlightedPrimitiveCount;
    summary->is_label = obol_summary.kind == BObolFeatureKind::Labels ||
	obol_summary.kind == BObolFeatureKind::HudLabel;
    summary->is_transient_preview =
	obol_summary.kind == BObolFeatureKind::EditPreview;
    summary->is_command_result = obol_summary.overlay.overlayClass ==
	BObolOverlayClass::CommandResult;
    summary->is_overlay = obol_summary.overlay.isOverlay ||
	(!summary->is_label && !summary->is_transient_preview);
    summary->owner_generation = obol_summary.owner.generation;
    if (obol_summary.owner.ownerId.getLength() > 0)
	snprintf(summary->owner_id, sizeof(summary->owner_id), "%s",
	    obol_summary.owner.ownerId.getString());
    if (obol_summary.owner.ownerRole.getLength() > 0)
	snprintf(summary->owner_role, sizeof(summary->owner_role), "%s",
	    obol_summary.owner.ownerRole.getString());

    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    BObolFeatureStyle style;
    if (lookup.controller && lookup.handle.isValid() &&
	lookup.controller->features().style(lookup.handle, style) &&
	style.hasColor)
	ged_bobol_rgb_from_color(style.color, summary->color);
    return 1;
}


size_t
ged_view_feature_metadata_count(struct ged_view_context *view_ctx, const char *name)
{
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    std::vector<BObolFeatureMetadata> metadata;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().metadata(lookup.handle, metadata))
	return 0;
    return metadata.size();
}


int
ged_view_feature_metadata_copy(
	struct ged_view_context *view_ctx,
	const char *name,
	size_t index,
	struct bu_vls *key,
	struct bu_vls *value)
{
    if (key)
	bu_vls_trunc(key, 0);
    if (value)
	bu_vls_trunc(value, 0);
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    std::vector<BObolFeatureMetadata> metadata;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().metadata(lookup.handle, metadata) ||
	index >= metadata.size())
	return 0;
    if (key)
	bu_vls_strcat(key, metadata[index].key.getString());
    if (value)
	bu_vls_strcat(value, metadata[index].value.getString());
    return 1;
}


size_t
ged_view_feature_primitive_metadata_count(
	struct ged_view_context *view_ctx,
	const char *name,
	int primitive)
{
    if (primitive < 0)
	return 0;
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    std::vector<BObolFeatureMetadata> metadata;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().primitiveMetadata(lookup.handle,
	    static_cast<int32_t>(primitive), metadata))
	return 0;
    return metadata.size();
}


int
ged_view_feature_primitive_metadata_copy(
	struct ged_view_context *view_ctx,
	const char *name,
	int primitive,
	size_t index,
	struct bu_vls *key,
	struct bu_vls *value)
{
    if (key)
	bu_vls_trunc(key, 0);
    if (value)
	bu_vls_trunc(value, 0);
    if (primitive < 0)
	return 0;
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    std::vector<BObolFeatureMetadata> metadata;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().primitiveMetadata(lookup.handle,
	    static_cast<int32_t>(primitive), metadata) ||
	index >= metadata.size())
	return 0;
    if (key)
	bu_vls_strcat(key, metadata[index].key.getString());
    if (value)
	bu_vls_strcat(value, metadata[index].value.getString());
    return 1;
}


int
ged_view_feature_pick_primitive_resolve(
	struct ged_view_context *view_ctx,
	const char *picked_feature_name,
	int picked_primitive,
	int select,
	int highlight,
	struct bu_vls *feature_name,
	int *feature_primitive)
{
    if (feature_name)
	bu_vls_trunc(feature_name, 0);
    if (feature_primitive)
	*feature_primitive = -1;
    if (!picked_feature_name || picked_primitive < 0 || !feature_primitive)
	return 0;

    BObolFeaturePrimitivePick pick;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    const BObolFeatureOwner owner = ged_bobol_feature_owner(view_ctx);
    int resolved = controller ? controller->features().resolvePrimitivePick(
	picked_feature_name, static_cast<int32_t>(picked_primitive), pick,
	BOBOL_FEATURE_SCOPE_LOCAL, &owner) : 0;
    if (!resolved && controller)
	resolved = controller->features().resolvePrimitivePick(
	    picked_feature_name, static_cast<int32_t>(picked_primitive), pick,
	    BOBOL_FEATURE_SCOPE_SHARED, nullptr);
    if (!resolved) {
	BObolViewController *shared_controller =
	    ged_bobol_shared_feature_controller(view_ctx);
	controller = shared_controller != controller ? shared_controller : nullptr;
	resolved = controller ? controller->features().resolvePrimitivePick(
	    picked_feature_name, static_cast<int32_t>(picked_primitive), pick,
	    BOBOL_FEATURE_SCOPE_SHARED, nullptr) : 0;
    }
    if (!resolved || !controller)
	return 0;

    const std::vector<int32_t> primitives(1, pick.primitiveIndex);
    if (select && !controller->features().replaceSelectedPrimitives(
	    pick.handle, primitives))
	return 0;
    if (highlight && !controller->features().replaceHighlightedPrimitives(
	    pick.handle, primitives))
	return 0;
    if (feature_name)
	bu_vls_strcat(feature_name, pick.featureName.getString());
    *feature_primitive = static_cast<int>(pick.primitiveIndex);
    return 1;
}


int
ged_view_feature_set_selection(
	struct ged_view_context *view_ctx,
	const char *name,
	const int *primitives,
	size_t primitive_count)
{
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    if (!lookup.controller || !lookup.handle.isValid())
	return 0;
    std::vector<int32_t> values;
    values.reserve(primitive_count);
    for (size_t i = 0; primitives && i < primitive_count; i++)
	if (primitives[i] >= 0)
	    values.push_back(static_cast<int32_t>(primitives[i]));
    return lookup.controller->features().replaceSelectedPrimitives(
	lookup.handle, values) ? 1 : 0;
}


int
ged_view_feature_set_highlights(
	struct ged_view_context *view_ctx,
	const char *name,
	const int *primitives,
	size_t primitive_count)
{
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    if (!lookup.controller || !lookup.handle.isValid())
	return 0;
    std::vector<int32_t> values;
    values.reserve(primitive_count);
    for (size_t i = 0; primitives && i < primitive_count; i++)
	if (primitives[i] >= 0)
	    values.push_back(static_cast<int32_t>(primitives[i]));
    return lookup.controller->features().replaceHighlightedPrimitives(
	lookup.handle, values) ? 1 : 0;
}


size_t
ged_view_feature_selection_count(
	struct ged_view_context *view_ctx,
	const char *name)
{
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    std::vector<int32_t> primitives;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().selectedPrimitives(lookup.handle,
	    primitives))
	return 0;
    return primitives.size();
}


size_t
ged_view_feature_highlight_count(
	struct ged_view_context *view_ctx,
	const char *name)
{
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    std::vector<int32_t> primitives;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().highlightedPrimitives(lookup.handle,
	    primitives))
	return 0;
    return primitives.size();
}


int
ged_view_feature_selection_at(
	struct ged_view_context *view_ctx,
	const char *name,
	size_t index,
	int *primitive)
{
    if (primitive)
	*primitive = -1;
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    std::vector<int32_t> primitives;
    if (!primitive || !lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().selectedPrimitives(lookup.handle,
	    primitives) || index >= primitives.size())
	return 0;
    *primitive = static_cast<int>(primitives[index]);
    return 1;
}


int
ged_view_feature_highlight_at(
	struct ged_view_context *view_ctx,
	const char *name,
	size_t index,
	int *primitive)
{
    if (primitive)
	*primitive = -1;
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    std::vector<int32_t> primitives;
    if (!primitive || !lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().highlightedPrimitives(lookup.handle,
	    primitives) || index >= primitives.size())
	return 0;
    *primitive = static_cast<int>(primitives[index]);
    return 1;
}

int
ged_view_edit_ref_is_null(ged_view_edit_ref ref)
{
    return ged_draw_obol_view_feature_ref_is_null(ref);
}


int
ged_view_edit_preview_publish_event(
	struct ged_view_context *view_ctx,
	ged_view_edit_ref feature,
	enum ged_view_edit_preview_event event,
	const char *source_path)
{
    return ged_draw_obol_view_context_edit_preview_publish_event(view_ctx,
	    feature, event, source_path);
}

int
ged_view_edit_transaction_apply(
	struct ged_view_context *view_ctx,
	const struct ged_view_edit_transaction *transaction,
	ged_view_edit_ref *feature_out)
{
    if (feature_out)
	*feature_out = GED_VIEW_EDIT_REF_NULL;
    if (!view_ctx || !transaction || !transaction->feature_name ||
	!transaction->feature_name[0])
	return 0;

    ged_view_edit_ref feature = transaction->feature;
    const int terminal =
	transaction->event == GED_VIEW_EDIT_PREVIEW_COMMIT ||
	transaction->event == GED_VIEW_EDIT_PREVIEW_CANCEL ||
	transaction->event == GED_VIEW_EDIT_PREVIEW_DISCARD;
    if (terminal) {
	int ret = ged_view_edit_ref_is_null(feature) ?
	    ged_view_feature_remove(view_ctx,
		transaction->feature_name) :
	    ged_view_edit_preview_publish_event(view_ctx, feature,
		transaction->event, transaction->source_path);
	return ret;
    }

    if (ged_view_edit_ref_is_null(feature))
	feature = ged_view_edit_overlay_ensure(view_ctx,
	    transaction->feature_name, transaction->source_path);
    if (ged_view_edit_ref_is_null(feature))
	return 0;

    if (transaction->color_valid)
	    ged_view_edit_color_set(view_ctx, feature, transaction->color[0],
	    transaction->color[1], transaction->color[2]);

    int updated = 1;
    if (transaction->internal) {
	updated = ged_view_edit_primitive_wireframe_replace(view_ctx, feature,
	    transaction->dbip, transaction->internal, transaction->matrix,
	    transaction->ttol, transaction->tol);
    } else if (transaction->point_count > 0) {
	updated = ged_view_edit_preview_replace(view_ctx, feature,
	    transaction->source_path, transaction->edit_intent_id,
	    transaction->edit_intent_role, transaction->points,
	    transaction->commands, transaction->point_count,
	    transaction->source_revision, transaction->inputs_revision);
    }
    if (!updated)
	return 0;

    (void)ged_view_edit_preview_publish_event(view_ctx,
	feature, transaction->event, transaction->source_path);
    if (feature_out)
	*feature_out = feature;
    return 1;
}

int
ged_view_edit_transaction_apply_all(
	struct ged *gedp,
	const struct ged_view_edit_transaction *transaction)
{
    if (!gedp || !transaction)
	return 0;

    struct ged_view_edit_transaction local = *transaction;
    local.feature = GED_VIEW_EDIT_REF_NULL;
    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    struct ged_view_context *active_view = ged_view_active_ctx(gedp);
    int active_seen = 0;
    int applied = 0;

    if (views) {
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    struct ged_view_context *view_ctx =
		(struct ged_view_context *)BU_PTBL_GET(views, i);
	    if (!view_ctx)
		continue;
	    if (view_ctx == active_view)
		active_seen = 1;
	    if (ged_view_edit_transaction_apply(view_ctx,
		    &local, NULL))
		applied++;
	}
    }

    if (active_view && !active_seen &&
	ged_view_edit_transaction_apply(active_view, &local, NULL))
	applied++;
    return applied;
}


ged_view_edit_ref
ged_view_edit_overlay_ensure(
	struct ged_view_context *view_ctx,
	const char *name,
	const char *source_path)
{
    return ged_draw_obol_view_context_feature_overlay_ensure(view_ctx, name,
	    source_path);
}


ged_view_edit_ref
ged_view_edit_label_ensure(struct ged_view_context *view_ctx,
					   const char *name)
{
    return ged_draw_obol_view_context_feature_label_ensure(view_ctx, name);
}


int
ged_view_edit_visible_set(struct ged_view_context *view_ctx,
	ged_view_edit_ref ref, int visible)
{
    return ged_draw_obol_view_feature_set_visible(view_ctx, ref, visible);
}


int
ged_view_edit_color_set(struct ged_view_context *view_ctx,
			ged_view_edit_ref ref,
				int r,
				int g,
				int b)
{
    return ged_draw_obol_view_feature_set_color(view_ctx, ref, r, g, b);
}


int
ged_view_edit_touch(struct ged_view_context *view_ctx,
	ged_view_edit_ref ref)
{
    return ged_draw_obol_view_feature_touch(view_ctx, ref);
}


int
ged_view_edit_labels_replace(
    struct ged_view_context *view_ctx,
    ged_view_edit_ref ref,
	const struct ged_view_feature_label *labels,
	size_t label_count)
{
    return ged_draw_obol_view_feature_labels_replace(view_ctx, ref, labels,
	    label_count);
}


int
ged_view_edit_points_replace(
    struct ged_view_context *view_ctx,
    ged_view_edit_ref ref,
	enum ged_view_edit_geometry_family family,
	const point_t *points,
	const int *cmds,
	size_t point_count)
{
    return ged_draw_obol_view_feature_points_replace(view_ctx, ref, family, points,
	    cmds, point_count);
}

int
ged_view_edit_preview_replace(
    struct ged_view_context *view_ctx,
    ged_view_edit_ref ref,
	const char *source_path,
	const char *edit_intent_id,
	const char *edit_intent_role,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	uint32_t source_revision,
	uint32_t inputs_revision)
{
    return ged_draw_obol_view_feature_edit_preview_replace(view_ctx, ref, source_path,
	    edit_intent_id, edit_intent_role, points, cmds, point_count,
	    source_revision, inputs_revision);
}


int
ged_view_edit_geometry_clear(struct ged_view_context *view_ctx,
	ged_view_edit_ref ref)
{
    return ged_draw_obol_view_feature_clear_geometry(view_ctx, ref);
}


int
ged_view_feature_remove_prefix(struct ged_view_context *view_ctx, const char *prefix)
{
    if (!prefix || !prefix[0])
	return 0;
    size_t removed = 0;
    BObolViewController *view_controller =
	ged_bobol_view_controller(view_ctx);
    if (view_controller)
	removed += view_controller->features().removePrefix(prefix,
	    BOBOL_FEATURE_SCOPE_LOCAL, nullptr);
    BObolViewController *shared_controller =
	ged_bobol_shared_feature_controller(view_ctx);
    if (shared_controller)
	removed += shared_controller->features().removePrefix(prefix,
	    BOBOL_FEATURE_SCOPE_SHARED, nullptr);
    return static_cast<int>(removed);
}


int
ged_view_feature_visible(struct ged_view_context *view_ctx, const char *name)
{
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    BObolFeatureStyle style;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().style(lookup.handle, style))
	return 0;
    return (!style.hasVisible || style.visible) ? 1 : 0;
}


int
ged_view_feature_visible_set(struct ged_view_context *view_ctx, const char *name, int visible)
{
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    return lookup.controller && lookup.handle.isValid() &&
	lookup.controller->features().setVisible(lookup.handle,
	    visible ? TRUE : FALSE) ? 1 : 0;
}


int
ged_view_feature_style_get(
	struct ged_view_context *view_ctx,
	const char *name,
	struct ged_view_feature_style *style)
{
    if (!style)
	return 0;
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    BObolFeatureStyle obol_style;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().style(lookup.handle, obol_style))
	return 0;
    ged_bobol_feature_style_to_ged(style, obol_style);
    return 1;
}


int
ged_view_feature_style_apply(
	struct ged_view_context *view_ctx,
	const char *name,
	const struct ged_view_feature_style *style,
	int recursive)
{
    if (!style)
	return 0;
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    if (!lookup.controller || !lookup.handle.isValid())
	return 0;
    const BObolFeatureStyle obol_style =
	ged_bobol_feature_style_from_ged(style);
    return lookup.controller->features().applyStyle(lookup.handle, obol_style,
	recursive ? TRUE : FALSE) ? 1 : 0;
}


int
ged_view_feature_realize(struct ged_view_context *view_ctx,
				      const char *name,
				      int recursive)
{
    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    return lookup.controller && lookup.handle.isValid() &&
	lookup.controller->features().realize(lookup.handle,
	    recursive ? TRUE : FALSE) ? 1 : 0;
}


size_t
ged_view_feature_label_count(struct ged_view_context *view_ctx, const char *name)
{
    return ged_draw_obol_view_context_label_count(view_ctx, name);
}


int
ged_view_feature_label_copy(struct ged_view_context *view_ctx,
				 const char *name,
				 size_t index,
				 struct bu_vls *text,
				 point_t point,
				 unsigned char rgb[3])
{
    return ged_draw_obol_view_context_label_copy(view_ctx, name, index,
	    text, point, rgb);
}


int
ged_view_feature_axes_copy(struct ged_view_context *view_ctx,
	const char *name, size_t index, point_t center, fastf_t *half_size)
{
    if (center)
	VSETALL(center, 0.0);
    if (half_size)
	*half_size = 0.0;
    if (!name || !name[0] || !center)
	return 0;

    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    std::vector<SbVec3f> centers;
    float size = 0.0f;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().axesCenters(lookup.handle, centers,
	    &size) || index >= centers.size())
	return 0;

    VSET(center, centers[index][0], centers[index][1], centers[index][2]);
    if (half_size)
	*half_size = static_cast<fastf_t>(size);
    return 1;
}


int
ged_view_feature_points_copy(struct ged_view_context *view_ctx,
					  const char *name,
					  point_t **points,
					  size_t *point_count)
{
    return ged_draw_obol_view_context_feature_points_copy(view_ctx, name,
	    points, point_count);
}

int
ged_view_feature_indexed_face_points_update(
	struct ged_view_context *view_ctx,
	const char *name,
	const int *point_indices,
	const point_t *points,
	size_t point_count)
{
    if (!view_ctx || !name || !name[0] || !point_indices || !points ||
	!point_count)
	return 0;

    const struct ged_bobol_feature_lookup lookup =
	ged_bobol_feature_find(view_ctx, name);
    if (!lookup.controller || !lookup.handle.isValid())
	return 0;

    std::vector<int32_t> indices;
    std::vector<SbVec3f> updatedPoints;
    indices.reserve(point_count);
    updatedPoints.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	if (point_indices[i] < 0)
	    return 0;
	indices.push_back(static_cast<int32_t>(point_indices[i]));
	updatedPoints.push_back(SbVec3f(
	    static_cast<float>(points[i][X]),
	    static_cast<float>(points[i][Y]),
	    static_cast<float>(points[i][Z])));
    }

    return lookup.controller->features().updateIndexedFaceSetPoints(
	lookup.handle, indices, updatedPoints) ? 1 : 0;
}


int
ged_view_feature_line_command_at(
	struct ged_view_context *view_ctx,
	const char *name,
	size_t index,
	int *out)
{
    return ged_draw_obol_view_context_feature_line_command_at(view_ctx, name,
	    index, out);
}


static int
_ged_view_polygon_edge_color(struct bu_color *edge_color,
	const struct ged_view_polygon_record *record)
{
    if (!edge_color || !record)
	return 0;
    return bu_color_from_rgb_chars(edge_color, record->edge_color) ? 1 : 0;
}


int
ged_view_polygon_ref_is_null(ged_view_polygon_ref ref)
{
    return (ref.owner && ref.id && ref.generation) ? 0 : 1;
}


ged_view_polygon_ref
ged_view_polygon_find(struct ged_view_context *view_ctx, const char *name)
{
    return ged_draw_obol_view_context_polygon_find(view_ctx, name);
}


ged_view_polygon_ref
ged_view_polygon_find_scoped(struct ged_view_context *view_ctx,
					  const char *name,
					  int local_only)
{
    return ged_draw_obol_view_context_polygon_find_scoped(view_ctx, name,
	    local_only);
}


ged_view_polygon_ref
ged_view_polygon_select(struct ged_view_context *view_ctx,
				     const point_t model_point)
{
    return ged_draw_obol_view_context_polygon_select(view_ctx, model_point);
}


ged_view_polygon_ref
ged_view_polygon_create(struct ged_view_context *view_ctx,
			const char *name, int local,
			enum ged_view_polygon_type type,
			const point_t model_point)
{
    return ged_draw_obol_view_context_polygon_create(view_ctx, name, local,
	    type, model_point);
}


ged_view_polygon_ref
ged_view_polygon_dup(struct ged_view_context *view_ctx,
				  const char *name,
				  const char *new_name)
{
    return ged_draw_obol_view_context_polygon_dup(view_ctx, name, new_name);
}


ged_view_polygon_ref
ged_view_polygon_import_sketch(struct ged_view_context *view_ctx,
	const char *name, struct db_i *dbip, struct directory *dp, int local)
{
    return ged_draw_obol_view_context_polygon_import_sketch(name, dbip, dp,
	    view_ctx, local);
}


void
ged_view_polygon_visit_records(
	struct ged_view_context *view_ctx,
	ged_view_polygon_record_cb callback,
	void *data)
{
    ged_draw_obol_view_context_polygon_visit_records(view_ctx, callback,
	    data);
}


size_t
ged_view_polygon_snap_count(struct ged_view_context *view_ctx,
	ged_view_polygon_ref exclude)
{
    return ged_draw_obol_view_context_polygon_snap_count(view_ctx, exclude);
}


int
ged_view_polygon_clear_point_selection(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_view_context_polygon_clear_point_selection(view_ctx);
}


int
ged_view_polygon_clear_selection(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_view_context_polygon_clear_selection(view_ctx);
}


int
ged_view_polygon_snap_exclude_set(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref)
{
    return ged_draw_obol_view_context_polygon_snap_exclude_set(view_ctx, ref);
}


struct directory *
ged_view_polygon_export_sketch(struct ged_view_context *view_ctx,
	struct db_i *dbip, const char *name, ged_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_export_sketch(view_ctx, dbip, name,
	    ref);
}


struct directory *
ged_view_polygon_update_sketch(struct ged_view_context *view_ctx,
	struct db_i *dbip, ged_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_update_sketch(view_ctx, dbip, ref);
}


int
ged_view_polygon_sketch_name_set(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const char *name)
{
    return ged_draw_obol_view_polygon_sketch_name_set(view_ctx, ref, name);
}


int
ged_view_polygon_sync_sketch(struct ged *gedp,
	struct ged_view_context *view_ctx, ged_view_polygon_ref ref)
{
    struct ged_view_polygon_record record;
    if (!view_ctx || ged_view_polygon_ref_is_null(ref) ||
	!ged_view_polygon_record_get(view_ctx, ref, &record))
	return -1;

    const char *name = record.sketch_name;
    if (!name || !name[0])
	return 0;
    if (!gedp || !gedp->dbip)
	return -1;

    const bool create = db_lookup(gedp->dbip, name, LOOKUP_QUIET) ==
	RT_DIR_NULL;
    struct directory *dp = create ?
	ged_view_polygon_export_sketch(view_ctx, gedp->dbip, name, ref) :
	ged_view_polygon_update_sketch(view_ctx, gedp->dbip, ref);
    if (!dp)
	return -1;

    if (create)
	(void)ged_event_notify_object_added(gedp, name, NULL);
    else
	(void)ged_event_notify_object_modified(gedp, name, 1, NULL);
    return create ? 1 : 2;
}


int
ged_view_polygon_record_get(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, struct ged_view_polygon_record *record)
{
    if (!record)
	return 0;

    return ged_draw_obol_view_polygon_record_get(view_ctx, ref, record);
}


int
ged_view_polygon_has_data(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref)
{
    struct ged_view_polygon_record ged_record;
    return ged_view_polygon_record_get(view_ctx, ref, &ged_record) ? 1 : 0;
}


int
ged_view_polygon_update(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, enum ged_view_polygon_update op)
{
    return ged_draw_obol_view_context_polygon_update(ref, view_ctx, op);
}


int
ged_view_polygon_update_screen_pt(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, int x, int y, enum ged_view_polygon_update op)
{
    return ged_draw_obol_view_context_polygon_update_screen_pt(ref,
	    view_ctx, x, y, op);
}


int
ged_view_polygon_update_model_pt(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const point_t model_point,
	enum ged_view_polygon_update op)
{
    return ged_draw_obol_view_context_polygon_update_model_pt(ref,
	    view_ctx, model_point, op);
}


int
ged_view_polygon_move(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const point_t current_point,
	const point_t previous_point)
{
    return ged_draw_obol_view_polygon_move(view_ctx, ref, current_point,
	    previous_point);
}


int
ged_view_polygon_set_name(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const char *name)
{
    return ged_draw_obol_view_polygon_set_name(view_ctx, ref, name);
}


int
ged_view_polygon_set_visual(struct ged_view_context *view_ctx,
				 ged_view_polygon_ref ref,
				 const struct bu_color *edge_color,
				 const struct bu_color *fill_color,
				 fastf_t fill_slope_x,
				 fastf_t fill_slope_y,
				 fastf_t fill_density,
				 fastf_t vZ,
				 int fill_flag)
{
    return ged_draw_obol_view_polygon_set_visual(view_ctx, ref, edge_color,
	    fill_color, fill_slope_x, fill_slope_y, fill_density, vZ,
	    fill_flag);
}


int
ged_view_polygon_set_current(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, long contour_i, long point_i)
{
    return ged_draw_obol_view_context_polygon_set_current(view_ctx, ref,
	    contour_i, point_i);
}


int
ged_view_polygon_set_selected(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, int selected)
{
    return ged_draw_obol_view_polygon_set_selected(view_ctx, ref, selected);
}


int
ged_view_polygon_set_contour_open(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, long contour_i, int open)
{
    return ged_draw_obol_view_context_polygon_set_contour_open(view_ctx,
	    ref, contour_i, open);
}


int
ged_view_polygon_set_all_contours_open(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, int open)
{
    return ged_draw_obol_view_polygon_set_all_contours_open(view_ctx, ref,
	    open);
}


int
ged_view_polygon_close(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_close(view_ctx, ref);
}


int
ged_view_polygon_clear_selected_point(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_clear_selected_point(view_ctx, ref);
}


int
ged_view_polygon_remove(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_remove(view_ctx, ref);
}


void *
ged_view_polygon_user_data(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_user_data(view_ctx, ref);
}


int
ged_view_polygon_user_data_set(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, void *user_data)
{
    return ged_draw_obol_view_polygon_user_data_set(view_ctx, ref, user_data);
}


int
ged_view_polygon_area(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, fastf_t *area)
{
    return ged_draw_obol_view_context_polygon_area(ref, view_ctx, area);
}


int
ged_view_polygon_overlap(struct ged_view_context *view_ctx,
				      ged_view_polygon_ref ref, const char *other_name,
				      const struct bn_tol *tol,
				      int *overlap)
{
    return ged_draw_obol_view_context_polygon_overlap(ref, view_ctx,
	    other_name, tol, overlap);
}


int
ged_view_polygon_set_fill(struct ged_view_context *view_ctx,
			       ged_view_polygon_ref ref,
			       int fill_flag,
			       fastf_t fill_slope_x,
			       fastf_t fill_slope_y,
			       fastf_t fill_density)
{
    struct ged_view_polygon_record record;
    struct bu_color edge_color;
    if (!ged_view_polygon_record_get(view_ctx, ref, &record) ||
	    !_ged_view_polygon_edge_color(&edge_color, &record))
	return 0;
    return ged_view_polygon_set_visual(view_ctx, ref, &edge_color,
	    &record.fill_color, fill_slope_x, fill_slope_y, fill_density,
	    record.vZ, fill_flag);
}


int
ged_view_polygon_fill_color_get(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, struct bu_color *fill_color)
{
    struct ged_view_polygon_record record;
    if (!fill_color || !ged_view_polygon_record_get(view_ctx, ref, &record))
	return 0;
    BU_COLOR_CPY(fill_color, &record.fill_color);
    return 1;
}


int
ged_view_polygon_fill_color_set(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const struct bu_color *fill_color)
{
    struct ged_view_polygon_record record;
    struct bu_color edge_color;
    if (!fill_color || !ged_view_polygon_record_get(view_ctx, ref, &record) ||
	    !_ged_view_polygon_edge_color(&edge_color, &record))
	return 0;
    return ged_view_polygon_set_visual(view_ctx, ref, &edge_color, fill_color,
	    record.fill_dir[0], record.fill_dir[1], record.fill_delta,
	    record.vZ, record.fill_flag);
}


int
ged_view_polygon_csg_name(struct ged_view_context *view_ctx,
				  ged_view_polygon_ref target, const char *other_name,
				  enum bg_polygon_boolean_op op)
{
    if (!view_ctx || !other_name)
	return 0;

    ged_view_polygon_ref other_ref =
	ged_view_polygon_find(view_ctx, other_name);
    if (ged_view_polygon_ref_is_null(other_ref))
	return 0;

    return ged_view_polygon_csg(view_ctx, target, other_ref, op);
}


int
ged_view_polygon_csg(struct ged_view_context *view_ctx,
			  ged_view_polygon_ref target, ged_view_polygon_ref stencil,
			  enum bg_polygon_boolean_op op)
{
    return ged_draw_obol_view_polygon_csg(view_ctx, target, stencil, op);
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
