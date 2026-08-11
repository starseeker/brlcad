/*                  G E D _ S C E N E _ A P I . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged_scene_api.cpp
 *
 * Typed semantic scene operations.  During the drawing-state consolidation
 * these requests feed the existing internal transaction reducer; clients do
 * not need to construct or interpret that migration representation.
 */

#include "common.h"

#include <climits>
#include <string>
#include <vector>

#include "bu/vls.h"
#include "ged/scene.h"
#include "ged/view.h"

#include "./ged_draw_private.h"
#include "./ged_private.h"
#include "./ged_scene_backend_private.h"

struct ged_scene_result {
    enum ged_scene_status status;
    int changed;
    size_t group_count;
    size_t shape_count;
    size_t conflict_count;
    uint64_t revision_before;
    uint64_t revision_after;
    std::vector<std::string> paths;
    std::string diagnostic;
};

static ged_scene_occurrence_ref ged_scene_occurrence_from_shape_ref(
    struct ged *gedp, ged_draw_shape_ref ref);
static ged_draw_shape_ref ged_scene_occurrence_to_shape_ref(
    struct ged *gedp, ged_scene_occurrence_ref occurrence);


static void
ged_scene_result_reset(struct ged_scene_result *result)
{
    if (!result)
	return;

    result->status = GED_SCENE_OK;
    result->changed = 0;
    result->group_count = 0;
    result->shape_count = 0;
    result->conflict_count = 0;
    result->revision_before = 0;
    result->revision_after = 0;
    result->paths.clear();
    result->diagnostic.clear();
}


static void
ged_scene_result_path_note(struct ged_scene_result *result, const char *path)
{
    if (!result || !path || !path[0])
	return;

    for (const std::string &existing : result->paths) {
	if (existing == path)
	    return;
    }
    result->paths.emplace_back(path);
}


static void
ged_scene_result_from_draw_result(
    struct ged_scene_result *result,
    const struct ged_draw_transaction_result *draw_result,
    int transaction_status)
{
    if (!result || !draw_result)
	return;

    result->status = transaction_status < 0 ? GED_SCENE_ERROR : GED_SCENE_OK;
    result->changed = transaction_status > 0 ||
	draw_result->scene_revision_before != draw_result->scene_revision_after;
    result->group_count = draw_result->affected_groups > 0 ?
	static_cast<size_t>(draw_result->affected_groups) : 0;
    result->shape_count = draw_result->affected_shapes > 0 ?
	static_cast<size_t>(draw_result->affected_shapes) : 0;
    result->revision_before = draw_result->scene_revision_before;
    result->revision_after = draw_result->scene_revision_after;
    if (bu_vls_strlen(&draw_result->errors))
	result->diagnostic = bu_vls_cstr(&draw_result->errors);
}


static int
ged_scene_draw_mode_internal(const struct ged *gedp,
			     enum ged_scene_draw_mode mode)
{
    if (mode == GED_SCENE_DRAW_DEFAULT)
	return ged_draw_default_mode(gedp);
    if (mode < GED_SCENE_DRAW_WIRE || mode > GED_SCENE_DRAW_EVALUATED_POINTS)
	return -1;
    return static_cast<int>(mode);
}


extern "C" void
ged_scene_style_init(struct ged_scene_style *style)
{
    if (!style)
	return;

    style->draw_mode = GED_SCENE_DRAW_DEFAULT;
    style->opacity = 1.0;
    style->color[0] = 255;
    style->color[1] = 0;
    style->color[2] = 0;
    style->color_override = 0;
    style->line_width = 1;
    style->mixed_modes = 0;
    style->flags = 0;
}


extern "C" void
ged_scene_realization_policy_init(
    struct ged_scene_realization_policy *policy)
{
    if (!policy)
	return;

    policy->mode = GED_SCENE_REALIZE_AUTO;
    policy->strict = 0;
}


extern "C" void
ged_scene_draw_request_init(struct ged_scene_draw_request *request)
{
    if (!request)
	return;

    request->view = NULL;
    request->paths = NULL;
    request->path_count = 0;
    ged_scene_style_init(&request->style);
    ged_scene_realization_policy_init(&request->realization);
    request->autoview = 0;
}


extern "C" void
ged_scene_erase_request_init(struct ged_scene_erase_request *request)
{
    if (!request)
	return;

    request->view = NULL;
    request->path = NULL;
    request->match = GED_SCENE_ERASE_EXACT;
    request->draw_mode = GED_SCENE_DRAW_DEFAULT;
}


extern "C" void
ged_scene_path_request_init(struct ged_scene_path_request *request)
{
    if (!request)
	return;

    request->view = NULL;
    request->path = NULL;
    request->draw_mode = GED_SCENE_DRAW_DEFAULT;
    request->match = GED_SCENE_PATH_MATCH_EXACT;
}


extern "C" void
ged_scene_redraw_request_init(struct ged_scene_redraw_request *request)
{
    if (!request)
	return;

    request->view = NULL;
    request->paths = NULL;
    request->path_count = 0;
}


extern "C" void
ged_scene_clear_request_init(struct ged_scene_clear_request *request)
{
    if (!request)
	return;

    request->view = NULL;
    request->scope = GED_SCENE_CLEAR_ALL;
}


extern "C" void
ged_scene_edit_request_init(struct ged_scene_edit_request *request)
{
    if (!request)
	return;

    request->view = NULL;
    request->path = NULL;
    request->draw_mode = GED_SCENE_DRAW_DEFAULT;
    request->occurrences = GED_SCENE_EDIT_EXACT_OCCURRENCE;
    request->draw_if_absent = 1;
    request->purpose = NULL;
}


extern "C" struct ged_scene_result *
ged_scene_result_create(void)
{
    struct ged_scene_result *result = new ged_scene_result;
    ged_scene_result_reset(result);
    return result;
}


extern "C" void
ged_scene_result_clear(struct ged_scene_result *result)
{
    ged_scene_result_reset(result);
}


extern "C" void
ged_scene_result_destroy(struct ged_scene_result *result)
{
    delete result;
}


extern "C" enum ged_scene_status
ged_scene_result_status(const struct ged_scene_result *result)
{
    return result ? result->status : GED_SCENE_INVALID;
}


extern "C" int
ged_scene_result_changed(const struct ged_scene_result *result)
{
    return result && result->changed ? 1 : 0;
}


extern "C" size_t
ged_scene_result_path_count(const struct ged_scene_result *result)
{
    return result ? result->paths.size() : 0;
}


extern "C" const char *
ged_scene_result_path_at(const struct ged_scene_result *result, size_t index)
{
    return result && index < result->paths.size() ?
	result->paths[index].c_str() : NULL;
}


extern "C" size_t
ged_scene_result_group_count(const struct ged_scene_result *result)
{
    return result ? result->group_count : 0;
}


extern "C" size_t
ged_scene_result_shape_count(const struct ged_scene_result *result)
{
    return result ? result->shape_count : 0;
}


extern "C" size_t
ged_scene_result_conflict_count(const struct ged_scene_result *result)
{
    return result ? result->conflict_count : 0;
}


extern "C" uint64_t
ged_scene_result_revision_before(const struct ged_scene_result *result)
{
    return result ? result->revision_before : 0;
}


extern "C" uint64_t
ged_scene_result_revision_after(const struct ged_scene_result *result)
{
    return result ? result->revision_after : 0;
}


extern "C" const char *
ged_scene_result_diagnostic(const struct ged_scene_result *result)
{
    return result ? result->diagnostic.c_str() : "";
}


extern "C" enum ged_scene_status
ged_scene_draw(struct ged *gedp,
	       const struct ged_scene_draw_request *request,
	       struct ged_scene_result *result)
{
    if (result)
	ged_scene_result_reset(result);
    if (!gedp || !request || !request->paths || !request->path_count ||
	request->path_count > static_cast<size_t>(INT_MAX) ||
	request->style.opacity < 0.0 || request->style.opacity > 1.0 ||
	request->style.line_width < 1) {
	if (result)
	    result->status = GED_SCENE_INVALID;
	return GED_SCENE_INVALID;
    }

    const int draw_mode = ged_scene_draw_mode_internal(gedp,
	request->style.draw_mode);
    if (draw_mode < 0) {
	if (result)
	    result->status = GED_SCENE_INVALID;
	return GED_SCENE_INVALID;
    }

    struct ged_view_context *view = request->view ? request->view :
	ged_view_active_ctx(gedp);
    const int evaluated = draw_mode == GED_DRAW_MODE_EVAL_WIRE ||
	draw_mode == GED_DRAW_MODE_EVAL_POINTS;
    int progressive = 0;
    switch (request->realization.mode) {
	case GED_SCENE_REALIZE_AUTO:
	    progressive = !evaluated && view &&
		ged_draw_obol_progressive_available(gedp, view) &&
		ged_draw_obol_view_lod_enabled(gedp, view, draw_mode);
	    break;
	case GED_SCENE_REALIZE_EAGER:
	    progressive = 0;
	    break;
	case GED_SCENE_REALIZE_PROGRESSIVE:
	    if (evaluated || !view ||
		!ged_draw_obol_progressive_available(gedp, view)) {
		if (result)
		    result->status = GED_SCENE_ERROR;
		return GED_SCENE_ERROR;
	    }
	    progressive = 1;
	    break;
	default:
	    if (result)
		result->status = GED_SCENE_INVALID;
	    return GED_SCENE_INVALID;
    }

    struct ged_draw_appearance_settings appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    appearance.draw_mode = draw_mode;
    appearance.mixed_modes = request->style.mixed_modes ? 1 : 0;
    appearance.transparency = request->style.opacity;
    appearance.color_override = request->style.color_override ? 1 : 0;
    appearance.color[0] = request->style.color[0];
    appearance.color[1] = request->style.color[1];
    appearance.color[2] = request->style.color[2];
    appearance.s_line_width = request->style.line_width;
    appearance.draw_solid_lines_only =
	(request->style.flags & GED_SCENE_STYLE_SOLID_LINES_ONLY) ? 1 : 0;
    appearance.draw_non_subtract_only =
	(request->style.flags & GED_SCENE_STYLE_NON_SUBTRACT_ONLY) ? 1 : 0;
    appearance.strict_fallback = request->realization.strict ? 1 : 0;
    appearance.defer_leaf_expansion = progressive;

    struct ged_draw_transaction transaction =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, NULL);
    transaction.view = view;
    transaction.paths = const_cast<const char **>(request->paths);
    transaction.path_count = static_cast<int>(request->path_count);
    transaction.appearance = &appearance;
    transaction.autoview = request->autoview ? 1 : 0;

    struct ged_draw_transaction_result draw_result;
    ged_draw_transaction_result_init(&draw_result);
    const int ret = ged_draw_apply_transaction(gedp, &transaction,
	&draw_result);
    if (result) {
	ged_scene_result_from_draw_result(result, &draw_result, ret);
	for (size_t i = 0; i < request->path_count; i++)
	    ged_scene_result_path_note(result, request->paths[i]);
    }
    ged_draw_transaction_result_free(&draw_result);
    return ret < 0 ? GED_SCENE_ERROR : GED_SCENE_OK;
}


extern "C" enum ged_scene_status
ged_scene_erase(struct ged *gedp,
		const struct ged_scene_erase_request *request,
		struct ged_scene_result *result)
{
    if (result)
	ged_scene_result_reset(result);
    if (!gedp || !request || !request->path || !request->path[0] ||
	(request->match != GED_SCENE_ERASE_EXACT &&
	 request->match != GED_SCENE_ERASE_SUBTREE)) {
	if (result)
	    result->status = GED_SCENE_INVALID;
	return GED_SCENE_INVALID;
    }

    int draw_mode = -1;
    if (request->draw_mode != GED_SCENE_DRAW_DEFAULT) {
	draw_mode = ged_scene_draw_mode_internal(gedp, request->draw_mode);
	if (draw_mode < 0) {
	    if (result)
		result->status = GED_SCENE_INVALID;
	    return GED_SCENE_INVALID;
	}
    }

    const ged_draw_transaction_kind kind =
	request->match == GED_SCENE_ERASE_SUBTREE ?
	GED_DRAW_TXN_ERASE_PREFIX : GED_DRAW_TXN_ERASE;
    struct ged_draw_transaction transaction =
	ged_draw_transaction_make(kind, request->path);
    transaction.view = request->view;
    transaction.mode = draw_mode;

    struct ged_draw_transaction_result draw_result;
    ged_draw_transaction_result_init(&draw_result);
    const int ret = ged_draw_apply_transaction(gedp, &transaction,
	&draw_result);
    if (result) {
	ged_scene_result_from_draw_result(result, &draw_result, ret);
	ged_scene_result_path_note(result, request->path);
    }
    ged_draw_transaction_result_free(&draw_result);
    return ret < 0 ? GED_SCENE_ERROR : GED_SCENE_OK;
}


static enum ged_scene_status
ged_scene_apply_typed_transaction(
    struct ged *gedp,
    struct ged_draw_transaction *transaction,
    const char *const *paths,
    size_t path_count,
    struct ged_scene_result *result)
{
    if (result)
	ged_scene_result_reset(result);
    if (!gedp || !transaction) {
	if (result)
	    result->status = GED_SCENE_INVALID;
	return GED_SCENE_INVALID;
    }

    struct ged_draw_transaction_result draw_result;
    ged_draw_transaction_result_init(&draw_result);
    const int ret = ged_draw_apply_transaction(gedp, transaction,
	&draw_result);
    if (result) {
	ged_scene_result_from_draw_result(result, &draw_result, ret);
	for (size_t i = 0; paths && i < path_count; i++)
	    ged_scene_result_path_note(result, paths[i]);
    }
    ged_draw_transaction_result_free(&draw_result);
    return ret < 0 ? GED_SCENE_ERROR : GED_SCENE_OK;
}


extern "C" enum ged_scene_status
ged_scene_redraw(struct ged *gedp,
		 const struct ged_scene_redraw_request *request,
		 struct ged_scene_result *result)
{
    if (!request || (request->path_count && !request->paths) ||
	request->path_count > static_cast<size_t>(INT_MAX)) {
	if (result) {
	    ged_scene_result_reset(result);
	    result->status = GED_SCENE_INVALID;
	}
	return GED_SCENE_INVALID;
    }

    struct ged_draw_transaction transaction =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    transaction.view = request->view;
    transaction.paths = const_cast<const char **>(request->paths);
    transaction.path_count = static_cast<int>(request->path_count);
    return ged_scene_apply_typed_transaction(gedp, &transaction,
	request->paths, request->path_count, result);
}


static enum ged_scene_status
ged_scene_path_value_set(struct ged *gedp,
			 const struct ged_scene_path_request *request,
			 enum ged_draw_transaction_kind kind,
			 double value,
			 struct ged_scene_result *result)
{
    if (!request || !request->path || !request->path[0] ||
	(request->match != GED_SCENE_PATH_MATCH_EXACT &&
	 request->match != GED_SCENE_PATH_MATCH_SUBTREE &&
	 request->match != GED_SCENE_PATH_MATCH_OBJECT)) {
	if (result) {
	    ged_scene_result_reset(result);
	    result->status = GED_SCENE_INVALID;
	}
	return GED_SCENE_INVALID;
    }

    int draw_mode = -1;
    if (request->draw_mode != GED_SCENE_DRAW_DEFAULT) {
	draw_mode = ged_scene_draw_mode_internal(gedp, request->draw_mode);
	if (draw_mode < 0) {
	    if (result) {
		ged_scene_result_reset(result);
		result->status = GED_SCENE_INVALID;
	    }
	    return GED_SCENE_INVALID;
	}
    }

    struct ged_draw_transaction transaction =
	ged_draw_transaction_make_value(kind, request->path, value);
    transaction.view = request->view;
    transaction.mode = draw_mode;
    transaction.match = request->match;
    const char *paths[] = {request->path};
    return ged_scene_apply_typed_transaction(gedp, &transaction, paths, 1,
	result);
}


extern "C" enum ged_scene_status
ged_scene_visibility_set(struct ged *gedp,
			 const struct ged_scene_path_request *request,
			 int visible,
			 struct ged_scene_result *result)
{
    return ged_scene_path_value_set(gedp, request, GED_DRAW_TXN_VISIBILITY,
	visible ? 1.0 : 0.0, result);
}


extern "C" enum ged_scene_status
ged_scene_opacity_set(struct ged *gedp,
		      const struct ged_scene_path_request *request,
		      double opacity,
		      struct ged_scene_result *result)
{
    if (opacity < 0.0 || opacity > 1.0) {
	if (result) {
	    ged_scene_result_reset(result);
	    result->status = GED_SCENE_INVALID;
	}
	return GED_SCENE_INVALID;
    }
    return ged_scene_path_value_set(gedp, request, GED_DRAW_TXN_TRANSPARENCY,
	1.0 - opacity, result);
}


extern "C" enum ged_scene_status
ged_scene_highlight_set(struct ged *gedp,
			const struct ged_scene_path_request *request,
			int highlighted,
			struct ged_scene_result *result)
{
    return ged_scene_path_value_set(gedp, request, GED_DRAW_TXN_HIGHLIGHT,
	highlighted ? 1.0 : 0.0, result);
}


extern "C" enum ged_scene_status
ged_scene_highlights_clear(struct ged *gedp, struct ged_scene_result *result)
{
    struct ged_draw_transaction transaction = ged_draw_transaction_make(
	GED_DRAW_TXN_HIGHLIGHTS_CLEAR, NULL);
    return ged_scene_apply_typed_transaction(gedp, &transaction, NULL, 0,
	result);
}


extern "C" enum ged_scene_status
ged_scene_occurrence_highlight_set(struct ged *gedp,
	ged_scene_occurrence_ref occurrence, int highlighted,
	struct ged_scene_result *result)
{
    ged_draw_shape_ref shape_ref = ged_scene_occurrence_to_shape_ref(gedp,
	occurrence);
    if (!gedp || (ged_scene_occurrence_ref_is_null(occurrence) && highlighted) ||
	(!ged_scene_occurrence_ref_is_null(occurrence) &&
	 ged_draw_shape_ref_is_null(shape_ref))) {
	if (result) {
	    ged_scene_result_reset(result);
	    result->status = GED_SCENE_INVALID;
	}
	return GED_SCENE_INVALID;
    }
    struct ged_draw_transaction transaction = ged_draw_transaction_make_value(
	GED_DRAW_TXN_HIGHLIGHT_OCCURRENCE, NULL, highlighted ? 1.0 : 0.0);
    transaction.shape_ref = shape_ref;
    return ged_scene_apply_typed_transaction(gedp, &transaction, NULL, 0,
	result);
}


extern "C" enum ged_scene_status
ged_scene_clear(struct ged *gedp,
		const struct ged_scene_clear_request *request,
		struct ged_scene_result *result)
{
    if (!request || (request->scope != GED_SCENE_CLEAR_ALL &&
	request->scope != GED_SCENE_CLEAR_VIEW)) {
	if (result) {
	    ged_scene_result_reset(result);
	    result->status = GED_SCENE_INVALID;
	}
	return GED_SCENE_INVALID;
    }

    struct ged_draw_transaction transaction = ged_draw_transaction_make(
	request->scope == GED_SCENE_CLEAR_ALL ? GED_DRAW_TXN_CLEAR :
	GED_DRAW_TXN_CLEAR_SCOPE, NULL);
    transaction.view = request->view;
    return ged_scene_apply_typed_transaction(gedp, &transaction, NULL, 0,
	result);
}


extern "C" enum ged_scene_status
ged_scene_default_draw_mode_set(struct ged *gedp,
				enum ged_scene_draw_mode mode,
				struct ged_scene_result *result)
{
    const int draw_mode = ged_scene_draw_mode_internal(gedp, mode);
    if (mode == GED_SCENE_DRAW_DEFAULT || draw_mode < 0) {
	if (result) {
	    ged_scene_result_reset(result);
	    result->status = GED_SCENE_INVALID;
	}
	return GED_SCENE_INVALID;
    }
    struct ged_draw_transaction transaction = ged_draw_transaction_make_value(
	GED_DRAW_TXN_DEFAULT_DRAW_MODE, NULL, static_cast<double>(draw_mode));
    return ged_scene_apply_typed_transaction(gedp, &transaction, NULL, 0,
	result);
}


extern "C" enum ged_scene_status
ged_scene_materials_changed(struct ged *gedp,
			    struct ged_scene_result *result)
{
    struct ged_draw_transaction transaction = ged_draw_transaction_make(
	GED_DRAW_TXN_MATERIAL_CHANGED, NULL);
    return ged_scene_apply_typed_transaction(gedp, &transaction, NULL, 0,
	result);
}


extern "C" uint64_t
ged_scene_material_revision(const struct ged *gedp)
{
    return ged_draw_material_revision(gedp);
}


extern "C" uint64_t
ged_scene_revision(const struct ged *gedp)
{
    return gedp ? ged_draw_scene_revision(const_cast<struct ged *>(gedp)) : 0;
}


extern "C" int
ged_scene_backend_apply_private(
    struct ged *gedp,
    const struct ged_draw_transaction *transaction,
    const struct ged_draw_transaction_result *result)
{
    if (!gedp || !transaction)
	return 0;
    Ged_Internal *state = gedp->i && gedp->i->i ? gedp->i->i : nullptr;
    if (state && state->scene_backend_ops) {
	return state->scene_backend_ops->apply ?
	    state->scene_backend_ops->apply(gedp, transaction, result,
		state->scene_backend_data) : 0;
    }
    return ged_draw_obol_backend_apply_transaction(gedp, transaction, result);
}


extern "C" int
ged_scene_backend_snapshot_private(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return 0;
    Ged_Internal *state = gedp->i->i;
    if (state->scene_backend_ops) {
	return state->scene_backend_ops->snapshot ?
	    state->scene_backend_ops->snapshot(gedp,
		state->scene_backend_data) : 0;
    }
    return ged_draw_obol_scene_sync_full_scene(gedp, nullptr, 0, nullptr);
}


extern "C" int
ged_scene_backend_selection_private(
    struct ged *gedp,
    const char *const *added_paths, size_t added_count,
    const char *const *removed_paths, size_t removed_count,
    const char *const *selected_paths, size_t selected_count)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return 0;
    Ged_Internal *state = gedp->i->i;
    if (state->scene_backend_ops) {
	return state->scene_backend_ops->selection ?
	    state->scene_backend_ops->selection(gedp, added_paths,
		added_count, removed_paths, removed_count, selected_paths,
		selected_count, state->scene_backend_data) : 0;
    }
    return ged_draw_obol_database_sources_apply_selection_delta(gedp,
	added_paths, added_count, removed_paths, removed_count,
	selected_paths, selected_count);
}


extern "C" void
ged_scene_backend_set_private(
    struct ged *gedp,
    const struct ged_scene_backend_ops *operations,
    void *client_data)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;
    gedp->i->i->scene_backend_ops = operations;
    gedp->i->i->scene_backend_data = operations ? client_data : nullptr;
    (void)ged_scene_backend_snapshot_private(gedp);
}


extern "C" int
ged_scene_available(const struct ged *gedp)
{
    return gedp ? ged_draw_scene_available(const_cast<struct ged *>(gedp)) : 0;
}


extern "C" enum ged_scene_draw_mode
ged_scene_default_draw_mode(const struct ged *gedp)
{
    const int mode = ged_draw_default_mode(gedp);
    if (mode < GED_DRAW_MODE_WIRE || mode > GED_DRAW_MODE_EVAL_POINTS)
	return GED_SCENE_DRAW_WIRE;
    return static_cast<enum ged_scene_draw_mode>(mode);
}


extern "C" int
ged_scene_has_paths(struct ged *gedp, struct ged_view_context *view,
		    enum ged_scene_draw_mode mode)
{
    const int internal_mode = mode == GED_SCENE_DRAW_DEFAULT ? -1 :
	ged_scene_draw_mode_internal(gedp, mode);
    if (!gedp || internal_mode < -1)
	return 0;
    return ged_draw_has_paths(gedp, view, internal_mode) ? 1 : 0;
}


extern "C" enum ged_scene_path_state
ged_scene_path_state_get(struct ged *gedp, struct ged_view_context *view,
			 const char *path, enum ged_scene_draw_mode mode)
{
    const int internal_mode = mode == GED_SCENE_DRAW_DEFAULT ? -1 :
	ged_scene_draw_mode_internal(gedp, mode);
    if (!gedp || !path || internal_mode < -1)
	return GED_SCENE_PATH_NOT_DRAWN;
    const int state = ged_draw_path_state(gedp, view, path, internal_mode);
    if (state == 1)
	return GED_SCENE_PATH_DRAWN;
    if (state == 2)
	return GED_SCENE_PATH_PARTIALLY_DRAWN;
    return GED_SCENE_PATH_NOT_DRAWN;
}


extern "C" size_t
ged_scene_paths_append(struct ged *gedp, struct ged_view_context *view,
		       enum ged_scene_draw_mode mode,
		       enum ged_scene_path_listing listing,
		       struct bu_vls *result)
{
    const int internal_mode = mode == GED_SCENE_DRAW_DEFAULT ? -1 :
	ged_scene_draw_mode_internal(gedp, mode);
    if (!gedp || !result || internal_mode < -1 ||
	(listing != GED_SCENE_PATHS_DRAW_INTENTS &&
	 listing != GED_SCENE_PATHS_REALIZED_OCCURRENCES))
	return 0;
    return ged_draw_list_paths(gedp, view, internal_mode,
	listing == GED_SCENE_PATHS_REALIZED_OCCURRENCES, result);
}


extern "C" size_t
ged_scene_occurrence_count(struct ged *gedp)
{
    const int count = gedp ? ged_draw_shape_count(gedp) : 0;
    if (count > 0)
	return static_cast<size_t>(count);
    /* Compact draw roots deliberately do not mirror every streamed leaf into
     * libged records.  Until the occurrence-query backend contract grows a
     * cheap compact count, report retained realized sources rather than
     * incorrectly claiming an attached compact scene is empty. */
    size_t source_count = 0;
    return gedp && ged_draw_obol_database_source_count(gedp, 0,
	&source_count) ? source_count : 0;
}


extern "C" int
ged_scene_has_occurrences(struct ged *gedp)
{
    return ged_scene_occurrence_count(gedp) ? 1 : 0;
}


static ged_scene_occurrence_ref
ged_scene_occurrence_from_shape_ref(struct ged *gedp, ged_draw_shape_ref ref)
{
    ged_scene_occurrence_ref occurrence = GED_SCENE_OCCURRENCE_REF_NULL;
    if (!gedp || ged_draw_shape_ref_is_null(ref))
	return occurrence;
    occurrence.owner = reinterpret_cast<uintptr_t>(gedp);
    occurrence.id = static_cast<uint64_t>(ref.token);
    occurrence.generation = 1;
    return occurrence;
}


static ged_draw_shape_ref
ged_scene_occurrence_to_shape_ref(struct ged *gedp,
	ged_scene_occurrence_ref occurrence)
{
    ged_draw_shape_ref ref = GED_DRAW_SHAPE_REF_NULL;
    if (!gedp || !occurrence.owner || !occurrence.id ||
	occurrence.owner != reinterpret_cast<uintptr_t>(gedp) ||
	occurrence.generation != 1)
	return ref;
    ref.token = static_cast<uintptr_t>(occurrence.id);
    ref.scene_revision = ged_draw_scene_revision(gedp);
    struct ged_draw_shape_record record;
    if (!ged_draw_shape_record_get(gedp, ref, &record))
	return GED_DRAW_SHAPE_REF_NULL;
    return ref;
}


extern "C" int
ged_scene_occurrence_ref_is_null(ged_scene_occurrence_ref occurrence)
{
    return !occurrence.owner && !occurrence.id && !occurrence.generation;
}


extern "C" int
ged_scene_occurrence_ref_equal(ged_scene_occurrence_ref a,
	ged_scene_occurrence_ref b)
{
    return a.owner == b.owner && a.id == b.id &&
	a.generation == b.generation;
}


static void
ged_scene_occurrence_info_fill(struct ged *gedp,
	const struct ged_draw_shape_record *record,
	struct ged_scene_occurrence_info *out)
{
    if (!record || !out)
	return;
    *out = ged_scene_occurrence_info();
    out->ref = ged_scene_occurrence_from_shape_ref(gedp, record->ref);
    out->fullpath = record->fullpath;
    out->path = record->display_name;
    out->leaf_name = record->leaf_name;
    out->path_hash = record->path_hash;
    out->view = ged_draw_shape_ref_view_context(gedp, record->ref);
    out->draw_mode = static_cast<enum ged_scene_draw_mode>(record->draw_mode);
    out->opacity = 1.0 - record->transparency;
    out->visible = record->visible;
    out->highlighted = record->highlighted;
    out->selected = record->selected;
    out->evaluated_region = record->evaluated_region;
    out->line_width = record->line_width;
    VMOVE(out->center, record->center);
}


struct ged_scene_occurrence_visit_context {
    struct ged *gedp;
    ged_scene_occurrence_func_t callback;
    void *client_data;
    size_t count;
};


static int
ged_scene_occurrence_visit_cb(const struct ged_draw_shape_record *record,
	void *client_data)
{
    struct ged_scene_occurrence_visit_context *context =
	static_cast<struct ged_scene_occurrence_visit_context *>(client_data);
    if (!context || !context->callback || !record)
	return 0;
    struct ged_scene_occurrence_info occurrence;
    ged_scene_occurrence_info_fill(context->gedp, record, &occurrence);
    context->count++;
    return context->callback(&occurrence, context->client_data);
}


extern "C" size_t
ged_scene_occurrences_visit(struct ged *gedp,
	ged_scene_occurrence_func_t callback, void *client_data)
{
    if (!gedp || !callback)
	return 0;
    struct ged_scene_occurrence_visit_context context = {
	gedp, callback, client_data, 0
    };
    ged_draw_foreach_shape_record(gedp, ged_scene_occurrence_visit_cb,
	&context);
    return context.count;
}


struct ged_scene_candidate_visit_context {
    ged_scene_occurrence_candidate_func_t callback;
    void *client_data;
    size_t count;
};


static int
ged_scene_candidate_visit_cb(const struct ged_draw_shape_candidate *candidate,
	void *client_data)
{
    struct ged_scene_candidate_visit_context *context =
	static_cast<struct ged_scene_candidate_visit_context *>(client_data);
    if (!context || !context->callback || !candidate)
	return 0;
    struct ged_scene_occurrence_candidate semantic;
    semantic.path = candidate->path;
    semantic.instance_key = candidate->instance_key;
    semantic.draw_mode =
	static_cast<enum ged_scene_draw_mode>(candidate->draw_mode);
    context->count++;
    return context->callback(&semantic, context->client_data);
}


extern "C" size_t
ged_scene_occurrence_candidates_visit(
    struct ged *gedp, ged_scene_occurrence_candidate_func_t callback,
    void *client_data)
{
    if (!gedp || !callback)
	return 0;
    struct ged_scene_candidate_visit_context context = {
	callback, client_data, 0
    };
    ged_draw_foreach_visible_shape_candidate(gedp,
	ged_scene_candidate_visit_cb, &context);
    return context.count;
}


extern "C" ged_scene_occurrence_ref
ged_scene_occurrence_candidate_resolve(
    struct ged *gedp, const struct ged_scene_occurrence_candidate *candidate)
{
    if (!gedp || !candidate)
	return GED_SCENE_OCCURRENCE_REF_NULL;
    struct ged_draw_shape_candidate internal;
    internal.path = candidate->path;
    internal.instance_key = candidate->instance_key;
    internal.draw_mode = static_cast<int>(candidate->draw_mode);
    return ged_scene_occurrence_from_shape_ref(gedp,
	ged_draw_shape_ref_for_candidate(gedp, &internal));
}


extern "C" int
ged_scene_occurrence_get(struct ged *gedp,
	ged_scene_occurrence_ref occurrence,
	struct ged_scene_occurrence_info *out)
{
    if (!out)
	return 0;
    ged_draw_shape_ref ref = ged_scene_occurrence_to_shape_ref(gedp,
	occurrence);
    struct ged_draw_shape_record record;
    if (ged_draw_shape_ref_is_null(ref) ||
	!ged_draw_shape_record_get(gedp, ref, &record))
	return 0;
    ged_scene_occurrence_info_fill(gedp, &record, out);
    return 1;
}


extern "C" ged_scene_occurrence_ref
ged_scene_occurrence_first(struct ged *gedp)
{
    return ged_scene_occurrence_from_shape_ref(gedp,
	ged_draw_first_shape_ref(gedp));
}


extern "C" ged_scene_occurrence_ref
ged_scene_occurrence_advance(struct ged *gedp,
	ged_scene_occurrence_ref occurrence, int delta)
{
    ged_draw_shape_ref ref = ged_scene_occurrence_to_shape_ref(gedp,
	occurrence);
    return ged_scene_occurrence_from_shape_ref(gedp,
	ged_draw_advance_shape_ref(gedp, ref, delta));
}


extern "C" int
ged_scene_bounds(struct ged *gedp, vect_t *min, vect_t *max,
		 enum ged_scene_bounds_scope scope)
{
    if (!gedp || !min || !max ||
	(scope != GED_SCENE_BOUNDS_DATABASE && scope != GED_SCENE_BOUNDS_ALL))
	return 1;
    return ged_draw_bounds(gedp, min, max,
	scope == GED_SCENE_BOUNDS_ALL ? 1 : 0);
}


extern "C" uint64_t
ged_scene_highlight_revision(const struct ged *gedp)
{
    return gedp ? ged_draw_highlight_revision(gedp) : 0;
}


static std::vector<std::string>
ged_scene_lines(const struct bu_vls *paths)
{
    std::vector<std::string> lines;
    if (!paths || !bu_vls_strlen(paths))
	return lines;
    const char *begin = bu_vls_cstr(paths);
    const char *cursor = begin;
    while (true) {
	if (*cursor != '\n' && *cursor != '\r' && *cursor != '\0') {
	    cursor++;
	    continue;
	}
	if (cursor > begin) {
	    std::string path(begin, static_cast<size_t>(cursor - begin));
	    lines.push_back(path);
	}
	if (*cursor == '\0')
	    break;
	begin = ++cursor;
    }
    return lines;
}


static void
ged_scene_result_note_lines(struct ged_scene_result *result,
			    const std::vector<std::string> &paths)
{
    if (!result)
	return;
    for (const std::string &path : paths)
	ged_scene_result_path_note(result, path.c_str());
}


extern "C" int
ged_scene_edit_scope_ref_is_null(ged_scene_edit_scope_ref scope)
{
    return (!scope.owner || !scope.id || !scope.generation) ? 1 : 0;
}


extern "C" enum ged_scene_status
ged_scene_edit_acquire(struct ged *gedp,
		       const struct ged_scene_edit_request *request,
		       ged_scene_edit_scope_ref *scope,
		       struct ged_scene_result *result)
{
    if (result)
	ged_scene_result_reset(result);
    if (scope)
	*scope = GED_SCENE_EDIT_SCOPE_REF_NULL;
    if (!gedp || !request || !scope || !request->path ||
	!request->path[0] ||
	(request->occurrences != GED_SCENE_EDIT_EXACT_OCCURRENCE &&
	 request->occurrences != GED_SCENE_EDIT_ALL_DRAWN_OCCURRENCES)) {
	if (result)
	    result->status = GED_SCENE_INVALID;
	return GED_SCENE_INVALID;
    }

    if (request->draw_mode != GED_SCENE_DRAW_DEFAULT &&
	ged_scene_draw_mode_internal(gedp, request->draw_mode) < 0) {
	    if (result)
		result->status = GED_SCENE_INVALID;
	    return GED_SCENE_INVALID;
    }

    struct ged_scene_edit_internal_result internal_result;
    ged_scene_edit_internal_result_init(&internal_result);
    const uint64_t revision_before = ged_draw_scene_revision(gedp);
    const int ret = ged_scene_edit_acquire_internal(gedp, request, scope,
	&internal_result);
    const uint64_t revision_after = ret > 0 ?
	ged_scene_revision_advance(gedp) : ged_draw_scene_revision(gedp);
    const std::vector<std::string> changed_paths =
	ged_scene_lines(&internal_result.paths);
    if (result) {
	result->status = ret < 0 ? GED_SCENE_ERROR :
	    (ret > 0 ? GED_SCENE_OK : GED_SCENE_NOT_FOUND);
	result->changed = ret > 0;
	result->group_count = internal_result.replaced_root_count;
	result->revision_before = revision_before;
	result->revision_after = revision_after;
	result->conflict_count = internal_result.conflict_count;
	result->diagnostic = bu_vls_cstr(&internal_result.errors);
	ged_scene_result_note_lines(result, changed_paths);
    }
    if (ret > 0) {
	std::vector<const char *> path_views;
	path_views.reserve(changed_paths.size());
	for (const std::string &path : changed_paths)
	    path_views.push_back(path.c_str());
	ged_scene_delta_dispatch_internal(gedp, GED_SCENE_DELTA_EDIT_SCOPE,
	    request->view, path_views.data(), path_views.size(),
	    path_views.empty(), internal_result.replaced_root_count, 0, 1,
	    revision_before, revision_after,
	    bu_vls_cstr(&internal_result.errors));
    }
    ged_scene_edit_internal_result_free(&internal_result);
    return ret < 0 ? GED_SCENE_ERROR :
	(ret > 0 ? GED_SCENE_OK : GED_SCENE_NOT_FOUND);
}


extern "C" enum ged_scene_status
ged_scene_edit_release(struct ged *gedp,
		       ged_scene_edit_scope_ref scope,
		       enum ged_scene_edit_outcome outcome,
		       struct ged_scene_result *result)
{
    if (result)
	ged_scene_result_reset(result);
    if (!gedp || ged_scene_edit_scope_ref_is_null(scope) ||
	(outcome != GED_SCENE_EDIT_CANCEL &&
	 outcome != GED_SCENE_EDIT_COMMIT)) {
	if (result)
	    result->status = GED_SCENE_INVALID;
	return GED_SCENE_INVALID;
    }

    struct ged_scene_edit_internal_result internal_result;
    ged_scene_edit_internal_result_init(&internal_result);
    const uint64_t revision_before = ged_draw_scene_revision(gedp);
    const int ret = ged_scene_edit_release_internal(gedp, scope, outcome,
	&internal_result);
    const uint64_t revision_after = ret >= 0 ?
	ged_scene_revision_advance(gedp) : ged_draw_scene_revision(gedp);
    const std::vector<std::string> changed_paths =
	ged_scene_lines(&internal_result.paths);
    if (result) {
	result->status = ret < 0 ? GED_SCENE_ERROR : GED_SCENE_OK;
	result->changed = ret >= 0;
	result->group_count = internal_result.replaced_root_count;
	result->revision_before = revision_before;
	result->revision_after = revision_after;
	result->conflict_count = internal_result.conflict_count;
	result->diagnostic = bu_vls_cstr(&internal_result.errors);
	ged_scene_result_note_lines(result, changed_paths);
    }
    if (ret >= 0) {
	std::vector<const char *> path_views;
	path_views.reserve(changed_paths.size());
	for (const std::string &path : changed_paths)
	    path_views.push_back(path.c_str());
	ged_scene_delta_dispatch_internal(gedp, GED_SCENE_DELTA_EDIT_SCOPE,
	    NULL, path_views.data(), path_views.size(), path_views.empty(),
	    internal_result.replaced_root_count, 0, 1, revision_before,
	    revision_after, bu_vls_cstr(&internal_result.errors));
    }
    ged_scene_edit_internal_result_free(&internal_result);
    return ret < 0 ? GED_SCENE_ERROR : GED_SCENE_OK;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
