/*                 G E D _ D R A W _ F R O N T I E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged_draw_frontier.cpp
 *
 * Backend-neutral draw-intent splitting used by nested erase and editing.
 *
 * A frontier replacement expands only the ancestors between an existing
 * drawn root and a target occurrence.  Sibling subtrees become independent
 * draw roots; the target is either retained (promotion) or omitted (erase).
 * Database-index occurrence ids, rather than directory pointers or path
 * strings, are the authority while planning.  This is important for duplicate
 * child instances and also keeps the operation proportional to path depth and
 * sibling count instead of total leaf count.
 */

#include "common.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "bu/malloc.h"
#include "bu/vls.h"
#include "ged/db_index.h"
#include "ged/scene.h"

#include "./ged_draw_private.h"
#include "./ged_private.h"


namespace {

using path_ids = std::vector<ged_db_index_id>;

struct visibility_rule {
    path_ids ids;
    bool visible = true;
};

struct frontier_root {
    std::string path;
    path_ids ids;
    struct ged_view_context *view = nullptr;
    int mode = -1;
    struct ged_draw_appearance_settings appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    std::vector<visibility_rule> rules;
    bool auto_draw = false;
    bool conflict = false;
};

struct promotion_record {
    uint64_t id = 0;
    uint64_t generation = 1;
    std::string role;
    std::vector<std::string> occurrences;
    std::vector<frontier_root> roots;
};

struct frontier_state {
    uint64_t owner = 0;
    uint64_t next_id = 1;
    std::map<uint64_t, promotion_record> promotions;
    std::vector<frontier_root> roots;
};


static frontier_state *
state_get(struct ged *gedp, bool create)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return nullptr;

    Ged_Internal *gi = gedp->i->i;
    frontier_state *state =
	static_cast<frontier_state *>(gi->draw_frontier_state);
    if (!state && create) {
	state = new frontier_state;
	state->owner = static_cast<uint64_t>(
	    reinterpret_cast<uintptr_t>(gedp));
	if (!state->owner)
	    state->owner = 1;
	gi->draw_frontier_state = state;
    }
    return state;
}


static bool
root_scope_equal(const frontier_root &left, const frontier_root &right)
{
    return left.path == right.path && left.view == right.view &&
	left.mode == right.mode;
}


static bool
root_in_scope(const frontier_root &root, struct ged_view_context *view,
	      int mode)
{
    if (view && root.view && view != root.view)
	return false;
    if (mode >= 0 && root.mode >= 0 && mode != root.mode)
	return false;
    return true;
}


static std::string
canonical_path(const char *path)
{
    if (!path)
	return std::string();
    while (*path == '/')
	path++;
    std::string result(path);
    while (!result.empty() && result.back() == '/')
	result.pop_back();
    return result;
}


static std::string
request_path(const struct ged_draw_promotion_request *request)
{
    if (!request)
	return std::string();
    if (request->path && request->path[0])
	return canonical_path(request->path);
    if (!request->dfp)
	return std::string();

    char *path = db_path_to_string(request->dfp);
    if (!path)
	return std::string();
    std::string result = canonical_path(path);
    bu_free(path, "draw promotion request path");
    return result;
}


static bool
path_prefix(const std::string &prefix, const std::string &path)
{
    if (prefix.empty() || path.empty())
	return false;
    return prefix == path ||
	(path.size() > prefix.size() &&
	 path.compare(0, prefix.size(), prefix) == 0 &&
	 path[prefix.size()] == '/');
}


static bool
ids_prefix(const path_ids &prefix, const path_ids &path)
{
    return prefix.size() <= path.size() &&
	std::equal(prefix.begin(), prefix.end(), path.begin());
}


static bool
resolve_path(struct ged *gedp, const std::string &path, path_ids &ids)
{
    ids.clear();
    if (!gedp || path.empty())
	return false;
    const size_t count =
	ged_db_index_path_resolve(gedp, path.c_str(), nullptr, 0);
    if (!count)
	return false;
    ids.resize(count);
    if (ged_db_index_path_resolve(gedp, path.c_str(), ids.data(),
	    ids.size()) != count) {
	ids.clear();
	return false;
    }
    return true;
}


static std::string
print_path(struct ged *gedp, const path_ids &ids)
{
    struct bu_vls path = BU_VLS_INIT_ZERO;
    if (ids.empty() ||
	!ged_db_index_path_print(gedp, &path, ids.data(), ids.size(), 0, 0)) {
	bu_vls_free(&path);
	return std::string();
    }
    std::string result = canonical_path(bu_vls_cstr(&path));
    bu_vls_free(&path);
    return result;
}


static std::set<std::string>
listed_frontier(struct ged *gedp, struct ged_view_context *view, int mode)
{
    std::set<std::string> paths;
    struct bu_vls listing = BU_VLS_INIT_ZERO;
    (void)ged_draw_list_paths(gedp, view, mode, 0, &listing);

    const char *start = bu_vls_cstr(&listing);
    for (const char *p = start; ; p++) {
	if (*p != '\n' && *p != '\r' && *p != '\0')
	    continue;
	if (p > start) {
	    std::string value(start, static_cast<size_t>(p - start));
	    value = canonical_path(value.c_str());
	    if (!value.empty())
		paths.insert(value);
	}
	if (*p == '\0')
	    break;
	start = p + 1;
    }
    bu_vls_free(&listing);
    return paths;
}


struct source_root_collect_context {
    std::set<std::string> *paths = nullptr;
};


static int
collect_source_root_cb(
    struct ged *UNUSED(gedp),
    const struct ged_draw_obol_database_source_record *record,
    void *data)
{
    source_root_collect_context *ctx =
	static_cast<source_root_collect_context *>(data);
    if (!ctx || !ctx->paths || !record || !record->valid ||
	!record->visible)
	return 1;
    const char *path = record->source_path && record->source_path[0] ?
	record->source_path : record->database_path;
    const std::string canonical = canonical_path(path);
    if (!canonical.empty())
	ctx->paths->insert(canonical);
    return 1;
}


static std::set<std::string>
listed_source_roots(struct ged *gedp)
{
    std::set<std::string> paths;
    source_root_collect_context ctx;
    ctx.paths = &paths;
    const int status =
	ged_draw_obol_visible_database_source_records_foreach_fast(
	    gedp, collect_source_root_cb, &ctx);
    if (status < 0)
	paths.clear();
    return paths;
}


struct root_collect_context {
    struct ged *gedp = nullptr;
    struct ged_view_context *view = nullptr;
    int mode = -1;
    const std::set<std::string> *listed = nullptr;
    std::vector<frontier_root> *roots = nullptr;
};


static int
collect_root_cb(const struct ged_draw_group_record *record, void *data)
{
    root_collect_context *ctx = static_cast<root_collect_context *>(data);
    if (!ctx || !record || !record->path || !ctx->listed || !ctx->roots)
	return 1;
    if (!record->visible || record->is_overlay || record->is_local_source ||
	!ged_draw_group_record_in_view(record, ctx->view))
	return 1;
    if (ctx->mode >= 0 && record->draw_mode != ctx->mode)
	return 1;

    const std::string path = canonical_path(record->path);
    if (path.empty() || ctx->listed->find(path) == ctx->listed->end())
	return 1;

    frontier_root root;
    root.path = path;
    root.view = record->view ? record->view : ctx->view;
    root.mode = record->draw_mode;
    root.appearance = record->appearance;
    if (!resolve_path(ctx->gedp, root.path, root.ids))
	return 1;

    const auto duplicate = std::find_if(ctx->roots->begin(), ctx->roots->end(),
	[&root](const frontier_root &candidate) {
	    return candidate.path == root.path &&
		candidate.view == root.view && candidate.mode == root.mode;
	});
    if (duplicate == ctx->roots->end())
	ctx->roots->push_back(root);
    return 1;
}


static std::vector<frontier_root>
current_roots(struct ged *gedp, struct ged_view_context *view, int mode)
{
    std::vector<frontier_root> roots;
    std::set<std::string> paths;
    if (ged_draw_obol_scene_controller_full_synced(gedp))
	paths = listed_source_roots(gedp);
    if (paths.empty())
	paths = listed_frontier(gedp, view, mode);
    if (!paths.empty()) {
	root_collect_context ctx;
	ctx.gedp = gedp;
	ctx.view = view;
	ctx.mode = mode;
	ctx.listed = &paths;
	ctx.roots = &roots;
	ged_draw_foreach_group_record(gedp, collect_root_cb, &ctx);
    }

    /* Once a source has a semantic frontier its physical owner root is
     * intentionally absent from the collapsed user-facing listing.  Merge the
     * retained owners back here so subsequent nested operations continue to
     * address the same source rather than manufacturing new sources. */
    frontier_state *state = state_get(gedp, false);
    if (state) {
	for (const frontier_root &retained : state->roots) {
	    if (!root_in_scope(retained, view, mode))
		continue;
	    auto found = std::find_if(roots.begin(), roots.end(),
		[&retained](const frontier_root &candidate) {
		    return root_scope_equal(candidate, retained);
		});
	    if (found == roots.end())
		roots.push_back(retained);
	    else
		*found = retained;
	}
    }
    return roots;
}


static bool
visibility_rule_less(const visibility_rule &left,
		     const visibility_rule &right)
{
    if (left.ids.size() != right.ids.size())
	return left.ids.size() < right.ids.size();
    return std::lexicographical_compare(left.ids.begin(), left.ids.end(),
	right.ids.begin(), right.ids.end());
}


static bool
frontier_path_visible(const frontier_root &root, const path_ids &target)
{
    bool visible = true;
    size_t matched_depth = 0;
    for (const visibility_rule &rule : root.rules) {
	if (rule.ids.size() >= matched_depth && ids_prefix(rule.ids, target)) {
	    visible = rule.visible;
	    matched_depth = rule.ids.size();
	}
    }
    return visible;
}


/* Replace the visibility state of a complete subtree.  Descendant rules are
 * discarded because the new draw/erase intent covers the whole subtree.  An
 * override is stored only when it differs from the inherited ancestor state,
 * keeping memory and update work proportional to user operations rather than
 * to the number of siblings in the database. */
static bool
frontier_set_subtree_visibility(frontier_root &root,
				const path_ids &target,
				bool visible)
{
    if (!ids_prefix(root.ids, target))
	return false;

    const std::vector<visibility_rule> previous = root.rules;
    root.rules.erase(std::remove_if(root.rules.begin(), root.rules.end(),
	[&target](const visibility_rule &rule) {
	    return ids_prefix(target, rule.ids);
	}), root.rules.end());

    bool inherited = true;
    size_t inherited_depth = 0;
    for (const visibility_rule &rule : root.rules) {
	if (rule.ids.size() >= inherited_depth &&
	    ids_prefix(rule.ids, target)) {
	    inherited = rule.visible;
	    inherited_depth = rule.ids.size();
	}
    }
    if (visible != inherited)
	root.rules.push_back(visibility_rule{target, visible});

    std::stable_sort(root.rules.begin(), root.rules.end(),
	visibility_rule_less);
    if (previous.size() != root.rules.size())
	return true;
    for (size_t i = 0; i < previous.size(); i++) {
	if (previous[i].ids != root.rules[i].ids ||
	    previous[i].visible != root.rules[i].visible)
	    return true;
    }
    return false;
}


/* Database hierarchy edits may retire occurrence ids recorded by visibility
 * rules.  Keeping such a rule would leave its owner permanently "partial"
 * even though the excluded occurrence no longer exists. */
static bool
frontier_reconcile_database_paths(struct ged *gedp, frontier_root &root)
{
    path_ids resolved_root;
    if (!resolve_path(gedp, root.path, resolved_root))
	return false;
    root.ids.swap(resolved_root);

    root.rules.erase(std::remove_if(root.rules.begin(), root.rules.end(),
	[gedp, &root](const visibility_rule &rule) {
	    const std::string path = print_path(gedp, rule.ids);
	    if (path.empty())
		return true;
	    path_ids resolved;
	    return !resolve_path(gedp, path, resolved) ||
		resolved != rule.ids || !ids_prefix(root.ids, resolved);
	}), root.rules.end());
    std::stable_sort(root.rules.begin(), root.rules.end(),
	visibility_rule_less);
    return true;
}


static bool
frontier_appearance_equal(
    const struct ged_draw_appearance_settings &left,
    const struct ged_draw_appearance_settings &right)
{
    return left.draw_mode == right.draw_mode &&
	std::fabs(left.transparency - right.transparency) <= 1.0e-12 &&
	left.color_override == right.color_override &&
	(!left.color_override ||
	 (left.color[0] == right.color[0] &&
	  left.color[1] == right.color[1] &&
	  left.color[2] == right.color[2])) &&
	left.s_line_width == right.s_line_width &&
	left.draw_solid_lines_only == right.draw_solid_lines_only &&
	left.draw_non_subtract_only == right.draw_non_subtract_only;
}


static int
retained_root_apply(struct ged *gedp, frontier_root &root)
{
    frontier_state *state = state_get(gedp, true);
    if (!state)
	return -1;

    std::vector<std::string> paths;
    std::vector<int> visible;
    paths.reserve(root.rules.size());
    visible.reserve(root.rules.size());
    for (const visibility_rule &rule : root.rules) {
	const std::string path = print_path(gedp, rule.ids);
	if (path.empty())
	    continue;
	paths.push_back(path);
	visible.push_back(rule.visible ? 1 : 0);
    }
    std::vector<const char *> args;
    args.reserve(paths.size());
    for (const std::string &path : paths)
	args.push_back(path.c_str());

    int renderer_changed = 0;
    if (paths.empty()) {
	renderer_changed = ged_draw_obol_source_visibility_frontier_clear(
	    gedp, root.path.c_str(), root.view, root.mode);
    } else {
	renderer_changed = ged_draw_obol_source_visibility_overrides_set(
	    gedp, root.path.c_str(), root.view, root.mode,
	    args.data(), visible.data(), args.size());
    }
    (void)renderer_changed;

    auto found = std::find_if(state->roots.begin(), state->roots.end(),
	[&root](const frontier_root &candidate) {
	    return root_scope_equal(candidate, root);
	});
    if (root.rules.empty()) {
	if (found != state->roots.end())
	    state->roots.erase(found);
	return 1;
    }

    if (found == state->roots.end())
	state->roots.push_back(root);
    else
	*found = root;
    return 1;
}


static int
draw_paths_direct(struct ged *gedp, const frontier_root &root,
		  const std::vector<std::string> &paths)
{
    if (paths.empty())
	return 0;

    std::vector<const char *> path_args;
    path_args.reserve(paths.size());
    for (const std::string &path : paths)
	path_args.push_back(path.c_str());

    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, nullptr);
    txn.view = root.view;
    txn.mode = root.mode;
    txn.paths = path_args.data();
    txn.path_count = static_cast<int>(path_args.size());
    txn.appearance = &root.appearance;
    txn.autoview = 0;
    return ged_draw_apply_draw_inner(gedp, &txn, nullptr, nullptr);
}


static std::vector<path_ids>
promotion_occurrences(struct ged *gedp, const path_ids &requested, int scope)
{
    std::vector<path_ids> occurrences;
    if (requested.empty())
	return occurrences;
    if (scope != GED_DRAW_PROMOTE_ALL_OBJECT_OCCURRENCES) {
	occurrences.push_back(requested);
	return occurrences;
    }

    struct ged_db_index_record leaf;
    if (!ged_db_index_record_get(gedp, requested.back(), &leaf) ||
	!leaf.valid)
	return occurrences;
    const size_t count = ged_db_index_affected_path_count(
	gedp, leaf.object_id, 0);
    for (size_t row = 0; row < count; row++) {
	const size_t depth = ged_db_index_affected_path_at(
	    gedp, leaf.object_id, row, nullptr, 0, 0);
	if (!depth)
	    continue;
	path_ids ids(depth);
	if (ged_db_index_affected_path_at(gedp, leaf.object_id, row,
		ids.data(), ids.size(), 0) == depth)
	    occurrences.push_back(ids);
    }
    return occurrences;
}


static void
append_result_path(struct ged_draw_promotion_result *result,
		   const std::string &path)
{
    if (!result || path.empty())
	return;
    if (bu_vls_strlen(&result->paths))
	bu_vls_putc(&result->paths, '\n');
    bu_vls_strcat(&result->paths, path.c_str());
}


static void
result_prepare(struct ged_draw_promotion_result *result, struct ged *gedp)
{
    if (!result)
	return;
    if (!BU_VLS_IS_INITIALIZED(&result->paths)) {
	BU_VLS_INIT(&result->paths);
    } else {
	bu_vls_trunc(&result->paths, 0);
    }
    if (!BU_VLS_IS_INITIALIZED(&result->errors)) {
	BU_VLS_INIT(&result->errors);
    } else {
	bu_vls_trunc(&result->errors, 0);
    }
    result->status = 0;
    result->promotion = GED_DRAW_PROMOTION_REF_NULL;
    result->occurrence_count = 0;
    result->replaced_root_count = 0;
    result->replacement_path_count = 0;
    result->conflict_count = 0;
    result->scene_revision_before = ged_draw_scene_revision(gedp);
    result->scene_revision_after = result->scene_revision_before;
}

} // namespace


extern "C" void
ged_draw_frontier_state_destroy(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;
    Ged_Internal *gi = gedp->i->i;
    frontier_state *state =
	static_cast<frontier_state *>(gi->draw_frontier_state);
    delete state;
    gi->draw_frontier_state = nullptr;
}


extern "C" int
ged_draw_promotion_ref_is_null(ged_draw_promotion_ref promotion)
{
    return (!promotion.owner || !promotion.id || !promotion.generation) ? 1 : 0;
}


extern "C" int
ged_draw_frontier_list_paths(struct ged *gedp,
			     struct ged_view_context *view_ctx,
			     int mode,
			     struct bu_vls *result,
			     size_t result_start)
{
    if (!gedp || !result || result_start > bu_vls_strlen(result))
	return -1;
    frontier_state *state = state_get(gedp, false);
    if (!state || state->roots.empty())
	return -1;

    std::set<std::string> paths;
    const char *all = bu_vls_cstr(result);
    const char *start = all + result_start;
    for (const char *p = start; ; p++) {
	if (*p != '\n' && *p != '\r' && *p != '\0')
	    continue;
	if (p > start)
	    paths.insert(canonical_path(
		std::string(start, static_cast<size_t>(p - start)).c_str()));
	if (*p == '\0')
	    break;
	start = p + 1;
    }

    for (const frontier_root &root : state->roots) {
	if (!root_in_scope(root, view_ctx, mode))
	    continue;
	for (auto it = paths.begin(); it != paths.end(); ) {
	    if (path_prefix(root.path, *it) ||
		path_prefix(*it, root.path))
		it = paths.erase(it);
	    else
		++it;
	}
	/* A retained source remains the user's draw root even when it has a
	 * compact set of hidden/re-included subtrees.  Detailed callers can ask
	 * ged_draw_path_state for full/partial/hidden state without expanding one
	 * exclusion into every visible sibling here. */
	paths.insert(root.path);
    }

    bu_vls_trunc(result, (ssize_t)result_start);
    for (const std::string &path : paths)
	bu_vls_printf(result, "%s\n", path.c_str());
    return paths.size() > static_cast<size_t>(INT_MAX) ? INT_MAX :
	static_cast<int>(paths.size());
}


extern "C" int
ged_draw_frontier_path_state(struct ged *gedp,
			     struct ged_view_context *view_ctx,
			     const char *path,
			     int mode)
{
    frontier_state *state = state_get(gedp, false);
    if (!state || !path || !path[0])
	return -1;

    path_ids target;
    if (!resolve_path(gedp, canonical_path(path), target))
	return -1;

    bool matched = false;
    bool fully_visible = false;
    bool partially_visible = false;
    for (const frontier_root &root : state->roots) {
	if (!root_in_scope(root, view_ctx, mode) ||
	    !ids_prefix(root.ids, target))
	    continue;
	matched = true;
	const bool visible = frontier_path_visible(root, target);
	bool partial = false;
	for (const visibility_rule &rule : root.rules) {
	    if (rule.ids.size() > target.size() &&
		ids_prefix(target, rule.ids) &&
		rule.visible != visible) {
		partial = true;
		break;
	    }
	}
	if (partial)
	    partially_visible = true;
	else if (visible)
	    fully_visible = true;
    }
    if (!matched)
	return -1;
    if (fully_visible)
	return 1;
    if (partially_visible)
	return 2;
    return 0;
}


extern "C" void
ged_draw_frontier_note_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const char *resolved_path)
{
    frontier_state *state = state_get(gedp, false);
    if (!state || !txn)
	return;

    bool structural = false;
    bool all_roots = false;
    switch (txn->kind) {
	case GED_DRAW_TXN_DRAW:
	case GED_DRAW_TXN_ERASE:
	case GED_DRAW_TXN_ERASE_PREFIX:
	case GED_DRAW_TXN_SOURCE_RENAMED:
	case GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
	    structural = true;
	    break;
	case GED_DRAW_TXN_SOURCE_UPDATED:
	    structural = txn->removed != 0 || txn->redraw != 0;
	    break;
	case GED_DRAW_TXN_CLEAR:
	case GED_DRAW_TXN_CLEAR_SCOPE:
	case GED_DRAW_TXN_TEARDOWN:
	    structural = true;
	    all_roots = true;
	    break;
	default:
	    break;
    }
    if (!structural)
	return;

    std::vector<std::string> paths;
    if (resolved_path && resolved_path[0])
	paths.push_back(canonical_path(resolved_path));
    if (txn->new_path && txn->new_path[0])
	paths.push_back(canonical_path(txn->new_path));
    if (txn->paths && txn->path_count > 0) {
	for (int i = 0; i < txn->path_count; i++) {
	    if (txn->paths[i] && txn->paths[i][0])
		paths.push_back(canonical_path(txn->paths[i]));
	}
    }

    for (auto &promotion_entry : state->promotions) {
	for (frontier_root &root : promotion_entry.second.roots) {
	    if (!root_in_scope(root, txn->view, txn->mode))
		continue;
	    if (all_roots) {
		root.conflict = true;
		continue;
	    }
	    for (const std::string &path : paths) {
		if (path_prefix(root.path, path) ||
		    path_prefix(path, root.path)) {
		    root.conflict = true;
		    break;
		}
	    }
	}
    }

    /* Reconcile compact occurrence rules against the updated database index
     * before a source redraw.  In particular, a nested erase issued as part
     * of a combination-member removal must not survive after that occurrence
     * has ceased to exist. */
    if (txn->kind == GED_DRAW_TXN_SOURCE_UPDATED && txn->redraw) {
	for (auto it = state->roots.begin(); it != state->roots.end(); ) {
	    if (!root_in_scope(*it, txn->view, txn->mode)) {
		++it;
		continue;
	    }
	    if (!frontier_reconcile_database_paths(gedp, *it) ||
		it->rules.empty()) {
		(void)ged_draw_obol_source_visibility_frontier_clear(gedp,
		    it->path.c_str(), it->view, it->mode);
		it = state->roots.erase(it);
		continue;
	    }
	    ++it;
	}
    }

    /* Removing an owning root or clearing its scope retires its presentation
     * frontier as well.  Descendant erases are handled by the frontier planner
     * after this notification. */
    for (auto it = state->roots.begin(); it != state->roots.end(); ) {
	if (!root_in_scope(*it, txn->view, txn->mode)) {
	    ++it;
	    continue;
	}
	bool retire = all_roots;
	if (!retire && (txn->kind == GED_DRAW_TXN_ERASE ||
		txn->kind == GED_DRAW_TXN_ERASE_PREFIX ||
		txn->kind == GED_DRAW_TXN_DRAW)) {
	    for (const std::string &path : paths) {
		/* Erasing an owner, or explicitly drawing a broader owner, retires
		 * the narrower source's presentation rules.  A descendant draw is
		 * absorbed by ged_draw_frontier_absorb_draw and keeps the owner. */
		if (path_prefix(path, it->path)) {
		    retire = true;
		    break;
		}
	    }
	}
	if (!retire) {
	    ++it;
	    continue;
	}
	(void)ged_draw_obol_source_visibility_frontier_clear(gedp,
	    it->path.c_str(), it->view, it->mode);
	it = state->roots.erase(it);
    }
}


extern "C" int
ged_draw_frontier_absorb_draw(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const char *resolved_path,
    struct ged_draw_transaction_result *result)
{
    if (!gedp || !txn || txn->kind != GED_DRAW_TXN_DRAW ||
	!ged_db_index_available(gedp))
	return 0;

    const struct ged_draw_appearance_settings default_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    const struct ged_draw_appearance_settings &appearance = txn->appearance ?
	*static_cast<const struct ged_draw_appearance_settings *>(
	    txn->appearance) : default_appearance;
    const int mode = txn->mode >= 0 ? txn->mode : appearance.draw_mode;

    std::vector<std::string> requested_paths;
    if (txn->paths && txn->path_count > 0) {
	for (int i = 0; i < txn->path_count; i++) {
	    const std::string path = canonical_path(txn->paths[i]);
	    if (!path.empty())
		requested_paths.push_back(path);
	}
    } else {
	const std::string path = canonical_path(resolved_path);
	if (!path.empty())
	    requested_paths.push_back(path);
    }
    if (requested_paths.empty())
	return 0;

    std::vector<frontier_root> roots =
	current_roots(gedp, txn->view, mode);
    std::vector<bool> changed(roots.size(), false);
    for (const std::string &path : requested_paths) {
	path_ids target;
	if (!resolve_path(gedp, path, target))
	    return 0;

	size_t owner = roots.size();
	for (size_t i = 0; i < roots.size(); i++) {
	    if (roots[i].mode != mode ||
		!frontier_appearance_equal(roots[i].appearance, appearance) ||
		!ids_prefix(roots[i].ids, target))
		continue;
	    if (owner == roots.size() ||
		roots[i].ids.size() > roots[owner].ids.size())
		owner = i;
	}
	if (owner == roots.size())
	    return 0;

	if (frontier_set_subtree_visibility(roots[owner], target, true))
	    changed[owner] = true;
    }

    int changed_roots = 0;
    for (size_t i = 0; i < roots.size(); i++) {
	if (!changed[i])
	    continue;
	if (retained_root_apply(gedp, roots[i]) < 0)
	    return -1;
	changed_roots++;
    }
    if (changed_roots)
	(void)ged_selection_draw_sync(gedp, nullptr);
    if (result) {
	result->affected_groups +=
	    static_cast<int>(requested_paths.size());
	result->affected_shapes +=
	    static_cast<int>(requested_paths.size());
	result->presentation_only = 1;
    }
    return static_cast<int>(requested_paths.size());
}


extern "C" void
ged_draw_promotion_result_init(struct ged_draw_promotion_result *result)
{
    if (!result)
	return;
    *result = {};
    BU_VLS_INIT(&result->paths);
    BU_VLS_INIT(&result->errors);
    result->promotion = GED_DRAW_PROMOTION_REF_NULL;
}


extern "C" void
ged_draw_promotion_result_free(struct ged_draw_promotion_result *result)
{
    if (!result)
	return;
    if (BU_VLS_IS_INITIALIZED(&result->paths))
	bu_vls_free(&result->paths);
    if (BU_VLS_IS_INITIALIZED(&result->errors))
	bu_vls_free(&result->errors);
    *result = {};
}


extern "C" int
ged_draw_promote_path(struct ged *gedp,
		      const struct ged_draw_promotion_request *request,
		      struct ged_draw_promotion_result *result)
{
    result_prepare(result, gedp);
    if (!gedp || !request || !ged_db_index_available(gedp)) {
	if (result)
	    bu_vls_strcat(&result->errors,
		"draw promotion requires a database index");
	return -1;
    }

    const std::string target_path = request_path(request);
    path_ids requested_ids;
    if (target_path.empty() ||
	!resolve_path(gedp, target_path, requested_ids)) {
	if (result)
	    bu_vls_printf(&result->errors,
		"cannot resolve promotion path: %s", target_path.c_str());
	return -1;
    }

    std::vector<path_ids> occurrences =
	promotion_occurrences(gedp, requested_ids, request->scope);
    if (occurrences.empty())
	occurrences.push_back(requested_ids);

    std::vector<frontier_root> roots =
	current_roots(gedp, request->view, request->mode);
    promotion_record promotion;
    promotion.role = request->role ? request->role : "";

    std::set<std::string> handled_occurrences;
    std::vector<bool> occurrence_drawn(occurrences.size(), false);
    for (frontier_root &root : roots) {
	bool root_has_target = false;
	for (size_t occurrence_index = 0;
	     occurrence_index < occurrences.size(); occurrence_index++) {
	    const path_ids &occurrence = occurrences[occurrence_index];
	    if (!ids_prefix(root.ids, occurrence) ||
		!frontier_path_visible(root, occurrence))
		continue;
	    root_has_target = true;
	    occurrence_drawn[occurrence_index] = true;
	    const std::string occurrence_path = print_path(gedp, occurrence);
	    if (!occurrence_path.empty())
		handled_occurrences.insert(occurrence_path);
	}
	if (!root_has_target)
	    continue;

	/* Promotion records an edit boundary without manufacturing sibling draw
	 * roots or changing the retained source's visibility.  The edit overlay
	 * addresses the returned occurrence paths while the source continues to
	 * own geometry, LoD state, and cache residency. */
	promotion.roots.push_back(root);
	if (result) {
	    result->replaced_root_count++;
	}
    }

    for (size_t occurrence_index = 0;
	 occurrence_index < occurrences.size(); occurrence_index++) {
	const path_ids &occurrence = occurrences[occurrence_index];
	/* "All occurrences" means all occurrences already represented by the
	 * user's draw intents.  It must not make every database occurrence
	 * visible.  Auto-draw applies only to the explicitly requested path. */
	if (!occurrence_drawn[occurrence_index] && request->auto_draw &&
	    occurrence == requested_ids) {
	    const std::string occurrence_path = print_path(gedp, occurrence);
	    if (occurrence_path.empty())
		continue;
	    frontier_root root;
	    root.path = occurrence_path;
	    root.ids = occurrence;
	    root.view = request->view;
	    root.mode = request->mode >= 0 ? request->mode :
		ged_draw_default_mode(gedp);
	    const struct ged_draw_appearance_settings default_appearance =
		GED_DRAW_APPEARANCE_SETTINGS_INIT;
	    root.appearance = default_appearance;
	    root.appearance.draw_mode = root.mode;
	    root.auto_draw = true;
	    std::vector<std::string> replacement(1, occurrence_path);
	    if (draw_paths_direct(gedp, root, replacement) < 0) {
		if (result)
		    bu_vls_printf(&result->errors,
			"failed to auto-draw promotion target %s",
			occurrence_path.c_str());
		return -1;
	    }
	    promotion.roots.push_back(root);
	    handled_occurrences.insert(occurrence_path);
	    if (result) {
		result->replaced_root_count++;
		result->replacement_path_count++;
	    }
	}
    }

    if (promotion.roots.empty()) {
	if (result)
	    bu_vls_printf(&result->errors,
		"no drawn owner contains promotion path %s",
		target_path.c_str());
	return 0;
    }

    frontier_state *state = state_get(gedp, true);
    if (!state)
	return -1;
    promotion.id = state->next_id++;
    if (!promotion.id)
	promotion.id = state->next_id++;
    promotion.occurrences.assign(handled_occurrences.begin(),
	handled_occurrences.end());
    if (result)
	result->replacement_path_count = promotion.occurrences.size();
    state->promotions[promotion.id] = promotion;

    ged_draw_promotion_ref ref = {
	state->owner, promotion.id, promotion.generation
    };
    (void)ged_selection_draw_sync(gedp, nullptr);
    if (result) {
	result->status = 1;
	result->promotion = ref;
	result->occurrence_count = promotion.occurrences.size();
	for (const std::string &path : promotion.occurrences)
	    append_result_path(result, path);
	result->scene_revision_after = ged_draw_scene_revision(gedp);
    }
    return 1;
}


extern "C" int
ged_draw_release_promotion(struct ged *gedp,
			   ged_draw_promotion_ref ref,
			   int outcome,
			   struct ged_draw_promotion_result *result)
{
    result_prepare(result, gedp);
    frontier_state *state = state_get(gedp, false);
    if (!state || ref.owner != state->owner || !ref.id ||
	!ref.generation) {
	if (result)
	    bu_vls_strcat(&result->errors, "unknown draw promotion");
	return -1;
    }
    auto found = state->promotions.find(ref.id);
    if (found == state->promotions.end() ||
	found->second.generation != ref.generation) {
	if (result)
	    bu_vls_strcat(&result->errors, "stale draw promotion");
	return -1;
    }

    promotion_record promotion = found->second;
    state->promotions.erase(found);
    size_t conflicts = 0;
    for (const frontier_root &root : promotion.roots) {
	if (root.conflict)
	    conflicts++;
    }

    int changed = 0;
    if (!conflicts) {
	for (frontier_root &root : promotion.roots) {
	    if (root.auto_draw) {
		if (ged_draw_erase_path_string_scoped(gedp,
			root.path.c_str(), root.view, root.mode) < 0) {
		    if (result)
			bu_vls_printf(&result->errors,
			    "failed to remove auto-drawn promotion root %s",
			root.path.c_str());
		    return -1;
		}
	    }
	    /* Non-auto-drawn promotions are logical edit boundaries.  Releasing
	     * one does not rebuild or remask its retained source. */
	    changed = 1;
	}
    }

    (void)outcome;
    (void)ged_selection_draw_sync(gedp, nullptr);
    if (result) {
	result->status = conflicts ? 0 : changed;
	result->promotion = ref;
	result->occurrence_count = promotion.occurrences.size();
	result->replaced_root_count = promotion.roots.size();
	result->conflict_count = conflicts;
	for (const std::string &path : promotion.occurrences)
	    append_result_path(result, path);
	if (conflicts)
	    bu_vls_printf(&result->errors,
		"promotion frontier changed under %zu root(s); current draw state preserved",
		conflicts);
	result->scene_revision_after = ged_draw_scene_revision(gedp);
    }
    return conflicts ? 0 : changed;
}


extern "C" int
ged_draw_frontier_erase_path(struct ged *gedp, const char *path,
			     struct ged_view_context *view_ctx,
			     int mode, int prefix,
			     struct ged_draw_transaction_result *result)
{
    if (!gedp || !path || !path[0] || !ged_db_index_available(gedp))
	return 0;

    const std::string target_path = canonical_path(path);
    path_ids target;
    if (!resolve_path(gedp, target_path, target))
	return 0;

    std::vector<frontier_root> roots =
	current_roots(gedp, view_ctx, mode);
    std::vector<frontier_root> plans;
    for (frontier_root &root : roots) {
	if (root.ids == target)
	    continue; /* the ordinary exact erase already handles this */
	if (!ids_prefix(root.ids, target))
	    continue;
	if (frontier_set_subtree_visibility(root, target, false))
	    plans.push_back(root);
    }
    if (plans.empty())
	return 0;

    int changed = 0;
    for (frontier_root &plan : plans) {
	if (retained_root_apply(gedp, plan) < 0)
	    return -1;
	changed++;
	if (result) {
	    result->affected_groups++;
	    result->affected_shapes++;
	}
    }
    (void)prefix;
    (void)ged_selection_draw_sync(gedp, nullptr);
    if (result)
	result->presentation_only = 1;
    return changed;
}
