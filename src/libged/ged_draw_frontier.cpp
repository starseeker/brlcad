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

struct presentation_rule {
    ged_scene_reducer_operation kind = GED_SCENE_REDUCER_NONE;
    std::string path;
    struct ged_view_context *view = nullptr;
    int mode = -1;
    enum ged_scene_path_match match = GED_SCENE_PATH_MATCH_EXACT;
    double value = 0.0;
};

struct frontier_root {
    std::string path;
    path_ids ids;
    struct ged_view_context *view = nullptr;
    int mode = -1;
    struct ged_draw_appearance_settings appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    std::vector<visibility_rule> rules;
    std::vector<presentation_rule> presentation;
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

struct frontier_backend_change {
    frontier_root root;
    bool clear = false;
};

struct frontier_state {
    uint64_t owner = 0;
    uint64_t next_id = 1;
    std::map<uint64_t, promotion_record> promotions;
    std::vector<frontier_root> roots;
    std::vector<frontier_backend_change> backend_changes;
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


static void
frontier_backend_change_note(frontier_state *state,
			     const frontier_root &root, bool clear)
{
    if (!state)
	return;
    auto found = std::find_if(state->backend_changes.begin(),
	state->backend_changes.end(),
	[&root](const frontier_backend_change &candidate) {
	    return root_scope_equal(candidate.root, root);
	});
    frontier_backend_change change;
    change.root = root;
    change.clear = clear;
    if (found == state->backend_changes.end())
	state->backend_changes.push_back(change);
    else
	*found = change;
}


static bool
root_in_scope(const frontier_root &root, struct ged_view_context *view,
	      int mode)
{
    /* A non-independent view is a viewport onto the shared scene, not a
     * distinct semantic draw scope.  Store and compare all shared scopes as
     * NULL so a draw issued from one shared view is visible in every other
     * shared view, while independent scopes remain strictly view-local. */
    struct ged_view_context *scope_view =
	(view && ged_view_context_is_independent(view)) ? view : nullptr;
    if (root.view != scope_view)
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
path_last_component(const std::string &path)
{
    const size_t slash = path.rfind('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}


static std::string
replace_path_component(const std::string &path, const std::string &old_name,
	const std::string &new_name, bool *changed)
{
    if (changed)
	*changed = false;
    if (path.empty() || old_name.empty() || new_name.empty())
	return path;

    std::string result;
    size_t offset = 0;
    while (offset < path.size()) {
	const size_t slash = path.find('/', offset);
	const size_t end = slash == std::string::npos ? path.size() : slash;
	const std::string component = path.substr(offset, end - offset);
	if (!result.empty())
	    result.push_back('/');
	if (component == old_name) {
	    result.append(new_name);
	    if (changed)
		*changed = true;
	} else {
	    result.append(component);
	}
	if (slash == std::string::npos)
	    break;
	offset = slash + 1;
    }
    return result;
}


static bool
path_has_component(const std::string &path, const std::string &component)
{
    if (path.empty() || component.empty())
	return false;
    size_t offset = 0;
    while (offset < path.size()) {
	const size_t slash = path.find('/', offset);
	const size_t end = slash == std::string::npos ? path.size() : slash;
	if (path.compare(offset, end - offset, component) == 0)
	    return true;
	if (slash == std::string::npos)
	    break;
	offset = slash + 1;
    }
    return false;
}


static std::string
request_path(const struct ged_scene_edit_request *request)
{
    if (!request)
	return std::string();
    return request->path && request->path[0] ?
	canonical_path(request->path) : std::string();
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


static bool resolve_path(struct ged *gedp, const std::string &path,
	path_ids &ids);


static bool
object_occurs_under_root(struct ged *gedp, const frontier_root &root,
			 const std::string &path)
{
    path_ids query;
    if (!resolve_path(gedp, path, query) || query.empty())
	return false;
    struct ged_db_index_record leaf;
    if (!ged_db_index_record_get(gedp, query.back(), &leaf) || !leaf.valid)
	return false;
    const size_t count = ged_db_index_affected_path_count(gedp,
	leaf.object_id, 0);
    for (size_t row = 0; row < count; row++) {
	const size_t depth = ged_db_index_affected_path_at(gedp,
	    leaf.object_id, row, nullptr, 0, 0);
	if (!depth || depth < root.ids.size())
	    continue;
	path_ids occurrence(depth);
	if (ged_db_index_affected_path_at(gedp, leaf.object_id, row,
		occurrence.data(), occurrence.size(), 0) == depth &&
	    ids_prefix(root.ids, occurrence))
	    return true;
    }
    return false;
}


static bool
presentation_root_matches(struct ged *gedp, const frontier_root &root,
			  const presentation_rule &rule)
{
    if (!root_in_scope(root, rule.view, rule.mode))
	return false;
    if (rule.match == GED_SCENE_PATH_MATCH_OBJECT)
	return object_occurs_under_root(gedp, root, rule.path);
    return path_prefix(root.path, rule.path) ||
	path_prefix(rule.path, root.path);
}


static bool
presentation_rule_equal(const presentation_rule &left,
			const presentation_rule &right)
{
    return left.kind == right.kind && left.path == right.path &&
	left.view == right.view && left.mode == right.mode &&
	left.match == right.match &&
	std::fabs(left.value - right.value) <= 1.0e-12;
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


static std::vector<frontier_root>
current_roots(struct ged *gedp, struct ged_view_context *view, int mode)
{
    std::vector<frontier_root> roots;
    frontier_state *state = state_get(gedp, false);
    if (state) {
	for (const frontier_root &retained : state->roots) {
	    if (root_in_scope(retained, view, mode))
		roots.push_back(retained);
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
	left.draw_non_subtract_only == right.draw_non_subtract_only &&
	left.mixed_modes == right.mixed_modes &&
	left.strict_fallback == right.strict_fallback &&
	left.defer_leaf_expansion == right.defer_leaf_expansion;
}


static int
retained_root_apply(struct ged *gedp, frontier_root &root)
{
    frontier_state *state = state_get(gedp, true);
    if (!state)
	return -1;

    auto found = std::find_if(state->roots.begin(), state->roots.end(),
	[&root](const frontier_root &candidate) {
	    return root_scope_equal(candidate, root);
	});
    if (found == state->roots.end())
	state->roots.push_back(root);
    else
	*found = root;
    frontier_backend_change_note(state, root, root.rules.empty());
    return 1;
}


static std::vector<path_ids>
promotion_occurrences(struct ged *gedp, const path_ids &requested, int scope)
{
    std::vector<path_ids> occurrences;
    if (requested.empty())
	return occurrences;
    if (scope != GED_SCENE_EDIT_ALL_DRAWN_OCCURRENCES) {
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
append_result_path(struct ged_scene_edit_internal_result *result,
		   const std::string &path)
{
    if (!result || path.empty())
	return;
    if (bu_vls_strlen(&result->paths))
	bu_vls_putc(&result->paths, '\n');
    bu_vls_strcat(&result->paths, path.c_str());
}


static void
result_prepare(struct ged_scene_edit_internal_result *result, struct ged *gedp)
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
    result->scope = GED_SCENE_EDIT_SCOPE_REF_NULL;
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
ged_draw_frontier_has_roots(struct ged *gedp,
			    struct ged_view_context *view_ctx,
			    int mode)
{
    frontier_state *state = state_get(gedp, false);
    if (!state)
	return 0;
    for (const frontier_root &root : state->roots) {
	if (root_in_scope(root, view_ctx, mode))
	    return 1;
    }
    return 0;
}


extern "C" int
ged_draw_frontier_root_count(struct ged *gedp)
{
    frontier_state *state = state_get(gedp, false);
    if (!state)
	return 0;
    return state->roots.size() > static_cast<size_t>(INT_MAX) ? INT_MAX :
	static_cast<int>(state->roots.size());
}


extern "C" int
ged_draw_frontier_clear(struct ged *gedp,
	struct ged_view_context *view_ctx, int mode)
{
    frontier_state *state = state_get(gedp, false);
    if (!state)
	return 0;

    int changed = 0;
    for (auto it = state->roots.begin(); it != state->roots.end(); ) {
	if (view_ctx && !root_in_scope(*it, view_ctx, mode)) {
	    ++it;
	    continue;
	}
	frontier_backend_change_note(state, *it, true);
	it = state->roots.erase(it);
	changed++;
    }
    return changed;
}


extern "C" int
ged_draw_frontier_source_affected(
    struct ged *gedp, const char *path, const char *const *paths,
    size_t path_count, struct ged_view_context *view_ctx, int mode)
{
    frontier_state *state = state_get(gedp, false);
    if (!state || state->roots.empty())
	return 0;

    std::vector<std::string> targets;
    if (path && path[0])
	targets.push_back(canonical_path(path));
    for (size_t i = 0; paths && i < path_count; i++) {
	if (paths[i] && paths[i][0])
	    targets.push_back(canonical_path(paths[i]));
    }
    if (targets.empty())
	return ged_draw_frontier_root_count(gedp);

    int affected = 0;
    for (const frontier_root &root : state->roots) {
	if ((view_ctx || mode >= 0) &&
	    !root_in_scope(root, view_ctx, mode))
	    continue;
	bool match = false;
	for (const std::string &target : targets) {
	    if (target.empty())
		continue;
	    if (path_prefix(root.path, target) ||
		path_prefix(target, root.path) ||
		path_has_component(root.path, path_last_component(target)) ||
		object_occurs_under_root(gedp, root, target)) {
		match = true;
		break;
	    }
	}
	if (match)
	    affected++;
    }
    return affected;
}


static bool
rename_frontier_root(struct ged *gedp, frontier_root &root,
	const std::string &old_name, const std::string &new_name)
{
    bool root_changed = false;
    const std::string updated_root = replace_path_component(root.path,
	old_name, new_name, &root_changed);
    if (root_changed) {
	root.path = updated_root;
	path_ids updated_ids;
	if (resolve_path(gedp, root.path, updated_ids))
	    root.ids.swap(updated_ids);
    }

    std::vector<visibility_rule> updated_rules;
    updated_rules.reserve(root.rules.size());
    for (const visibility_rule &rule : root.rules) {
	const std::string rule_path = print_path(gedp, rule.ids);
	bool rule_changed = false;
	const std::string updated_path = replace_path_component(rule_path,
	    old_name, new_name, &rule_changed);
	if (!rule_changed) {
	    updated_rules.push_back(rule);
	    continue;
	}
	path_ids updated_ids;
	if (resolve_path(gedp, updated_path, updated_ids))
	    updated_rules.push_back(visibility_rule{updated_ids, rule.visible});
	root_changed = true;
    }
    root.rules.swap(updated_rules);
    std::stable_sort(root.rules.begin(), root.rules.end(),
	visibility_rule_less);

    for (presentation_rule &rule : root.presentation) {
	bool rule_changed = false;
	rule.path = replace_path_component(rule.path, old_name, new_name,
	    &rule_changed);
	root_changed = root_changed || rule_changed;
    }
    return root_changed;
}


extern "C" int
ged_draw_frontier_source_rename(struct ged *gedp, const char *old_path,
	const char *new_path)
{
    frontier_state *state = state_get(gedp, false);
    if (!state || !old_path || !old_path[0] || !new_path || !new_path[0])
	return 0;

    const std::string old_name = path_last_component(canonical_path(old_path));
    const std::string new_name = path_last_component(canonical_path(new_path));
    int changed = 0;
    for (frontier_root &root : state->roots) {
	if (!rename_frontier_root(gedp, root, old_name, new_name))
	    continue;
	frontier_backend_change_note(state, root, root.rules.empty());
	changed++;
    }
    for (auto &entry : state->promotions) {
	promotion_record &promotion = entry.second;
	for (std::string &occurrence : promotion.occurrences) {
	    bool occurrence_changed = false;
	    occurrence = replace_path_component(occurrence, old_name,
		new_name, &occurrence_changed);
	    if (occurrence_changed) {
		for (frontier_root &root : promotion.roots)
		    root.conflict = true;
	    }
	}
	for (frontier_root &root : promotion.roots) {
	    if (rename_frontier_root(gedp, root, old_name, new_name))
		root.conflict = true;
	}
    }
    return changed;
}


extern "C" int
ged_draw_frontier_foreach_root(struct ged *gedp,
			       struct ged_view_context *view_ctx,
			       int mode,
			       ged_draw_frontier_root_cb callback,
			       void *userdata)
{
    if (!callback)
	return -1;
    frontier_state *state = state_get(gedp, false);
    if (!state)
	return 0;

    int count = 0;
    for (const frontier_root &root : state->roots) {
	if (!root_in_scope(root, view_ctx, mode))
	    continue;
	struct ged_draw_frontier_root_record record;
	record.path = root.path.c_str();
	record.view = root.view;
	record.mode = root.mode;
	record.appearance = &root.appearance;
	count++;
	if (!callback(&record, userdata))
	    break;
    }
    return count;
}


static int
frontier_visibility_emit(
    struct ged *gedp,
    const frontier_backend_change &change,
    ged_draw_frontier_visibility_cb callback,
    void *userdata)
{
    std::vector<std::string> paths;
    std::vector<int> visible;
    if (!change.clear) {
	paths.reserve(change.root.rules.size());
	visible.reserve(change.root.rules.size());
	for (const visibility_rule &rule : change.root.rules) {
	    const std::string path = print_path(gedp, rule.ids);
	    if (path.empty())
		continue;
	    paths.push_back(path);
	    visible.push_back(rule.visible ? 1 : 0);
	}
    }
    std::vector<const char *> path_views;
    path_views.reserve(paths.size());
    for (const std::string &path : paths)
	path_views.push_back(path.c_str());

    struct ged_draw_frontier_visibility_record record;
    record.root_path = change.root.path.c_str();
    record.view = change.root.view;
    record.mode = change.root.mode;
    record.paths = path_views.empty() ? nullptr : path_views.data();
    record.visible = visible.empty() ? nullptr : visible.data();
    record.rule_count = path_views.size();
    record.clear = change.clear || path_views.empty();
    return callback(&record, userdata);
}


extern "C" int
ged_draw_frontier_visibility_changes_foreach(
    struct ged *gedp,
    ged_draw_frontier_visibility_cb callback,
    void *userdata)
{
    if (!gedp || !callback)
	return -1;
    frontier_state *state = state_get(gedp, false);
    if (!state || state->backend_changes.empty())
	return 0;

    std::vector<frontier_backend_change> changes;
    changes.swap(state->backend_changes);
    int count = 0;
    for (size_t i = 0; i < changes.size(); i++) {
	count++;
	if (frontier_visibility_emit(gedp, changes[i], callback, userdata))
	    continue;
	for (; i < changes.size(); i++)
	    frontier_backend_change_note(state, changes[i].root,
		changes[i].clear);
	break;
    }
    return count;
}


extern "C" void
ged_draw_frontier_visibility_changes_discard(struct ged *gedp)
{
    frontier_state *state = state_get(gedp, false);
    if (state)
	state->backend_changes.clear();
}


extern "C" int
ged_draw_frontier_visibility_snapshot_foreach(
    struct ged *gedp,
    ged_draw_frontier_visibility_cb callback,
    void *userdata)
{
    if (!gedp || !callback)
	return -1;
    frontier_state *state = state_get(gedp, false);
    if (!state)
	return 0;

    int count = 0;
    for (const frontier_root &root : state->roots) {
	if (root.rules.empty())
	    continue;
	frontier_backend_change change;
	change.root = root;
	if (!frontier_visibility_emit(gedp, change, callback, userdata))
	    return count;
	count++;
    }
    state->backend_changes.clear();
    return count;
}


extern "C" int
ged_draw_frontier_presentation_snapshot_foreach(
    struct ged *gedp,
    ged_draw_frontier_presentation_cb callback,
    void *userdata)
{
    if (!gedp || !callback)
	return -1;
    frontier_state *state = state_get(gedp, false);
    if (!state)
	return 0;

    int count = 0;
    for (const frontier_root &root : state->roots) {
	for (const presentation_rule &rule : root.presentation) {
	    struct ged_draw_frontier_presentation_record record;
	    record.root_path = root.path.c_str();
	    record.path = rule.path.c_str();
	    record.view = rule.view;
	    record.mode = rule.mode;
	    record.kind = rule.kind;
	    record.match = rule.match;
	    record.value = rule.value;
	    count++;
	    if (!callback(&record, userdata))
		return count;
	}
    }
    return count;
}


extern "C" int
ged_draw_frontier_presentation_set(
    struct ged *gedp,
    const struct ged_scene_reducer_request *txn,
    const char *resolved_path)
{
    if (!gedp || !txn || !resolved_path || !resolved_path[0])
	return 0;
    if (txn->kind != GED_SCENE_REDUCER_VISIBILITY &&
	txn->kind != GED_SCENE_REDUCER_HIGHLIGHT &&
	txn->kind != GED_SCENE_REDUCER_TRANSPARENCY)
	return 0;

    frontier_state *state = state_get(gedp, false);
    if (!state)
	return 0;

    presentation_rule rule;
    rule.kind = txn->kind;
    rule.path = canonical_path(resolved_path);
    rule.view = (txn->view && ged_view_context_is_independent(txn->view)) ?
	txn->view : nullptr;
    rule.mode = txn->mode;
    rule.match = txn->match;
    rule.value = txn->value;
    if (rule.path.empty())
	return 0;

    int changed = 0;
    for (frontier_root &root : state->roots) {
	if (!presentation_root_matches(gedp, root, rule))
	    continue;
	const std::vector<presentation_rule> previous = root.presentation;
	root.presentation.erase(std::remove_if(root.presentation.begin(),
	    root.presentation.end(), [&rule](const presentation_rule &candidate) {
		return candidate.kind == rule.kind &&
		    candidate.path == rule.path &&
		    candidate.view == rule.view &&
		    candidate.mode == rule.mode &&
		    candidate.match == rule.match;
	    }), root.presentation.end());
	root.presentation.push_back(rule);
	if (previous.size() != root.presentation.size() ||
	    !std::equal(previous.begin(), previous.end(),
		root.presentation.begin(), presentation_rule_equal))
	    changed++;
    }
    return changed;
}


extern "C" int
ged_draw_frontier_highlights_clear(struct ged *gedp)
{
    frontier_state *state = state_get(gedp, false);
    if (!state)
	return 0;
    int changed = 0;
    for (frontier_root &root : state->roots) {
	const size_t previous = root.presentation.size();
	root.presentation.erase(std::remove_if(root.presentation.begin(),
	    root.presentation.end(), [](const presentation_rule &rule) {
		return rule.kind == GED_SCENE_REDUCER_HIGHLIGHT;
	    }), root.presentation.end());
	if (previous != root.presentation.size())
	    changed++;
    }
    return changed;
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
    if (!gedp || !path || !path[0])
	return -1;

    const std::string target_path = canonical_path(path);
    std::vector<frontier_root> roots =
	current_roots(gedp, view_ctx, mode);
    path_ids target;
    /* Removal notifications arrive after librt has removed the directory
     * entry.  Preserve exact semantic lookup by retained path long enough for
     * the reducer to retire that draw root; requiring a fresh database-index
     * resolution here left killed top-level roots permanently drawn. */
    if (!resolve_path(gedp, target_path, target)) {
	for (const frontier_root &root : roots) {
	    if (root.path == target_path)
		return 1;
	}
	return -1;
    }

    /* The retained database source is itself the semantic visibility
     * frontier even before it acquires an erase/edit override.  Consulting
     * only state->roots here missed that overwhelmingly common case because
     * state->roots intentionally stores only roots with active rules.  The
     * caller then fell back to enumerating every compact leaf occurrence for
     * every visible Qt tree row.
     *
     * Query cost is proportional to semantic draw roots and path depth, not
     * to realized leaf count. */
    if (roots.empty())
	return -1;

    bool matched = false;
    bool fully_visible = false;
    bool partially_visible = false;
    for (const frontier_root &root : roots) {
	if (!root_in_scope(root, view_ctx, mode) ||
	    (!ids_prefix(root.ids, target) &&
	     !ids_prefix(target, root.ids)))
	    continue;
	matched = true;
	/* One or more draw roots below the queried hierarchy item establish a
	 * partial state.  This is the common Qt-tree case after drawing a child
	 * path directly; it must not require enumerating all of the parent's
	 * leaves to answer the row-paint query. */
	if (ids_prefix(target, root.ids) && target != root.ids) {
	    partially_visible = true;
	    continue;
	}
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
    const struct ged_scene_reducer_request *txn,
    const char *resolved_path)
{
    frontier_state *state = state_get(gedp, false);
    if (!state || !txn)
	return;

    bool structural = false;
    bool all_roots = false;
    switch (txn->kind) {
	case GED_SCENE_REDUCER_DRAW:
	case GED_SCENE_REDUCER_ERASE:
	case GED_SCENE_REDUCER_ERASE_PREFIX:
	case GED_SCENE_REDUCER_SOURCE_RENAMED:
	case GED_SCENE_REDUCER_SOURCE_REFERENCES_REMOVED:
	    structural = true;
	    break;
	case GED_SCENE_REDUCER_SOURCE_UPDATED:
	    structural = txn->removed != 0 || txn->redraw != 0;
	    break;
	case GED_SCENE_REDUCER_CLEAR:
	case GED_SCENE_REDUCER_CLEAR_SCOPE:
	case GED_SCENE_REDUCER_TEARDOWN:
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
    if (txn->kind == GED_SCENE_REDUCER_SOURCE_UPDATED && txn->redraw) {
	for (auto it = state->roots.begin(); it != state->roots.end(); ) {
	    if (!root_in_scope(*it, txn->view, txn->mode)) {
		++it;
		continue;
	    }
	    bool affected = paths.empty();
	    for (const std::string &path : paths) {
		if (path_prefix(it->path, path) ||
		    path_prefix(path, it->path)) {
		    affected = true;
		    break;
		}
	    }
	    if (!affected) {
		++it;
		continue;
	    }
	    if (!frontier_reconcile_database_paths(gedp, *it)) {
		frontier_backend_change_note(state, *it, true);
		it = state->roots.erase(it);
		continue;
	    }
	    frontier_backend_change_note(state, *it, it->rules.empty());
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
	if (!retire && (txn->kind == GED_SCENE_REDUCER_ERASE ||
		txn->kind == GED_SCENE_REDUCER_ERASE_PREFIX ||
		txn->kind == GED_SCENE_REDUCER_DRAW ||
		(txn->kind == GED_SCENE_REDUCER_SOURCE_UPDATED &&
		 txn->removed))) {
	    for (const std::string &path : paths) {
		/* Erasing an owner, or explicitly drawing a broader owner, retires
		 * the narrower source's presentation rules.  A descendant draw is
		 * absorbed by ged_draw_frontier_absorb_draw and keeps the owner. */
		const bool erase_owner = txn->kind != GED_SCENE_REDUCER_DRAW &&
		    path == it->path;
		const bool broader_owner = path != it->path &&
		    path_prefix(path, it->path);
		if (erase_owner || broader_owner) {
		    retire = true;
		    break;
		}
	    }
	}
	if (!retire) {
	    ++it;
	    continue;
	}
	frontier_backend_change_note(state, *it, true);
	it = state->roots.erase(it);
    }
}


extern "C" int
ged_draw_frontier_absorb_draw(
    struct ged *gedp,
    const struct ged_scene_reducer_request *txn,
    const char *resolved_path,
    struct ged_scene_reducer_result *result)
{
    if (!gedp || !txn || txn->kind != GED_SCENE_REDUCER_DRAW ||
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
    std::vector<bool> backend_refresh(roots.size(), false);
    std::vector<bool> added(roots.size(), false);
    frontier_state *state = state_get(gedp, true);
    if (!state)
	return -1;
    int new_roots = 0;
    for (const std::string &path : requested_paths) {
	path_ids target;
	if (!resolve_path(gedp, path, target))
	    return 0;

	size_t owner = roots.size();
	for (size_t i = 0; i < roots.size(); i++) {
	    if (roots[i].mode != mode || !ids_prefix(roots[i].ids, target))
		continue;
	    if (owner == roots.size() ||
		roots[i].ids.size() > roots[owner].ids.size())
		owner = i;
	}
	if (owner != roots.size() && roots[owner].ids == target) {
	    if (!frontier_appearance_equal(roots[owner].appearance,
		    appearance)) {
		roots[owner].appearance = appearance;
		changed[owner] = true;
		backend_refresh[owner] = true;
	    }
	    continue;
	}
	if (owner != roots.size() &&
	    !frontier_appearance_equal(roots[owner].appearance, appearance))
	    owner = roots.size();
	if (owner == roots.size()) {
	    frontier_root root;
	    root.path = path;
	    root.ids = target;
	    root.view = (txn->view &&
		ged_view_context_is_independent(txn->view)) ? txn->view : nullptr;
	    root.mode = mode;
	    root.appearance = appearance;
	    roots.push_back(root);
	    changed.push_back(false);
	    backend_refresh.push_back(false);
	    added.push_back(true);
	    new_roots++;
	    continue;
	}

	if (frontier_set_subtree_visibility(roots[owner], target, true))
	    changed[owner] = true;
    }

    int changed_roots = 0;
    int refresh_roots = 0;
    for (size_t i = 0; i < roots.size(); i++) {
	if (added[i]) {
	    state->roots.push_back(roots[i]);
	    continue;
	}
	if (!changed[i])
	    continue;
	if (retained_root_apply(gedp, roots[i]) < 0)
	    return -1;
	changed_roots++;
	if (backend_refresh[i])
	    refresh_roots++;
    }
    if (result && changed_roots && !new_roots && !refresh_roots) {
	result->affected_groups +=
	    static_cast<int>(requested_paths.size());
	result->affected_shapes +=
	    static_cast<int>(requested_paths.size());
	result->presentation_only = 1;
    }
    return (new_roots || refresh_roots) ? 0 :
	static_cast<int>(requested_paths.size());
}


extern "C" void
ged_scene_edit_internal_result_init(
    struct ged_scene_edit_internal_result *result)
{
    if (!result)
	return;
    *result = {};
    BU_VLS_INIT(&result->paths);
    BU_VLS_INIT(&result->errors);
    result->scope = GED_SCENE_EDIT_SCOPE_REF_NULL;
}


extern "C" void
ged_scene_edit_internal_result_free(
    struct ged_scene_edit_internal_result *result)
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
ged_scene_edit_acquire_internal(
    struct ged *gedp,
    const struct ged_scene_edit_request *request,
    ged_scene_edit_scope_ref *scope,
    struct ged_scene_edit_internal_result *result)
{
    result_prepare(result, gedp);
    if (scope)
	*scope = GED_SCENE_EDIT_SCOPE_REF_NULL;
    if (!gedp || !request || !scope || !ged_db_index_available(gedp)) {
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

    const int requested_mode =
	request->draw_mode == GED_SCENE_DRAW_DEFAULT ? -1 :
	static_cast<int>(request->draw_mode);
    std::vector<path_ids> occurrences =
	promotion_occurrences(gedp, requested_ids, request->occurrences);
    if (occurrences.empty())
	occurrences.push_back(requested_ids);

    std::vector<frontier_root> roots =
	current_roots(gedp, request->view, requested_mode);
    promotion_record promotion;
    promotion.role = request->purpose ? request->purpose : "";

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
	if (!occurrence_drawn[occurrence_index] && request->draw_if_absent &&
	    occurrence == requested_ids) {
	    const std::string occurrence_path = print_path(gedp, occurrence);
	    if (occurrence_path.empty())
		continue;
	    frontier_root root;
	    root.path = occurrence_path;
	    root.ids = occurrence;
	    root.view = (request->view &&
		ged_view_context_is_independent(request->view)) ?
		request->view : nullptr;
	    root.mode = requested_mode >= 0 ? requested_mode :
		ged_draw_default_mode(gedp);
	    const struct ged_draw_appearance_settings default_appearance =
		GED_DRAW_APPEARANCE_SETTINGS_INIT;
	    root.appearance = default_appearance;
	    root.appearance.draw_mode = root.mode;
	    root.auto_draw = true;
	    const char *draw_path = occurrence_path.c_str();
	    struct ged_scene_draw_request draw_request;
	    ged_scene_draw_request_init(&draw_request);
	    draw_request.view = root.view;
	    draw_request.paths = &draw_path;
	    draw_request.path_count = 1;
	    draw_request.style.draw_mode =
		static_cast<enum ged_scene_draw_mode>(root.mode);
	    if (ged_scene_draw(gedp, &draw_request, nullptr) != GED_SCENE_OK) {
		if (result)
		    bu_vls_printf(&result->errors,
			"failed to draw edit target %s",
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

    ged_scene_edit_scope_ref ref = {
	state->owner, promotion.id, promotion.generation
    };
    *scope = ref;
    (void)ged_selection_present_private(gedp);
    if (result) {
	result->status = 1;
	result->scope = ref;
	result->occurrence_count = promotion.occurrences.size();
	for (const std::string &path : promotion.occurrences)
	    append_result_path(result, path);
	result->scene_revision_after = ged_draw_scene_revision(gedp);
    }
    return 1;
}


extern "C" int
ged_scene_edit_release_internal(
    struct ged *gedp,
    ged_scene_edit_scope_ref ref,
    enum ged_scene_edit_outcome outcome,
    struct ged_scene_edit_internal_result *result)
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
		struct ged_scene_erase_request erase_request;
		ged_scene_erase_request_init(&erase_request);
		erase_request.view = root.view;
		erase_request.path = root.path.c_str();
		erase_request.draw_mode =
		    static_cast<enum ged_scene_draw_mode>(root.mode);
		if (ged_scene_erase(gedp, &erase_request, nullptr) !=
		    GED_SCENE_OK) {
		    if (result)
			bu_vls_printf(&result->errors,
			    "failed to remove automatically drawn edit root %s",
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
    (void)ged_selection_present_private(gedp);
    if (result) {
	result->status = conflicts ? 0 : changed;
	result->scope = ref;
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
			     struct ged_scene_reducer_result *result)
{
    if (!gedp || !path || !path[0] || !ged_db_index_available(gedp))
	return 0;

    const std::string target_path = canonical_path(path);
    path_ids target;
    if (!resolve_path(gedp, target_path, target))
	return 0;

    /* An owning root is semantic state even when no legacy per-leaf records
     * were materialized (the normal headless/progressive case).  Depending
     * on the old record eraser here made an exact erase a no-op, leaving the
     * draw intent alive forever.  Retire matching owners directly and let
     * the one committed ERASE delta remove their backend sources. */
    frontier_state *state = state_get(gedp, false);
    int retired_roots = 0;
    if (state) {
	for (auto it = state->roots.begin(); it != state->roots.end(); ) {
	    const bool owns_exact = it->ids == target;
	    const bool below_prefix = prefix && ids_prefix(target, it->ids);
	    if (!root_in_scope(*it, view_ctx, mode) ||
		(!owns_exact && !below_prefix)) {
		++it;
		continue;
	    }
	    frontier_backend_change_note(state, *it, true);
	    it = state->roots.erase(it);
	    retired_roots++;
	}
    }
    if (retired_roots) {
	if (result) {
	    result->affected_groups += retired_roots;
	    result->affected_shapes += retired_roots;
	}
	return retired_roots;
    }

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
    if (result)
	result->presentation_only = 1;
    return changed;
}
