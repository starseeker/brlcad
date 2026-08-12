/*                         E V E N T . C P P
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file libged/event.cpp
 *
 * GED semantic database mutation/reconciliation events.
 */

#include "common.h"

#include <algorithm>
#include <cstring>
#include <string>
#include <unordered_set>
#include <vector>

#include "ged/draw.h"
#include "ged/db_index.h"
#include "ged/event.h"
#include "ged/selection.h"
#include "rt/db5.h"

#include "./ged_private.h"


struct ged_event_owned {
    enum ged_event_kind kind = GED_EVENT_NONE;
    std::string name;
    std::string new_name;
    std::string parent_name;
    std::string child_name;
    std::string path;
    int redraw = 0;
    int librt = 0;
};


struct ged_event_observer_entry {
    ged_event_observer_token token = 0;
    enum ged_event_observer_phase phase = GED_EVENT_OBSERVER_POST_RECONCILE;
    ged_event_observer_func_t func = nullptr;
    void *client_data = nullptr;
    int removed = 0;
};


struct ged_event_service {
    struct ged *gedp = nullptr;
    struct db_i *callback_dbip = nullptr;
    ged_event_observer_token next_token = 1;
    int suspended = 0;
    int bulk_depth = 0;
    int bulk_dirty = 0;
    int bulk_callbacks_were_enabled = 0;
    int bulk_refresh_pending = 0;
    int batch_depth = 0;
    int dispatch_depth = 0;
    std::vector<ged_event_owned> queued_events;
    std::vector<ged_event_owned> followup_events;
    std::vector<ged_event_observer_entry> observers;
};


struct ged_event_result {
    enum ged_event_status status = GED_EVENT_OK;
    size_t event_count = 0;
    size_t coalesced_count = 0;
    size_t draw_change_count = 0;
    size_t db_index_change_count = 0;
    int selection_changed = 0;
    uint64_t scene_revision_before = 0;
    uint64_t scene_revision_after = 0;
    std::vector<std::string> paths;
    std::string diagnostic;
};


static struct ged_event_service *
ged_event_service_get(struct ged *gedp)
{
    return (gedp && gedp->i) ? gedp->i->ged_event_servicep : nullptr;
}


static const char *
ged_event_string_or_null(const std::string &s)
{
    return s.empty() ? nullptr : s.c_str();
}


static struct ged_event
ged_event_view(const ged_event_owned &owned)
{
    struct ged_event ev;
    ev.kind = owned.kind;
    ev.name = ged_event_string_or_null(owned.name);
    ev.new_name = ged_event_string_or_null(owned.new_name);
    ev.parent_name = ged_event_string_or_null(owned.parent_name);
    ev.child_name = ged_event_string_or_null(owned.child_name);
    ev.path = ged_event_string_or_null(owned.path);
    ev.redraw = owned.redraw;
    return ev;
}


static ged_event_owned
ged_event_own(const struct ged_event &ev)
{
    ged_event_owned owned;
    owned.kind = ev.kind;
    if (ev.name)
	owned.name = ev.name;
    if (ev.new_name)
	owned.new_name = ev.new_name;
    if (ev.parent_name)
	owned.parent_name = ev.parent_name;
    if (ev.child_name)
	owned.child_name = ev.child_name;
    if (ev.path)
	owned.path = ev.path;
    owned.redraw = ev.redraw;
    return owned;
}


static void
ged_event_queue(std::vector<ged_event_owned> &queue,
		const struct ged_event *events,
		size_t event_count,
		int librt = 0)
{
    if (!events || !event_count)
	return;
    queue.reserve(queue.size() + event_count);
    for (size_t i = 0; i < event_count; i++) {
	if (events[i].kind == GED_EVENT_NONE)
	    continue;
	ged_event_owned owned = ged_event_own(events[i]);
	owned.librt = librt ? 1 : 0;
	queue.push_back(owned);
    }
}


static int
ged_event_same(const ged_event_owned &a, const ged_event_owned &b)
{
    return a.kind == b.kind &&
	a.name == b.name &&
	a.new_name == b.new_name &&
	a.parent_name == b.parent_name &&
	a.child_name == b.child_name &&
	a.path == b.path &&
	a.redraw == b.redraw;
}


static int
ged_event_contains_batch_rebuild(const std::vector<ged_event_owned> &events)
{
    for (const ged_event_owned &event : events) {
	if (event.kind == GED_EVENT_BATCH_REBUILD)
	    return 1;
    }
    return 0;
}


static std::unordered_set<std::string>
ged_event_librt_fallback_cover_names(const std::vector<ged_event_owned> &events)
{
    std::unordered_set<std::string> names;
    for (const ged_event_owned &event : events) {
	if (event.librt || event.name.empty())
	    continue;
	if (event.kind == GED_EVENT_ATTRIBUTE_CHANGED ||
		event.kind == GED_EVENT_MATERIAL_CHANGED ||
		event.kind == GED_EVENT_OBJECT_VISIBILITY_CHANGED)
	    names.insert(event.name);
    }
    return names;
}


static std::unordered_set<std::string>
ged_event_added_cover_names(const std::vector<ged_event_owned> &events)
{
    std::unordered_set<std::string> names;
    for (const ged_event_owned &event : events) {
	if (event.name.empty())
	    continue;
	if (event.kind == GED_EVENT_OBJECT_ADDED)
	    names.insert(event.name);
    }
    return names;
}


static std::unordered_set<std::string>
ged_event_rename_cover_names(const std::vector<ged_event_owned> &events)
{
    std::unordered_set<std::string> names;
    for (const ged_event_owned &event : events) {
	if (event.librt || event.new_name.empty())
	    continue;
	if (event.kind == GED_EVENT_OBJECT_RENAMED)
	    names.insert(event.new_name);
    }
    return names;
}


static std::unordered_set<std::string>
ged_event_comb_instance_cover_parent_names(const std::vector<ged_event_owned> &events)
{
    std::unordered_set<std::string> names;
    for (const ged_event_owned &event : events) {
	if (event.librt || event.parent_name.empty())
	    continue;
	if (event.kind == GED_EVENT_COMB_INSTANCE_REMOVED ||
		event.kind == GED_EVENT_OBJECT_REFERENCES_REMOVED)
	    names.insert(event.parent_name);
    }
    return names;
}


static int
ged_event_database_metadata_covers_librt(const std::vector<ged_event_owned> &events)
{
    for (const ged_event_owned &event : events) {
	if (!event.librt &&
		event.kind == GED_EVENT_DATABASE_METADATA_CHANGED)
	    return 1;
    }
    return 0;
}


static int
ged_event_librt_named_structural_fallback_covered(
	const ged_event_owned &event,
	const std::unordered_set<std::string> &names)
{
    return event.librt && !event.name.empty() &&
	(event.kind == GED_EVENT_OBJECT_MODIFIED ||
	 event.kind == GED_EVENT_COMB_TREE_CHANGED) &&
	names.find(event.name) != names.end();
}


static int
ged_event_librt_fallback_covered(const ged_event_owned &event,
				 const std::unordered_set<std::string> &names)
{
    return ged_event_librt_named_structural_fallback_covered(event, names);
}


static int
ged_event_librt_database_metadata_covered(const ged_event_owned &event,
					  int metadata_cover)
{
    return metadata_cover && event.librt &&
	event.kind == GED_EVENT_OBJECT_MODIFIED &&
	event.name == DB5_GLOBAL_OBJECT_NAME;
}


static int
ged_event_librt_add_final_state_covered(const ged_event_owned &event,
					const std::unordered_set<std::string> &names)
{
    return ged_event_librt_named_structural_fallback_covered(event, names);
}


static std::vector<ged_event_owned>
ged_event_coalesce(const std::vector<ged_event_owned> &events)
{
    std::unordered_set<std::string> fallback_cover_names =
	ged_event_librt_fallback_cover_names(events);
    std::unordered_set<std::string> added_cover_names =
	ged_event_added_cover_names(events);
    std::unordered_set<std::string> rename_cover_names =
	ged_event_rename_cover_names(events);
    std::unordered_set<std::string> comb_instance_cover_parent_names =
	ged_event_comb_instance_cover_parent_names(events);
    int metadata_cover = ged_event_database_metadata_covers_librt(events);
    std::vector<ged_event_owned> coalesced;
    for (const ged_event_owned &event : events) {
	if (ged_event_librt_fallback_covered(event, fallback_cover_names))
	    continue;
	if (ged_event_librt_database_metadata_covered(event, metadata_cover))
	    continue;
	if (ged_event_librt_add_final_state_covered(event, added_cover_names))
	    continue;
	if (ged_event_librt_named_structural_fallback_covered(event,
		rename_cover_names))
	    continue;
	if (ged_event_librt_named_structural_fallback_covered(event,
		comb_instance_cover_parent_names))
	    continue;

	bool found = false;
	for (ged_event_owned &cevent : coalesced) {
	    if (ged_event_same(event, cevent)) {
		cevent.librt = cevent.librt && event.librt;
		found = true;
		break;
	    }
	    if ((event.kind == GED_EVENT_OBJECT_MODIFIED ||
		 event.kind == GED_EVENT_COMB_TREE_CHANGED ||
		 event.kind == GED_EVENT_ATTRIBUTE_CHANGED) &&
		    event.kind == cevent.kind &&
		    event.name == cevent.name &&
		    event.new_name == cevent.new_name &&
		    event.parent_name == cevent.parent_name &&
		     event.child_name == cevent.child_name &&
		     event.path == cevent.path) {
		cevent.redraw = cevent.redraw || event.redraw;
		cevent.librt = cevent.librt && event.librt;
		found = true;
		break;
	    }
	}
	if (!found)
	    coalesced.push_back(event);
    }
    return coalesced;
}


static std::vector<struct ged_event>
ged_event_views(const std::vector<ged_event_owned> &events)
{
    std::vector<struct ged_event> views;
    views.reserve(events.size());
    for (const ged_event_owned &event : events)
	views.push_back(ged_event_view(event));
    return views;
}


static void
ged_event_result_note_name(struct ged_event_result *result,
			   const char *name)
{
    if (!result || !name || !name[0])
	return;
    if (std::find(result->paths.begin(), result->paths.end(), name) ==
	    result->paths.end())
	result->paths.emplace_back(name);
}


static void
ged_event_result_note_draw(struct ged_event_result *result,
			   const struct ged_draw_transaction_result *draw_result)
{
    if (!result || !draw_result)
	return;
    for (size_t i = 0; i < draw_result->path_count; i++)
	ged_event_result_note_name(result, draw_result->paths[i]);
    if (bu_vls_strlen(&draw_result->errors)) {
	if (!result->diagnostic.empty())
	    result->diagnostic.push_back('\n');
	result->diagnostic.append(bu_vls_cstr(&draw_result->errors));
    }
}


static void
ged_event_collect_index_affected_paths(struct ged *gedp,
				       const std::string &name,
				       std::vector<std::string> &paths)
{
    if (!gedp || name.empty() || !ged_db_index_available(gedp))
	return;

    size_t path_len = ged_db_index_path_resolve(gedp, name.c_str(), nullptr,
	    0);
    if (!path_len)
	return;

    std::vector<ged_db_index_id> ids(path_len);
    if (ged_db_index_path_resolve(gedp, name.c_str(), ids.data(),
		ids.size()) != path_len || ids.empty())
	return;

    ged_db_index_id object_id = ids.back();
    size_t affected_count = ged_db_index_affected_path_count(gedp, object_id,
	    0);
    for (size_t row = 0; row < affected_count; row++) {
	size_t affected_len = ged_db_index_affected_path_at(gedp, object_id,
		row, nullptr, 0, 0);
	if (!affected_len)
	    continue;

	std::vector<ged_db_index_id> affected_ids(affected_len);
	if (ged_db_index_affected_path_at(gedp, object_id, row,
		    affected_ids.data(), affected_ids.size(), 0) !=
		affected_len)
	    continue;

	struct bu_vls printed = BU_VLS_INIT_ZERO;
	if (ged_db_index_path_print(gedp, &printed, affected_ids.data(),
		    affected_ids.size(), 0, 0) &&
		bu_vls_strlen(&printed)) {
	    std::string path = bu_vls_cstr(&printed);
	    if (std::find(paths.begin(), paths.end(), path) == paths.end())
		paths.push_back(path);
	}
	bu_vls_free(&printed);
    }
}


static void
ged_event_result_prepare(struct ged_event_result *result,
			 struct ged *gedp,
			 size_t event_count,
			 size_t coalesced_event_count)
{
    if (!result)
	return;
    ged_event_result_clear(result);
    result->status = GED_EVENT_OK;
    result->event_count = event_count;
    result->coalesced_count = coalesced_event_count;
    result->scene_revision_before = ged_draw_scene_revision(gedp);
    result->scene_revision_after = result->scene_revision_before;
}


static void
ged_event_prune_observers(struct ged_event_service *state)
{
    if (!state || state->dispatch_depth > 0)
	return;
    state->observers.erase(
	    std::remove_if(state->observers.begin(), state->observers.end(),
		[](const ged_event_observer_entry &entry) {
		    return entry.removed || !entry.func;
		}),
	    state->observers.end());
}


static int
ged_event_observer_count(struct ged_event_service *state)
{
    if (!state)
	return 0;

    int count = 0;
    for (const ged_event_observer_entry &entry : state->observers) {
	if (!entry.removed && entry.func)
	    count++;
    }
    return count;
}


static int
ged_event_service_has_live_consumers(struct ged_event_service *state)
{
    if (!state || !state->gedp)
	return 0;

    if (ged_event_observer_count(state) > 0)
	return 1;

    struct ged *gedp = state->gedp;
    if (ged_db_index_available(gedp))
	return 1;

    if (gedp->ged_refresh_handler != GED_REFRESH_FUNC_NULL)
	return 1;

    if (ged_draw_has_shapes(gedp) || ged_draw_has_groups(gedp))
	return 1;

    if (ged_selection_count(gedp, nullptr) > 0)
	return 1;

    return 0;
}


static void
ged_event_dispatch_phase(struct ged_event_service *state,
			 enum ged_event_observer_phase phase,
			 const std::vector<struct ged_event> &events,
			 const struct ged_event_result *result)
{
    if (!state || events.empty())
	return;

    size_t len = state->observers.size();
    for (size_t i = 0; i < len; i++) {
	ged_event_observer_entry &entry = state->observers[i];
	if (entry.removed || entry.phase != phase || !entry.func)
	    continue;
	entry.func(state->gedp, events.data(), events.size(), result,
		entry.client_data);
    }
}


static int
ged_event_apply_draw_txn(struct ged *gedp,
			 struct ged_draw_transaction *txn,
			 struct ged_event_result *result)
{
    struct ged_draw_transaction_result draw_result;
    ged_draw_transaction_result_init(&draw_result);
    int ret = ged_draw_apply_transaction(gedp, txn, &draw_result);
    ged_event_result_note_draw(result, &draw_result);
    ged_draw_transaction_result_free(&draw_result);
    return ret;
}


static int
ged_event_reconcile_draw(struct ged *gedp,
			 const ged_event_owned &event,
			 struct ged_event_result *result)
{
    int ret = 0;

    if (event.librt)
	return 0;

    switch (event.kind) {
	case GED_EVENT_OBJECT_REMOVED: {
	    if (event.name.empty())
		return 0;
	    struct ged_draw_transaction txn =
		ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED,
			event.name.c_str());
	    txn.removed = 1;
	    ret = ged_event_apply_draw_txn(gedp, &txn, result);
	    break;
	}
	case GED_EVENT_OBJECT_RENAMED: {
	    if (event.name.empty() || event.new_name.empty())
		return 0;
	    struct ged_draw_transaction txn =
		ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_RENAMED,
			event.name.c_str());
	    txn.new_path = event.new_name.c_str();
	    ret = ged_event_apply_draw_txn(gedp, &txn, result);
	    break;
	}
	case GED_EVENT_OBJECT_ADDED:
	case GED_EVENT_OBJECT_MODIFIED:
	case GED_EVENT_COMB_TREE_CHANGED:
	case GED_EVENT_ATTRIBUTE_CHANGED: {
	    if (event.name.empty())
		return 0;
	    std::vector<std::string> affected_paths;
	    ged_event_collect_index_affected_paths(gedp, event.name,
		    affected_paths);
	    struct ged_draw_transaction txn =
		ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED,
		    event.name.c_str());
	    std::vector<const char *> affected_path_views;
	    affected_path_views.reserve(affected_paths.size());
	    for (const std::string &affected_path : affected_paths)
		affected_path_views.push_back(affected_path.c_str());
	    txn.paths = affected_path_views.data();
	    txn.path_count = static_cast<int>(affected_path_views.size());
	    txn.redraw = event.redraw ? 1 : 0;
	    ret = ged_event_apply_draw_txn(gedp, &txn, result);
	    for (const std::string &affected_path : affected_paths)
		ged_event_result_note_name(result, affected_path.c_str());
	    break;
	}
	case GED_EVENT_OBJECT_VISIBILITY_CHANGED: {
	    if (event.name.empty())
		return 0;
	    std::vector<std::string> affected_paths;
	    ged_event_collect_index_affected_paths(gedp, event.name,
		    affected_paths);
	    struct directory *dp = db_lookup(gedp->dbip, event.name.c_str(),
		    LOOKUP_QUIET);
	    int visible = (!dp || !(dp->d_flags & RT_DIR_HIDDEN)) ? 1 : 0;
	    struct ged_draw_transaction txn =
		ged_draw_transaction_make_value(GED_DRAW_TXN_VISIBILITY,
			event.name.c_str(), visible ? 1.0 : 0.0);
	    ret = ged_event_apply_draw_txn(gedp, &txn, result);
	    for (const std::string &affected_path : affected_paths)
		ged_event_result_note_name(result, affected_path.c_str());
	    break;
	}
	case GED_EVENT_COMB_INSTANCE_REMOVED: {
	    if (!event.path.empty()) {
		struct ged_draw_transaction erase_txn =
		    ged_draw_transaction_make(GED_DRAW_TXN_ERASE_PREFIX,
			    event.path.c_str());
		ret += ged_event_apply_draw_txn(gedp, &erase_txn, result);
	    }
	    if (!event.parent_name.empty()) {
		struct ged_draw_transaction parent_txn =
		    ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED,
			    event.parent_name.c_str());
		parent_txn.redraw = 1;
		int parent_ret = ged_event_apply_draw_txn(gedp, &parent_txn,
			result);
		if (parent_ret < 0)
		    ret = parent_ret;
		else if (ret >= 0)
		    ret += parent_ret;
	    }
	    break;
	}
	case GED_EVENT_OBJECT_REFERENCES_REMOVED: {
	    if (event.name.empty())
		return 0;
	    if (!event.parent_name.empty()) {
		std::string removed_path = event.path;
		if (removed_path.empty() && !event.child_name.empty())
		    removed_path = event.parent_name + "/" + event.child_name;
		if (!removed_path.empty()) {
		    struct ged_draw_transaction erase_txn =
			ged_draw_transaction_make(GED_DRAW_TXN_ERASE_PREFIX,
				removed_path.c_str());
		    ret += ged_event_apply_draw_txn(gedp, &erase_txn,
			    result);
		}
		struct ged_draw_transaction parent_txn =
		    ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED,
			    event.parent_name.c_str());
		parent_txn.redraw = 1;
		int parent_ret = ged_event_apply_draw_txn(gedp, &parent_txn,
			result);
		if (parent_ret < 0)
		    ret = parent_ret;
		else if (ret >= 0)
		    ret += parent_ret;
		break;
	    }
	    struct ged_draw_transaction txn =
		ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED,
			event.name.c_str());
	    ret = ged_event_apply_draw_txn(gedp, &txn, result);
	    break;
	}
	case GED_EVENT_MATERIAL_CHANGED: {
	    std::vector<std::string> affected_paths;
	    if (!event.name.empty())
		ged_event_collect_index_affected_paths(gedp, event.name,
		    affected_paths);
	    std::vector<const char *> affected_path_views;
	    affected_path_views.reserve(affected_paths.size());
	    for (const std::string &affected_path : affected_paths)
		affected_path_views.push_back(affected_path.c_str());
	    struct ged_draw_transaction changed =
		ged_draw_transaction_make(GED_DRAW_TXN_MATERIAL_CHANGED,
		    nullptr);
	    changed.paths = affected_path_views.empty() ? nullptr :
		affected_path_views.data();
	    changed.path_count = (int)affected_path_views.size();
	    ret = ged_event_apply_draw_txn(gedp, &changed, result);
	    for (const std::string &affected_path : affected_paths)
		ged_event_result_note_name(result, affected_path.c_str());
	    break;
	}
	case GED_EVENT_BATCH_REBUILD: {
	    struct ged_draw_transaction txn =
		ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, nullptr);
	    ret = ged_event_apply_draw_txn(gedp, &txn, result);
	    break;
	}
	case GED_EVENT_DATABASE_METADATA_CHANGED:
	case GED_EVENT_NONE:
	default:
	    ret = 0;
	    break;
    }

    return ret;
}


static void
ged_event_reconcile_selection(struct ged *gedp,
			      struct ged_event_result *result)
{
    if (!result || !ged_selection_state_available(gedp))
	return;

    int recomputed = ged_selection_recompute(gedp, nullptr);
    int draw_synced = ged_selection_present_private(gedp);
    result->selection_changed = (recomputed || draw_synced) ? 1 : 0;
}


static int
ged_event_events_affect_selection(const std::vector<ged_event_owned> &events)
{
    for (const ged_event_owned &event : events) {
	if (event.kind != GED_EVENT_DATABASE_METADATA_CHANGED)
	    return 1;
    }
    return 0;
}


static int
ged_event_reconcile_db_index(struct ged *gedp, const ged_event_owned &event)
{
    int ret = 0;

    switch (event.kind) {
	case GED_EVENT_OBJECT_ADDED:
	    if (!event.name.empty()) {
		ret += ged_db_index_note_object_name_change(gedp,
			event.name.c_str(), GED_DB_INDEX_OBJECT_ADDED);
	    }
	    break;
	case GED_EVENT_OBJECT_REMOVED:
	    if (!event.name.empty()) {
		ret += ged_db_index_note_object_name_change(gedp,
			event.name.c_str(), GED_DB_INDEX_OBJECT_REMOVED);
	    }
	    break;
	case GED_EVENT_OBJECT_RENAMED:
	    if (!event.name.empty()) {
		ret += ged_db_index_note_object_name_change(gedp,
			event.name.c_str(), GED_DB_INDEX_OBJECT_REMOVED);
	    }
	    if (!event.new_name.empty()) {
		ret += ged_db_index_note_object_name_change(gedp,
			event.new_name.c_str(), GED_DB_INDEX_OBJECT_ADDED);
	    }
	    break;
	case GED_EVENT_OBJECT_MODIFIED:
	case GED_EVENT_OBJECT_VISIBILITY_CHANGED:
	case GED_EVENT_COMB_TREE_CHANGED:
	    if (!event.name.empty()) {
		ret += ged_db_index_note_object_name_change(gedp,
			event.name.c_str(), GED_DB_INDEX_OBJECT_CHANGED);
	    }
	    break;
	case GED_EVENT_ATTRIBUTE_CHANGED:
	case GED_EVENT_MATERIAL_CHANGED:
	    break;
	case GED_EVENT_COMB_INSTANCE_REMOVED:
	    if (!event.parent_name.empty()) {
		ret += ged_db_index_note_object_name_change(gedp,
			event.parent_name.c_str(), GED_DB_INDEX_OBJECT_CHANGED);
	    }
	    break;
	case GED_EVENT_OBJECT_REFERENCES_REMOVED:
	    if (!event.parent_name.empty()) {
		ret += ged_db_index_note_object_name_change(gedp,
			event.parent_name.c_str(), GED_DB_INDEX_OBJECT_CHANGED);
	    }
	    break;
	case GED_EVENT_BATCH_REBUILD:
	    if (ged_db_index_available(gedp) && ged_db_index_refresh(gedp))
		ret += 1;
	    break;
	case GED_EVENT_DATABASE_METADATA_CHANGED:
	case GED_EVENT_NONE:
	default:
	    break;
    }

    return ret;
}


static int ged_event_process_owned(struct ged_event_service *state,
				   const std::vector<ged_event_owned> &events,
				   struct ged_event_result *result);


static enum ged_event_status ged_event_publish_impl(
    struct ged *gedp, const struct ged_event *events, size_t event_count,
    struct ged_event_result *result, int librt);


static void
ged_event_process_followups(struct ged_event_service *state)
{
    if (!state || state->batch_depth > 0 || state->dispatch_depth > 0)
	return;

    int guard = 0;
    while (!state->followup_events.empty() && guard++ < 1024) {
	std::vector<ged_event_owned> followups;
	followups.swap(state->followup_events);
	ged_event_process_owned(state, followups, nullptr);
    }
}


static int
ged_event_process_owned(struct ged_event_service *state,
			const std::vector<ged_event_owned> &events,
			struct ged_event_result *result)
{
    if (!state || events.empty())
	return 0;

    std::vector<ged_event_owned> coalesced = ged_event_coalesce(events);
    std::vector<struct ged_event> views = ged_event_views(coalesced);

    struct ged_event_result local_result;
    if (!result) {
	result = &local_result;
    }

    ged_event_result_prepare(result, state->gedp, events.size(),
	    coalesced.size());

    state->dispatch_depth++;
    ged_event_dispatch_phase(state, GED_EVENT_OBSERVER_INTERNAL, views,
	    result);

    int status = 0;
    int index_status = 0;
    std::unordered_set<std::string> material_names;
    for (const ged_event_owned &event : coalesced) {
	if (event.kind == GED_EVENT_MATERIAL_CHANGED && !event.name.empty())
	    material_names.insert(event.name);
    }
    for (const ged_event_owned &event : coalesced) {
	index_status += ged_event_reconcile_db_index(state->gedp, event);

	/* A paired material event refreshes retained path metadata as well as
	 * effective color.  Do not first treat the attribute notification as a
	 * geometry mutation and re-expand the same combination. */
	const int covered_material_attribute =
	    event.kind == GED_EVENT_ATTRIBUTE_CHANGED &&
	    material_names.find(event.name) != material_names.end();
	int ret = covered_material_attribute ? 0 :
	    ged_event_reconcile_draw(state->gedp, event, result);
	if (ret < 0) {
	    status = ret;
	    if (!result->diagnostic.empty())
		result->diagnostic.push_back('\n');
	    result->diagnostic.append("event draw reconciliation failed");
	    if (!event.name.empty())
		result->diagnostic.append(": ").append(event.name);
	} else if (status >= 0) {
	    status += ret;
	}

	ged_event_result_note_name(result, event.name.c_str());
	ged_event_result_note_name(result, event.path.c_str());
	ged_event_result_note_name(result, event.new_name.c_str());
	ged_event_result_note_name(result, event.parent_name.c_str());
	ged_event_result_note_name(result, event.child_name.c_str());
    }

    result->db_index_change_count = index_status > 0 ?
	static_cast<size_t>(index_status) : 0;
    result->draw_change_count = status > 0 ?
	static_cast<size_t>(status) : 0;
    if (ged_event_events_affect_selection(coalesced))
	ged_event_reconcile_selection(state->gedp, result);
    result->status = status < 0 ? GED_EVENT_ERROR : GED_EVENT_OK;
    result->scene_revision_after = ged_draw_scene_revision(state->gedp);

    ged_event_dispatch_phase(state, GED_EVENT_OBSERVER_POST_RECONCILE, views,
	    result);
    state->dispatch_depth--;
    ged_event_prune_observers(state);

    if (state->bulk_refresh_pending &&
	    ged_event_contains_batch_rebuild(coalesced)) {
	state->bulk_refresh_pending = 0;
	if (state->gedp &&
		state->gedp->ged_refresh_handler != GED_REFRESH_FUNC_NULL)
	    ged_refresh_cb(state->gedp);
    }

    ged_event_process_followups(state);
    return result->status;
}


static void
ged_event_librt_changed_cb(struct db_i *UNUSED(dbip),
			   struct directory *dp,
			   int mode,
			   void *u_data)
{
    struct ged *gedp = (struct ged *)u_data;
    if (!gedp || !dp || !dp->d_namep)
	return;

    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.name = dp->d_namep;

    switch (mode) {
	case 0:
	    ev.kind = (dp->d_flags & RT_DIR_COMB) ?
		GED_EVENT_COMB_TREE_CHANGED : GED_EVENT_OBJECT_MODIFIED;
	    ev.redraw = 0;
	    ged_event_publish_impl(gedp, &ev, 1, nullptr, 1);
	    break;
	case 1:
	    ev.kind = GED_EVENT_OBJECT_ADDED;
	    ged_event_publish_impl(gedp, &ev, 1, nullptr, 1);
	    break;
	case 2:
	    ev.kind = GED_EVENT_OBJECT_REMOVED;
	    ged_event_publish_impl(gedp, &ev, 1, nullptr, 1);
	    break;
	default:
	    break;
    }
}


struct ged_event_service *
ged_event_service_create(struct ged *gedp)
{
    struct ged_event_service *state = new ged_event_service;
    state->gedp = gedp;
    return state;
}


void
ged_event_service_destroy(struct ged_event_service *state)
{
    if (!state)
	return;
    if (state->callback_dbip)
	db_rm_changed_clbk(state->callback_dbip, ged_event_librt_changed_cb,
	    (void *)state->gedp);
    delete state;
}


int
ged_event_available(struct ged *gedp)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    return (state && !state->suspended) ? 1 : 0;
}


enum ged_event_status
ged_event_disable(struct ged *gedp)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state)
	return GED_EVENT_UNAVAILABLE;

    state->suspended++;
    state->bulk_depth = 0;
    state->bulk_dirty = 0;
    state->bulk_callbacks_were_enabled = 0;
    state->bulk_refresh_pending = 0;
    state->batch_depth = 0;
    state->queued_events.clear();
    state->followup_events.clear();
    if (state->callback_dbip)
	ged_event_librt_callbacks_disable(gedp);
    return GED_EVENT_OK;
}


enum ged_event_status
ged_event_bulk_begin(struct ged *gedp)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state || state->suspended)
	return GED_EVENT_UNAVAILABLE;

    if (state->bulk_depth == 0) {
	state->bulk_callbacks_were_enabled = (state->callback_dbip != nullptr);
	if (state->callback_dbip)
	    ged_event_librt_callbacks_disable(gedp);
	state->bulk_dirty = 1;
	state->bulk_refresh_pending = 0;
    }

    state->bulk_depth++;
    return GED_EVENT_OK;
}


enum ged_event_status
ged_event_bulk_end(struct ged *gedp, struct ged_event_result *result)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state || state->bulk_depth <= 0)
	return state ? GED_EVENT_INVALID : GED_EVENT_UNAVAILABLE;

    state->bulk_depth--;
    if (state->bulk_depth > 0)
	return GED_EVENT_OK;

    int dirty = state->bulk_dirty;
    int live = ged_event_service_has_live_consumers(state);
    int restore_callbacks = state->bulk_callbacks_were_enabled;

    state->bulk_dirty = 0;
    state->bulk_callbacks_were_enabled = 0;

    if (restore_callbacks && !state->suspended && gedp && gedp->dbip)
	ged_event_librt_callbacks_enable(gedp);

    if (!dirty || !live)
	return GED_EVENT_OK;

    if (gedp && gedp->ged_refresh_handler != GED_REFRESH_FUNC_NULL)
	state->bulk_refresh_pending = 1;

    return ged_event_notify_batch_rebuild(gedp, result);
}


int
ged_event_bulk_active(struct ged *gedp)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    return (state && !state->suspended && state->bulk_depth > 0) ? 1 : 0;
}


int
ged_event_has_live_consumers(struct ged *gedp)
{
    return ged_event_service_has_live_consumers(ged_event_service_get(gedp));
}


enum ged_event_status
ged_event_enable(struct ged *gedp)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state || state->suspended <= 0)
	return state ? GED_EVENT_INVALID : GED_EVENT_UNAVAILABLE;

    state->suspended--;
    if (!state->suspended && gedp->dbip)
	ged_event_librt_callbacks_enable(gedp);
    return GED_EVENT_OK;
}


enum ged_event_status
ged_event_librt_callbacks_enable(struct ged *gedp)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state || !gedp || !gedp->dbip || state->suspended ||
	    state->bulk_depth > 0)
	return GED_EVENT_UNAVAILABLE;

    if (state->callback_dbip == gedp->dbip)
	return GED_EVENT_OK;

    if (state->callback_dbip)
	db_rm_changed_clbk(state->callback_dbip, ged_event_librt_changed_cb,
	    (void *)gedp);

    if (db_add_changed_clbk(gedp->dbip, ged_event_librt_changed_cb,
	    (void *)gedp) != 0) {
	state->callback_dbip = nullptr;
	return GED_EVENT_ERROR;
    }

    state->callback_dbip = gedp->dbip;
    return GED_EVENT_OK;
}


enum ged_event_status
ged_event_librt_callbacks_disable(struct ged *gedp)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state)
	return GED_EVENT_UNAVAILABLE;
    if (!state->callback_dbip)
	return GED_EVENT_OK;

    int ret = db_rm_changed_clbk(state->callback_dbip,
	ged_event_librt_changed_cb, (void *)gedp);
    state->callback_dbip = nullptr;
    return ret > 0 ? GED_EVENT_OK : GED_EVENT_ERROR;
}


void
ged_event_init(struct ged_event *event)
{
    if (!event)
	return;
    std::memset(event, 0, sizeof(*event));
}


struct ged_event_result *
ged_event_result_create(void)
{
    return new ged_event_result;
}


void
ged_event_result_clear(struct ged_event_result *result)
{
    if (!result)
	return;
    result->status = GED_EVENT_OK;
    result->event_count = 0;
    result->coalesced_count = 0;
    result->draw_change_count = 0;
    result->db_index_change_count = 0;
    result->selection_changed = 0;
    result->scene_revision_before = 0;
    result->scene_revision_after = 0;
    result->paths.clear();
    result->diagnostic.clear();
}


void
ged_event_result_destroy(struct ged_event_result *result)
{
    delete result;
}


enum ged_event_status
ged_event_result_status(const struct ged_event_result *result)
{
    return result ? result->status : GED_EVENT_INVALID;
}


size_t
ged_event_result_event_count(const struct ged_event_result *result)
{
    return result ? result->event_count : 0;
}


size_t
ged_event_result_coalesced_count(const struct ged_event_result *result)
{
    return result ? result->coalesced_count : 0;
}


size_t
ged_event_result_draw_change_count(const struct ged_event_result *result)
{
    return result ? result->draw_change_count : 0;
}


size_t
ged_event_result_db_index_change_count(
    const struct ged_event_result *result)
{
    return result ? result->db_index_change_count : 0;
}


int
ged_event_result_selection_changed(const struct ged_event_result *result)
{
    return result ? result->selection_changed : 0;
}


uint64_t
ged_event_result_scene_revision_before(const struct ged_event_result *result)
{
    return result ? result->scene_revision_before : 0;
}


uint64_t
ged_event_result_scene_revision_after(const struct ged_event_result *result)
{
    return result ? result->scene_revision_after : 0;
}


size_t
ged_event_result_path_count(const struct ged_event_result *result)
{
    return result ? result->paths.size() : 0;
}


const char *
ged_event_result_path_at(const struct ged_event_result *result, size_t index)
{
    return (result && index < result->paths.size()) ?
	result->paths[index].c_str() : nullptr;
}


const char *
ged_event_result_diagnostic(const struct ged_event_result *result)
{
    return result ? result->diagnostic.c_str() : "";
}


enum ged_event_status
ged_event_batch_begin(struct ged *gedp)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state || state->suspended)
	return GED_EVENT_UNAVAILABLE;
    state->batch_depth++;
    return GED_EVENT_OK;
}


enum ged_event_status
ged_event_batch_end(struct ged *gedp, struct ged_event_result *result)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state || state->batch_depth <= 0)
	return state ? GED_EVENT_INVALID : GED_EVENT_UNAVAILABLE;

    state->batch_depth--;
    if (state->batch_depth > 0)
	return GED_EVENT_OK;

    if (state->bulk_depth > 0) {
	if (!state->queued_events.empty()) {
	    state->bulk_dirty = 1;
	    state->queued_events.clear();
	}
	state->followup_events.clear();
	return GED_EVENT_OK;
    }

    std::vector<ged_event_owned> events;
    events.swap(state->queued_events);
    enum ged_event_status ret = GED_EVENT_OK;
    if (!events.empty())
	ret = static_cast<enum ged_event_status>(
	    ged_event_process_owned(state, events, result));
    else if (result)
	ged_event_result_clear(result);
    ged_event_process_followups(state);
    return ret;
}


static enum ged_event_status
ged_event_publish_impl(struct ged *gedp,
		       const struct ged_event *events,
		       size_t event_count,
		       struct ged_event_result *result,
		       int librt)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state || state->suspended)
	return GED_EVENT_UNAVAILABLE;
    if (!events || !event_count)
	return GED_EVENT_INVALID;

    if (state->bulk_depth > 0) {
	state->bulk_dirty = 1;
	return GED_EVENT_OK;
    }

    if (state->dispatch_depth > 0) {
	ged_event_queue(state->followup_events, events, event_count, librt);
	return GED_EVENT_OK;
    }

    if (state->batch_depth > 0) {
	ged_event_queue(state->queued_events, events, event_count, librt);
	return GED_EVENT_OK;
    }

    std::vector<ged_event_owned> owned;
    ged_event_queue(owned, events, event_count, librt);
    enum ged_event_status ret = static_cast<enum ged_event_status>(
	ged_event_process_owned(state, owned, result));
    ged_event_process_followups(state);
    return ret;
}


enum ged_event_status
ged_event_publish(struct ged *gedp,
		  const struct ged_event *events,
		  size_t event_count,
		  struct ged_event_result *result)
{
    return ged_event_publish_impl(gedp, events, event_count, result, 0);
}


enum ged_event_status
ged_event_notify_object_added(struct ged *gedp,
			      const char *name,
			      struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_OBJECT_ADDED;
    ev.name = name;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_object_removed(struct ged *gedp,
				const char *name,
				struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_OBJECT_REMOVED;
    ev.name = name;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_object_renamed(struct ged *gedp,
				const char *old_name,
				const char *new_name,
				struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_OBJECT_RENAMED;
    ev.name = old_name;
    ev.new_name = new_name;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_object_modified(struct ged *gedp,
				 const char *name,
				 int redraw,
				 struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_OBJECT_MODIFIED;
    ev.name = name;
    ev.redraw = redraw ? 1 : 0;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_object_visibility_changed(struct ged *gedp,
					   const char *name,
					   struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_OBJECT_VISIBILITY_CHANGED;
    ev.name = name;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_comb_tree_changed(struct ged *gedp,
				   const char *name,
				   int redraw,
				   struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_COMB_TREE_CHANGED;
    ev.name = name;
    ev.redraw = redraw ? 1 : 0;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_comb_instance_removed(struct ged *gedp,
				       const char *parent_name,
				       const char *child_name,
				       const char *path,
				       struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_COMB_INSTANCE_REMOVED;
    ev.parent_name = parent_name;
    ev.child_name = child_name;
    ev.path = path;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_object_references_removed(struct ged *gedp,
					   const char *name,
					   struct ged_event_result *result)
{
    return ged_event_notify_object_reference_removed_from_parent(gedp, name,
	    nullptr, nullptr, result);
}


enum ged_event_status
ged_event_notify_object_reference_removed_from_parent(struct ged *gedp,
						      const char *name,
						      const char *parent_name,
						      const char *path,
						      struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_OBJECT_REFERENCES_REMOVED;
    ev.name = name;
    ev.parent_name = parent_name;
    ev.child_name = name;
    ev.path = path;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_attribute_changed(struct ged *gedp,
				   const char *name,
				   int redraw,
				   struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_ATTRIBUTE_CHANGED;
    ev.name = name;
    ev.redraw = redraw ? 1 : 0;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_material_changed(struct ged *gedp,
				  struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_MATERIAL_CHANGED;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_object_material_changed(struct ged *gedp,
					 const char *name,
					 struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_MATERIAL_CHANGED;
    ev.name = name;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_batch_rebuild(struct ged *gedp,
			       struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_BATCH_REBUILD;
    return ged_event_publish(gedp, &ev, 1, result);
}


enum ged_event_status
ged_event_notify_database_metadata_changed(struct ged *gedp,
					   struct ged_event_result *result)
{
    struct ged_event ev;
    std::memset(&ev, 0, sizeof(ev));
    ev.kind = GED_EVENT_DATABASE_METADATA_CHANGED;
    return ged_event_publish(gedp, &ev, 1, result);
}


ged_event_observer_token
ged_event_observer_add(struct ged *gedp,
		       enum ged_event_observer_phase phase,
		       ged_event_observer_func_t func,
		       void *client_data)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state || !func)
	return 0;
    if (phase != GED_EVENT_OBSERVER_INTERNAL &&
	phase != GED_EVENT_OBSERVER_POST_RECONCILE)
	return 0;

    ged_event_observer_entry entry;
    entry.token = state->next_token++;
    if (!entry.token)
	entry.token = state->next_token++;
    entry.phase = phase;
    entry.func = func;
    entry.client_data = client_data;
    state->observers.push_back(entry);
    return entry.token;
}


int
ged_event_observer_remove(struct ged *gedp, ged_event_observer_token token)
{
    struct ged_event_service *state = ged_event_service_get(gedp);
    if (!state || !token)
	return 0;

    for (ged_event_observer_entry &entry : state->observers) {
	if (entry.token != token || entry.removed)
	    continue;
	entry.removed = 1;
	if (state->dispatch_depth == 0)
	    ged_event_prune_observers(state);
	return 1;
    }

    return 0;
}


/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
