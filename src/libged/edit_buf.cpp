/*                    E D I T _ B U F . C P P
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
/** @file edit_buf.cpp
 *
 * Temporary edit buffer API.
 *
 * Maintains a per-gedp map of (path_string -> rt_edit *) that persists across
 * command invocations within a session.  Allows multi-step CLI edits without
 * a disk write per step.
 *
 * The buffer is stored inside Ged_Internal (ged_private.h) and is destroyed
 * (entries abandoned) when the Ged_Internal destructor runs.
 */

#include "common.h"

#include <string>
#include <vector>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "raytrace.h"
#include "rt/db_fullpath.h"
#include "rt/edit.h"
#include "rt/func.h"

#include "ged/db_index.h"
#include "ged/edit.h"
#include "ged/event.h"
#include "ged/scene.h"
#include "ged/view_edit.h"
#include "ged/view_types.h"
#include "./ged_private.h"
#include "./ged_scene_backend_private.h"


/* ------------------------------------------------------------------ *
 * Ged_Internal destructor — abandons all pending edits on close.
 * ------------------------------------------------------------------ */

Ged_Internal::~Ged_Internal()
{
    for (auto &kv : edit_buf) {
	/* Normal GED teardown closes sessions while the semantic scene is still
	 * alive.  A partially initialized or exceptional owner may reach this
	 * fallback after draw teardown, so only release owned memory here. */
	rt_edit_destroy(kv.second.s);
	db_free_full_path(&kv.second.dfp);
    }
    edit_buf.clear();
    /* Normal ged_free teardown detaches earlier, while view and scene state
     * are intact.  Keep this idempotent fallback for partially initialized
     * owners and exceptional destruction paths. */
    ged_scene_backend_detach_private(gedp);
    ged_draw_frontier_state_destroy(gedp);
}


/* ------------------------------------------------------------------ *
 * Internal helper: compute the canonical path key from a db_full_path.
 * Uses the full path string so there is zero collision risk.
 * ------------------------------------------------------------------ */

static std::string
_edit_buf_key(const struct db_full_path *dfp)
{
    char *pstr = db_path_to_string(dfp);
    std::string key(pstr);
    bu_free(pstr, "edit_buf path key");
    return key;
}


static int
_edit_resolve_key(std::string &key, struct ged *gedp, const char *path,
	struct db_full_path *resolved)
{
    if (!gedp || !gedp->dbip || !path || !path[0])
	return 0;
    struct db_full_path local_path;
    struct db_full_path *dfp = resolved ? resolved : &local_path;
    db_full_path_init(dfp);
    if (db_string_to_path(dfp, gedp->dbip, path) < 0) {
	db_free_full_path(dfp);
	return 0;
    }
    key = _edit_buf_key(dfp);
    if (!resolved)
	db_free_full_path(dfp);
    return 1;
}


/* A primitive edit ultimately commits the terminal database object, not one
 * occurrence of it.  Distinct occurrence paths may need distinct presentation
 * matrices, so they cannot safely own independent writable copies of the same
 * leaf.  Until a multi-occurrence presentation attachment is introduced, make
 * that conflict explicit rather than permitting last-commit-wins data loss. */
static bool
_edit_same_terminal(const Ged_Internal::ged_edit_buf_entry &entry,
	const struct db_full_path *dfp)
{
    if (!dfp || entry.path_names.empty())
	return false;
    struct directory *terminal = DB_FULL_PATH_CUR_DIR(dfp);
    return terminal && terminal->d_namep &&
	entry.path_names.back() == terminal->d_namep;
}


static std::atomic<uint64_t> edit_owner_seed(1);


static uint64_t
_edit_owner(Ged_Internal *gi)
{
    if (!gi)
	return 0;
    if (!gi->edit_owner) {
	gi->edit_owner = edit_owner_seed.fetch_add(1);
	if (!gi->edit_owner)
	    gi->edit_owner = edit_owner_seed.fetch_add(1);
    }
    return gi->edit_owner;
}


static ged_edit_session_ref
_edit_ref(Ged_Internal *gi, const Ged_Internal::ged_edit_buf_entry &entry)
{
    ged_edit_session_ref ref = {
	_edit_owner(gi), entry.id, entry.generation
    };
    return ref;
}


/* Publish the client-independent edit wireframe.  GUI plugins may add richer
 * labels and manipulators, but every GED host must show the same authoritative
 * intermediate geometry without installing a toolkit-specific observer. */
static void
_edit_present(struct ged *gedp, enum ged_edit_session_event_kind kind,
	const std::string &path,
	const Ged_Internal::ged_edit_buf_entry &entry)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;

    char feature_name[96] = {0};
    const ged_edit_session_ref ref = _edit_ref(gedp->i->i, entry);
    snprintf(feature_name, sizeof(feature_name),
	"_ged_edit_preview_%016llx_%016llx",
	(unsigned long long)ref.owner, (unsigned long long)ref.id);

    struct ged_view_edit_transaction transaction =
	ged_view_edit_transaction_default();
    transaction.feature_name = feature_name;
    transaction.source_path = path.c_str();

    switch (kind) {
	case GED_EDIT_SESSION_BEGIN:
	    transaction.event = GED_VIEW_EDIT_PREVIEW_BEGIN;
	    break;
	case GED_EDIT_SESSION_UPDATE:
	case GED_EDIT_SESSION_REVERT:
	    transaction.event = GED_VIEW_EDIT_PREVIEW_UPDATE;
	    break;
	case GED_EDIT_SESSION_COMMIT:
	    transaction.event = GED_VIEW_EDIT_PREVIEW_COMMIT;
	    (void)ged_view_edit_transaction_apply_all(gedp, &transaction);
	    return;
	case GED_EDIT_SESSION_CANCEL:
	    transaction.event = GED_VIEW_EDIT_PREVIEW_CANCEL;
	    (void)ged_view_edit_transaction_apply_all(gedp, &transaction);
	    return;
	case GED_EDIT_SESSION_INVALIDATE:
	    transaction.event = GED_VIEW_EDIT_PREVIEW_DISCARD;
	    (void)ged_view_edit_transaction_apply_all(gedp, &transaction);
	    return;
	case GED_EDIT_SESSION_CHECKPOINT:
	default:
	    return;
    }

    if (!entry.s)
	return;

    mat_t path_matrix;
    MAT_IDN(path_matrix);
    const fastf_t *matrix = NULL;
    if (entry.dfp.fp_len && gedp->dbip &&
	    db_path_to_mat(gedp->dbip,
		(struct db_full_path *)&entry.dfp, path_matrix,
		(int)entry.dfp.fp_len - 1))
	matrix = path_matrix;

    struct rt_wdb *wdbp = gedp->dbip ?
	wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT) : NULL;
    struct bn_tol tolerance = BN_TOL_INIT_TOL;
    transaction.dbip = gedp->dbip;
    transaction.internal = &entry.s->es_int;
    transaction.matrix = matrix;
    transaction.ttol = wdbp ? &wdbp->wdb_ttol : NULL;
    transaction.tol = &tolerance;
    transaction.color_valid = 1;
    transaction.color[0] = 255;
    transaction.color[1] = 255;
    transaction.color[2] = 255;
    transaction.source_revision = (uint32_t)entry.revision;
    transaction.inputs_revision = (uint32_t)entry.revision;
    (void)ged_view_edit_transaction_apply_all(gedp, &transaction);
}


static void
_edit_notify(struct ged *gedp, enum ged_edit_session_event_kind kind,
	const std::string &path,
	const Ged_Internal::ged_edit_buf_entry &entry,
	const std::string *replacement_path = nullptr,
	enum ged_edit_session_invalidation_reason invalidation_reason =
	    GED_EDIT_INVALIDATION_NONE)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;

    Ged_Internal *gi = gedp->i->i;

    _edit_present(gedp, kind, path, entry);

    Ged_Internal::ged_edit_pending_event pending;
    pending.kind = kind;
    pending.session = _edit_ref(gi, entry);
    pending.path = path;
    if (replacement_path)
	pending.replacement_path = *replacement_path;
    pending.invalidation_reason = invalidation_reason;
    pending.revision = entry.revision;
    pending.command_id = entry.last_cmd;
    gi->edit_pending_events.push_back(pending);

    if (gi->edit_dispatch_depth)
	return;

    gi->edit_dispatch_depth++;
    while (!gi->edit_pending_events.empty()) {
	Ged_Internal::ged_edit_pending_event current =
	    gi->edit_pending_events.front();
	gi->edit_pending_events.erase(gi->edit_pending_events.begin());

	std::vector<ged_edit_observer_token> observers;
	observers.reserve(gi->edit_observers.size());
	for (const auto &observer : gi->edit_observers)
	    observers.push_back(observer.first);

	struct ged_edit_session_event event;
	event.kind = current.kind;
	event.session = current.session;
	event.path = current.path.c_str();
	event.replacement_path = current.replacement_path.empty() ? nullptr :
	    current.replacement_path.c_str();
	event.invalidation_reason = current.invalidation_reason;
	event.revision = current.revision;
	event.command_id = current.command_id;

	for (ged_edit_observer_token token : observers) {
	    auto found = gi->edit_observers.find(token);
	    if (found == gi->edit_observers.end() || !found->second.first)
		continue;
	    found->second.first(gedp, &event, found->second.second);
	}
    }
    gi->edit_dispatch_depth--;
}


static std::unordered_map<std::string,
    Ged_Internal::ged_edit_buf_entry>::iterator
_edit_find_ref(Ged_Internal *gi, ged_edit_session_ref ref)
{
    if (!gi || ref.owner != _edit_owner(gi) || !ref.id || !ref.generation)
	return gi ? gi->edit_buf.end() :
	    std::unordered_map<std::string,
		Ged_Internal::ged_edit_buf_entry>::iterator();
    for (auto it = gi->edit_buf.begin(); it != gi->edit_buf.end(); ++it) {
	if (it->second.id == ref.id &&
	    it->second.generation == ref.generation)
	    return it;
    }
    return gi->edit_buf.end();
}


static void
_edit_buf_release_scope(struct ged *gedp,
			ged_scene_edit_scope_ref scope,
			enum ged_scene_edit_outcome outcome)
{
    if (!gedp || ged_scene_edit_scope_ref_is_null(scope))
	return;
    (void)ged_scene_edit_release(gedp, scope, outcome, NULL);
}


static bool
_edit_entry_has_name(const Ged_Internal::ged_edit_buf_entry &entry,
	const char *name, size_t first = 0)
{
    if (!name || !name[0])
	return false;
    for (size_t i = first; i < entry.path_names.size(); ++i) {
	if (entry.path_names[i] == name)
	    return true;
    }
    return false;
}


static std::string
_edit_canonical_event_path(const char *path)
{
    if (!path || !path[0])
	return std::string();
    std::string canonical(path);
    if (canonical[0] != '/')
	canonical.insert(canonical.begin(), '/');
    while (canonical.size() > 1 && canonical.back() == '/')
	canonical.pop_back();
    return canonical;
}


static bool
_edit_path_is_at_or_below(const std::string &session_path,
	const char *changed_path)
{
    const std::string changed = _edit_canonical_event_path(changed_path);
    if (changed.empty() || session_path.size() < changed.size() ||
	session_path.compare(0, changed.size(), changed) != 0)
	return false;
    return session_path.size() == changed.size() ||
	session_path[changed.size()] == '/';
}


static std::string
_edit_renamed_path(const Ged_Internal::ged_edit_buf_entry &entry,
	const char *old_name, const char *new_name)
{
    if (!old_name || !old_name[0] || !new_name || !new_name[0])
	return std::string();
    std::string result;
    for (size_t i = 0; i < entry.path_names.size(); ++i) {
	result.push_back('/');
	result.append(entry.path_names[i] == old_name ? new_name :
	    entry.path_names[i]);
	if (i < entry.path_instances.size() && entry.path_instances[i]) {
	    result.push_back('@');
	    result.append(std::to_string(entry.path_instances[i]));
	}
    }
    return result;
}


static void
_edit_release_key(struct ged *gedp, const std::string &key,
	enum ged_edit_session_event_kind kind,
	const std::string *replacement_path = nullptr,
	enum ged_edit_session_invalidation_reason invalidation_reason =
	    GED_EDIT_INVALIDATION_NONE)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;
    Ged_Internal *gi = gedp->i->i;
    auto it = gi->edit_buf.find(key);
    if (it == gi->edit_buf.end())
	return;

    const ged_scene_edit_scope_ref scope = it->second.edit_scope;
    const Ged_Internal::ged_edit_buf_entry final_entry = it->second;
    rt_edit_destroy(it->second.s);
    db_free_full_path(&it->second.dfp);
    gi->edit_buf.erase(it);
    _edit_buf_release_scope(gedp, scope, GED_SCENE_EDIT_CANCEL);
    _edit_notify(gedp, kind, key, final_entry, replacement_path,
	invalidation_reason);
}


struct ged_edit_invalidation {
    std::string replacement_path;
    enum ged_edit_session_invalidation_reason reason =
	GED_EDIT_INVALIDATION_NONE;
};


static void
_edit_database_event_cb(struct ged *gedp,
	const struct ged_event *events, size_t event_count,
	const struct ged_event_result *UNUSED(result), void *client_data)
{
    Ged_Internal *gi = static_cast<Ged_Internal *>(client_data);
    if (!gedp || !gedp->i || !gi || !events || !event_count ||
	gedp->i->i != gi ||
	gi->edit_buf.empty())
	return;

    std::unordered_map<std::string, ged_edit_invalidation> invalidations;
    for (size_t ei = 0; ei < event_count; ++ei) {
	const struct ged_event &event = events[ei];
	for (const auto &session : gi->edit_buf) {
	    const Ged_Internal::ged_edit_buf_entry &entry = session.second;
	    if (entry.committing)
		continue;

	    bool affected = false;
	    enum ged_edit_session_invalidation_reason reason =
		GED_EDIT_INVALIDATION_NONE;
	    std::string replacement;
	    switch (event.kind) {
		case GED_EVENT_OBJECT_REMOVED:
		    affected = _edit_entry_has_name(entry, event.name);
		    reason = GED_EDIT_INVALIDATION_SOURCE_REMOVED;
		    break;
		case GED_EVENT_OBJECT_RENAMED:
		    affected = _edit_entry_has_name(entry, event.name);
		    reason = GED_EDIT_INVALIDATION_SOURCE_RENAMED;
		    if (affected)
			replacement = _edit_renamed_path(entry, event.name,
			    event.new_name);
		    break;
		case GED_EVENT_OBJECT_MODIFIED:
		case GED_EVENT_COMB_TREE_CHANGED:
		    affected = _edit_entry_has_name(entry, event.name);
		    reason = GED_EDIT_INVALIDATION_SOURCE_CHANGED;
		    break;
		case GED_EVENT_COMB_INSTANCE_REMOVED: {
		    const bool path_removed = _edit_path_is_at_or_below(
			session.first, event.path);
		    affected = path_removed || _edit_entry_has_name(entry,
			event.parent_name);
		    reason = path_removed ? GED_EDIT_INVALIDATION_PATH_REMOVED :
			GED_EDIT_INVALIDATION_SOURCE_CHANGED;
		    break;
		}
		case GED_EVENT_OBJECT_REFERENCES_REMOVED: {
		    const bool path_removed = _edit_path_is_at_or_below(
			session.first, event.path) || (!event.parent_name &&
			_edit_entry_has_name(entry, event.name, 1));
		    affected = path_removed || _edit_entry_has_name(entry,
			event.parent_name);
		    reason = path_removed ? GED_EDIT_INVALIDATION_PATH_REMOVED :
			GED_EDIT_INVALIDATION_SOURCE_CHANGED;
		    break;
		}
		case GED_EVENT_BATCH_REBUILD:
		    affected = true;
		    reason = GED_EDIT_INVALIDATION_DATABASE_REBUILT;
		    break;
		case GED_EVENT_NONE:
		case GED_EVENT_OBJECT_ADDED:
		case GED_EVENT_OBJECT_VISIBILITY_CHANGED:
		case GED_EVENT_ATTRIBUTE_CHANGED:
		case GED_EVENT_MATERIAL_CHANGED:
		case GED_EVENT_DATABASE_METADATA_CHANGED:
		default:
		    break;
	    }
	    if (!affected)
		continue;
	    ged_edit_invalidation &pending = invalidations[session.first];
	    const bool removes_path =
		reason == GED_EDIT_INVALIDATION_SOURCE_REMOVED ||
		reason == GED_EDIT_INVALIDATION_PATH_REMOVED;
	    const bool pending_removes_path =
		pending.reason == GED_EDIT_INVALIDATION_SOURCE_REMOVED ||
		pending.reason == GED_EDIT_INVALIDATION_PATH_REMOVED;
	    if (removes_path) {
		pending.reason = reason;
		pending.replacement_path.clear();
	    } else if (!pending_removes_path &&
		reason == GED_EDIT_INVALIDATION_SOURCE_RENAMED) {
		pending.reason = reason;
		pending.replacement_path = replacement;
	    } else if (pending.reason == GED_EDIT_INVALIDATION_NONE) {
		pending.reason = reason;
	    }
	}
    }

    if (invalidations.empty())
	return;
    const bool already_closing = gi->edit_closing;
    gi->edit_closing = true;
    for (const auto &pending : invalidations) {
	const std::string *replacement = pending.second.replacement_path.empty() ?
	    nullptr : &pending.second.replacement_path;
	_edit_release_key(gedp, pending.first, GED_EDIT_SESSION_INVALIDATE,
	    replacement, pending.second.reason);
    }
    if (!already_closing)
	gi->edit_closing = false;
}


static void
_edit_database_event_observer_ensure(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;
    Ged_Internal *gi = gedp->i->i;
    if (!gi->edit_event_observer)
	gi->edit_event_observer = ged_event_observer_add(gedp,
	    GED_EVENT_OBSERVER_INTERNAL, _edit_database_event_cb, gi);
}


/* ------------------------------------------------------------------ *
 * Public API
 * ------------------------------------------------------------------ */

struct rt_edit *
ged_edit_buf_get(struct ged *gedp, const struct db_full_path *dfp)
{
    if (!gedp || !dfp || !gedp->i || !gedp->i->i)
	return NULL;

    Ged_Internal *gi = gedp->i->i;
    std::string key = _edit_buf_key(dfp);
    auto it = gi->edit_buf.find(key);
    if (it == gi->edit_buf.end())
	return NULL;
    return it->second.s;
}


void
ged_edit_buf_set(struct ged *gedp, const struct db_full_path *dfp, struct rt_edit *s)
{
    if (!gedp || !dfp || !s || !gedp->i || !gedp->i->i)
	return;

    Ged_Internal *gi = gedp->i->i;
    std::string key = _edit_buf_key(dfp);

    /* If an entry already exists, free it first */
    auto it = gi->edit_buf.find(key);
    if (it != gi->edit_buf.end()) {
	Ged_Internal::ged_edit_buf_entry old_entry = it->second;
	_edit_buf_release_scope(gedp, it->second.edit_scope,
	    GED_SCENE_EDIT_CANCEL);
	rt_edit_destroy(it->second.s);
	db_free_full_path(&it->second.dfp);
	gi->edit_buf.erase(it);
	_edit_notify(gedp, GED_EDIT_SESSION_CANCEL, key, old_entry);
    }

    Ged_Internal::ged_edit_buf_entry entry;
    db_full_path_init(&entry.dfp);
    db_dup_full_path(&entry.dfp, dfp);
    entry.s = s;
    entry.path_names.reserve(dfp->fp_len);
    entry.path_instances.reserve(dfp->fp_len);
    for (size_t i = 0; i < dfp->fp_len; ++i) {
	struct directory *component = DB_FULL_PATH_GET(dfp, i);
	entry.path_names.emplace_back(component && component->d_namep ?
	    component->d_namep : "");
	entry.path_instances.push_back(DB_FULL_PATH_GET_COMB_INST(dfp, i));
    }
    entry.id = gi->edit_next_id++;
    if (!entry.id)
	entry.id = gi->edit_next_id++;
    /* Primitive edits affect every displayed occurrence of the database
     * object, but must not draw occurrences the user did not already request.
     * The promotion only changes the logical presentation frontier; retained
     * Obol sources continue to own geometry, cache, and PoP residency. */
    struct ged_scene_edit_request request;
    ged_scene_edit_request_init(&request);
    request.path = key.c_str();
    request.view = ged_draw_active_view_ctx(gedp);
    request.occurrences = GED_SCENE_EDIT_ALL_DRAWN_OCCURRENCES;
    request.draw_if_absent = 0;
    request.purpose = "libged-edit-buffer";
    (void)ged_scene_edit_acquire(gedp, &request, &entry.edit_scope, NULL);
    gi->edit_buf[key] = entry;
    _edit_notify(gedp, GED_EDIT_SESSION_BEGIN, key, gi->edit_buf[key]);
}


static void
_edit_buf_notify_update(struct ged *gedp, const struct db_full_path *dfp,
	int command_id, bool geometry_changed)
{
    if (!gedp || !dfp || !gedp->i || !gedp->i->i)
	return;
    Ged_Internal *gi = gedp->i->i;
    const std::string key = _edit_buf_key(dfp);
    auto it = gi->edit_buf.find(key);
    if (it == gi->edit_buf.end())
	return;
    it->second.revision++;
    it->second.last_cmd = command_id;
    if (geometry_changed)
	it->second.dirty = true;
    _edit_notify(gedp, GED_EDIT_SESSION_UPDATE, key, it->second);
}

void
ged_edit_buf_notify_update(struct ged *gedp, const struct db_full_path *dfp,
	int command_id)
{
    _edit_buf_notify_update(gedp, dfp, command_id, true);
}


void
ged_edit_buf_notify_checkpoint(struct ged *gedp,
	const struct db_full_path *dfp)
{
    if (!gedp || !dfp || !gedp->i || !gedp->i->i)
	return;
    Ged_Internal *gi = gedp->i->i;
    const std::string key = _edit_buf_key(dfp);
    auto it = gi->edit_buf.find(key);
    if (it == gi->edit_buf.end())
	return;
    it->second.checkpoint_dirty = it->second.dirty;
    _edit_notify(gedp, GED_EDIT_SESSION_CHECKPOINT, key, it->second);
}


void
ged_edit_buf_notify_revert(struct ged *gedp, const struct db_full_path *dfp)
{
    if (!gedp || !dfp || !gedp->i || !gedp->i->i)
	return;
    Ged_Internal *gi = gedp->i->i;
    const std::string key = _edit_buf_key(dfp);
    auto it = gi->edit_buf.find(key);
    if (it == gi->edit_buf.end())
	return;
    it->second.revision++;
    it->second.last_cmd = 0;
    it->second.dirty = it->second.checkpoint_dirty;
    _edit_notify(gedp, GED_EDIT_SESSION_REVERT, key, it->second);
}


int
ged_edit_buf_promote(struct ged *gedp, const struct db_full_path *dfp)
{
    if (!gedp || !dfp || !gedp->i || !gedp->i->i)
	return BRLCAD_ERROR;

    Ged_Internal *gi = gedp->i->i;
    std::string key = _edit_buf_key(dfp);
    auto it = gi->edit_buf.find(key);
    if (it == gi->edit_buf.end())
	return BRLCAD_ERROR;

    struct rt_edit *s = it->second.s;
    struct directory *dp = DB_FULL_PATH_CUR_DIR(&it->second.dfp);
    const ged_scene_edit_scope_ref edit_scope = it->second.edit_scope;
    it->second.committing = true;
    const Ged_Internal::ged_edit_buf_entry final_entry = it->second;

    int ret = BRLCAD_ERROR;
    if (dp && gedp->dbip) {
	/* rt_db_put_internal frees es_int internals and re-initializes it
	 * via RT_DB_INTERNAL_INIT, so the subsequent rt_edit_destroy is safe. */
	int put_ret = rt_db_put_internal(dp, gedp->dbip, &s->es_int);
	ret = (put_ret < 0) ? BRLCAD_ERROR : BRLCAD_OK;
    }

    rt_edit_destroy(s);
    db_free_full_path(&it->second.dfp);
    gi->edit_buf.erase(it);

    if (ret == BRLCAD_OK) {
	ged_db_index_refresh(gedp);
	ged_event_notify_object_modified(gedp, dp->d_namep, 1, NULL);
    }
    _edit_buf_release_scope(gedp, edit_scope,
	ret == BRLCAD_OK ? GED_SCENE_EDIT_COMMIT :
	GED_SCENE_EDIT_CANCEL);
    _edit_notify(gedp,
	ret == BRLCAD_OK ? GED_EDIT_SESSION_COMMIT : GED_EDIT_SESSION_CANCEL,
	key, final_entry);

    return ret;
}


void
ged_edit_buf_abandon(struct ged *gedp, const struct db_full_path *dfp)
{
    if (!gedp || !dfp || !gedp->i || !gedp->i->i)
	return;
    std::string key = _edit_buf_key(dfp);
    _edit_release_key(gedp, key, GED_EDIT_SESSION_CANCEL);
}


void
ged_edit_buf_flush(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;

    Ged_Internal *gi = gedp->i->i;
    if (gi->edit_buf.empty())
	return;

    int event_batch_opened = (ged_event_batch_begin(gedp) == GED_EVENT_OK);

    /* Collect keys up front to avoid iterator invalidation during erase */
    std::vector<std::string> keys;
    keys.reserve(gi->edit_buf.size());
    for (auto &kv : gi->edit_buf)
	keys.push_back(kv.first);

    bool any_written = false;
    for (const std::string &key : keys) {
	auto it = gi->edit_buf.find(key);
	if (it == gi->edit_buf.end())
	    continue;

	struct rt_edit *s = it->second.s;
	struct directory *dp = DB_FULL_PATH_CUR_DIR(&it->second.dfp);
	const ged_scene_edit_scope_ref edit_scope = it->second.edit_scope;
	it->second.committing = true;
	const Ged_Internal::ged_edit_buf_entry final_entry = it->second;
	int committed = 0;

	if (dp && gedp->dbip) {
	    if (rt_db_put_internal(dp, gedp->dbip, &s->es_int) >= 0) {
		any_written = true;
		committed = 1;
		(void)ged_event_notify_object_modified(gedp, dp->d_namep, 1, NULL);
	    }
	}

	rt_edit_destroy(s);
	db_free_full_path(&it->second.dfp);
	gi->edit_buf.erase(it);
	_edit_buf_release_scope(gedp, edit_scope,
	    committed ? GED_SCENE_EDIT_COMMIT : GED_SCENE_EDIT_CANCEL);
	_edit_notify(gedp,
	    committed ? GED_EDIT_SESSION_COMMIT : GED_EDIT_SESSION_CANCEL,
	    key, final_entry);
    }

    if (any_written)
	ged_db_index_refresh(gedp);

    if (event_batch_opened)
	ged_event_batch_end(gedp, NULL);
}


static void
_edit_sessions_cancel_all(struct ged *gedp, bool permanent)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;
    Ged_Internal *gi = gedp->i->i;
    if (gi->edit_closing)
	return;
    gi->edit_closing = true;

    std::vector<std::string> keys;
    keys.reserve(gi->edit_buf.size());
    for (const auto &entry : gi->edit_buf)
	keys.push_back(entry.first);

    for (const std::string &key : keys)
	_edit_release_key(gedp, key, GED_EDIT_SESSION_INVALIDATE, nullptr,
	    GED_EDIT_INVALIDATION_DATABASE_CLOSED);
    if (!permanent)
	gi->edit_closing = false;
}


void
ged_edit_sessions_database_close_private(struct ged *gedp)
{
    _edit_sessions_cancel_all(gedp, false);
}


void
ged_edit_sessions_close_private(struct ged *gedp)
{
    _edit_sessions_cancel_all(gedp, true);
}


int
ged_edit_session_ref_is_null(ged_edit_session_ref session)
{
    return !session.owner && !session.id && !session.generation;
}


enum ged_edit_status
ged_edit_session_find(struct ged *gedp, const char *path,
	ged_edit_session_ref *session)
{
    if (session)
	*session = GED_EDIT_SESSION_REF_NULL;
    if (!gedp || !path || !path[0] || !session ||
	!gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;

    Ged_Internal *gi = gedp->i->i;
    std::string key;
    if (!_edit_resolve_key(key, gedp, path, NULL))
	return GED_EDIT_NOT_FOUND;
    auto it = gi->edit_buf.find(key);
    if (it == gi->edit_buf.end())
	return GED_EDIT_NOT_FOUND;
    *session = _edit_ref(gi, it->second);
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_begin(struct ged *gedp, const char *path,
	struct ged_view_context *view, ged_edit_session_ref *session)
{
    if (session)
	*session = GED_EDIT_SESSION_REF_NULL;
    if (!gedp || !path || !path[0] || !session ||
	!gedp->i || !gedp->i->i || !gedp->dbip)
	return GED_EDIT_INVALID;

    _edit_database_event_observer_ensure(gedp);

    struct db_full_path dfp;
    std::string key;
    if (!_edit_resolve_key(key, gedp, path, &dfp))
	return GED_EDIT_NOT_FOUND;

    Ged_Internal *gi = gedp->i->i;
    if (gi->edit_closing) {
	db_free_full_path(&dfp);
	return GED_EDIT_CONFLICT;
    }
    auto existing = gi->edit_buf.find(key);
    if (existing != gi->edit_buf.end()) {
	*session = _edit_ref(gi, existing->second);
	db_free_full_path(&dfp);
	return GED_EDIT_OK;
    }

    for (const auto &active : gi->edit_buf) {
	if (_edit_same_terminal(active.second, &dfp)) {
	    db_free_full_path(&dfp);
	    return GED_EDIT_CONFLICT;
	}
    }

    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct rt_edit *state = rt_edit_create_context(&dfp, gedp->dbip, &tol,
	view ? ged_view_context_bv(view) : NULL);
    if (!state) {
	db_free_full_path(&dfp);
	return GED_EDIT_ERROR;
    }
    ged_edit_buf_set(gedp, &dfp, state);
    auto added = gi->edit_buf.find(key);
    if (added == gi->edit_buf.end()) {
	/* A synchronous BEGIN observer may deliberately cancel the session.  In
	 * that case the buffer has already destroyed state. */
	db_free_full_path(&dfp);
	return GED_EDIT_ERROR;
    }
    *session = _edit_ref(gi, added->second);
    db_free_full_path(&dfp);
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_info_get(struct ged *gedp, ged_edit_session_ref session,
	struct ged_edit_session_info *info)
{
    if (!gedp || !info || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;

    info->session = session;
    info->revision = it->second.revision;
    info->command_id = it->second.last_cmd;
    info->primitive_type = it->second.s ?
	it->second.s->es_int.idb_type : -1;
    info->dirty = it->second.dirty ? 1 : 0;
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_path_get(struct ged *gedp, ged_edit_session_ref session,
	struct bu_vls *path)
{
    if (!gedp || !path || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    bu_vls_strcat(path, it->first.c_str());
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_descriptor_get(struct ged *gedp,
	ged_edit_session_ref session,
	const struct rt_edit_prim_desc **descriptor)
{
    if (descriptor)
	*descriptor = NULL;
    if (!gedp || !descriptor || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    const int type = it->second.s ? it->second.s->es_int.idb_type : -1;
    if (type < 0 || EDOBJ[type].magic != RT_FUNCTAB_MAGIC ||
	!EDOBJ[type].ft_edit_desc)
	return GED_EDIT_UNSUPPORTED;
    *descriptor = EDOBJ[type].ft_edit_desc();
    return *descriptor ? GED_EDIT_OK : GED_EDIT_UNSUPPORTED;
}


enum ged_edit_status
ged_edit_session_command_values_get(struct ged *gedp,
	ged_edit_session_ref session, int command_id,
	struct rt_edit_cmd_values *result)
{
    if (result)
	rt_edit_cmd_values_init(result);
    if (!gedp || !result || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    const int type = it->second.s ? it->second.s->es_int.idb_type : -1;
    if (type < 0 || EDOBJ[type].magic != RT_FUNCTAB_MAGIC)
	return GED_EDIT_UNSUPPORTED;
    const enum rt_edit_value_status status = rt_edit_cmd_values_get(
	it->second.s, command_id, result);
    if (status == RT_EDIT_VALUE_ERROR)
	return GED_EDIT_ERROR;
    if (status == RT_EDIT_VALUE_UNAVAILABLE)
	return GED_EDIT_UNSUPPORTED;
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_parameter_bounds_get(struct ged *gedp,
	ged_edit_session_ref session, int command_id, int parameter_index,
	struct rt_edit_param_bounds *bounds)
{
    if (bounds)
	memset(bounds, 0, sizeof(*bounds));
    if (!gedp || !bounds || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    return rt_edit_param_bounds_get(it->second.s, command_id,
	parameter_index, bounds) == BRLCAD_OK ? GED_EDIT_OK :
	GED_EDIT_UNSUPPORTED;
}


static const struct rt_edit_cmd_desc *
_edit_command(const struct rt_edit_prim_desc *descriptor, int command_id)
{
    if (!descriptor)
	return NULL;
    for (int i = 0; i < descriptor->ncmd; i++) {
	if (descriptor->cmds[i].cmd_id == command_id)
	    return &descriptor->cmds[i];
    }
    return NULL;
}


static size_t
_edit_parameter_width(const struct rt_edit_param_desc *param)
{
    if (!param || param->type == RT_EDIT_PARAM_STRING)
	return 0;
    if (param->type == RT_EDIT_PARAM_POINT2 ||
	param->type == RT_EDIT_PARAM_VECTOR2)
	return 2;
    if (param->type == RT_EDIT_PARAM_POINT ||
	param->type == RT_EDIT_PARAM_VECTOR ||
	param->type == RT_EDIT_PARAM_COLOR)
	return 3;
    if (param->type == RT_EDIT_PARAM_MATRIX)
	return 16;
    if (param->type == RT_EDIT_PARAM_SCALAR_LIST ||
	param->type == RT_EDIT_PARAM_INTEGER_LIST)
	return (size_t)param->nenum;
    return 1;
}


static size_t
_edit_command_value_count(const struct rt_edit_cmd_desc *command)
{
    size_t count = 0;
    if (!command)
	return 0;
    for (int i = 0; i < command->nparam; i++) {
	const struct rt_edit_param_desc *param = &command->params[i];
	if (param->type == RT_EDIT_PARAM_STRING)
	    continue;
	const size_t width = _edit_parameter_width(param);
	const size_t end = (size_t)param->index + width;
	if (end > count)
	    count = end;
    }
    return count;
}


static const struct rt_edit_param_desc *
_edit_command_list_param(const struct rt_edit_cmd_desc *command)
{
    if (!command)
	return NULL;
    for (int i = 0; i < command->nparam; i++) {
	const struct rt_edit_param_desc *param = &command->params[i];
	if (param->type == RT_EDIT_PARAM_SCALAR_LIST ||
	    param->type == RT_EDIT_PARAM_INTEGER_LIST)
	    return param;
    }
    return NULL;
}


static int
_edit_command_values_valid(const struct rt_edit *state,
	const struct rt_edit_cmd_desc *command,
	const struct ged_edit_command_input *input)
{
    if (!state || !command || !input)
	return 0;
    const size_t required = _edit_command_value_count(command);

    const struct rt_edit_param_desc *list_param =
	_edit_command_list_param(command);
    if (required > RT_EDIT_MAXPARA || input->value_count > RT_EDIT_MAXPARA ||
	required > input->value_count ||
	(required && !input->values))
	return 0;
    if (!list_param && input->value_count != required)
	return 0;
    if (list_param) {
	if (list_param->index < 0 ||
	    (size_t)list_param->index > input->value_count)
	    return 0;
	const size_t list_count = input->value_count -
	    (size_t)list_param->index;
	if (list_count < (size_t)list_param->nenum ||
	    list_count > (size_t)(RT_EDIT_MAXPARA - list_param->index))
	    return 0;
    }

    for (int i = 0; i < command->nparam; i++) {
	const struct rt_edit_param_desc *param = &command->params[i];
	if (param->type == RT_EDIT_PARAM_STRING) {
	    if (param->index < 0 || (size_t)param->index >= input->string_count ||
		!input->strings || !input->strings[param->index])
		return 0;
	    continue;
	}
	if (param->index < 0 || (size_t)param->index >= input->value_count)
	    return 0;
	size_t value_end = (size_t)param->index +
	    _edit_parameter_width(param);
	if (param->type == RT_EDIT_PARAM_SCALAR_LIST ||
	    param->type == RT_EDIT_PARAM_INTEGER_LIST)
	    value_end = input->value_count;
	if (value_end > input->value_count)
	    return 0;
	struct rt_edit_param_bounds bounds = {};
	if (rt_edit_param_bounds_get(state, command->cmd_id, i, &bounds) !=
	    BRLCAD_OK)
	    return 0;
	for (size_t vi = (size_t)param->index; vi < value_end; vi++) {
	    const fastf_t value = input->values[vi];
	    if ((bounds.has_minimum && value < bounds.minimum) ||
		(bounds.has_maximum && value > bounds.maximum))
		return 0;
	    if ((param->type == RT_EDIT_PARAM_INTEGER ||
		 param->type == RT_EDIT_PARAM_INTEGER_LIST ||
		 param->type == RT_EDIT_PARAM_ENUM ||
		 param->type == RT_EDIT_PARAM_BOOLEAN) &&
		fabs(value - nearbyint(value)) > SMALL_FASTF)
		return 0;
	    if (param->type == RT_EDIT_PARAM_BOOLEAN &&
		(value < 0.0 || value > 1.0))
		return 0;
	}
    }
    return 1;
}


enum ged_edit_status
ged_edit_session_apply(struct ged *gedp, ged_edit_session_ref session,
	const struct ged_edit_command_input *input)
{
    if (!gedp || !input || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;

    struct rt_edit *state = it->second.s;
    const int type = state ? state->es_int.idb_type : -1;
    if (type < 0 || EDOBJ[type].magic != RT_FUNCTAB_MAGIC ||
	!EDOBJ[type].ft_edit_desc)
	return GED_EDIT_UNSUPPORTED;
    const struct rt_edit_prim_desc *descriptor = EDOBJ[type].ft_edit_desc();
    const struct rt_edit_cmd_desc *command =
	_edit_command(descriptor, input->command_id);
    if (!command)
	return GED_EDIT_UNSUPPORTED;
    if (rt_edit_cmd_control_class(descriptor, command) ==
	    RT_EDIT_CONTROL_UNSUPPORTED)
	return GED_EDIT_UNSUPPORTED;
    if (!_edit_command_values_valid(state, command, input))
	return GED_EDIT_INVALID;

    memset(state->e_para, 0, sizeof(state->e_para));
    memset(state->e_str, 0, sizeof(state->e_str));
    state->e_inpara = (int)input->value_count;
    state->e_nstr = 0;
    if (state->e_inpara)
	memcpy(state->e_para, input->values,
	    (size_t)state->e_inpara * sizeof(fastf_t));
    for (int i = 0; i < command->nparam; i++) {
	const struct rt_edit_param_desc *param = &command->params[i];
	if (param->type == RT_EDIT_PARAM_STRING)
	    rt_edit_set_str(state, param->index,
		input->strings[param->index]);
    }
    if (input->view) {
	struct rt_edit_view edit_view;
	rt_edit_view_from_context(&edit_view,
	    ged_view_context_bv(input->view));
	rt_edit_set_view(state, &edit_view);
    }

    rt_edit_set_edflag(state, command->cmd_id);
    if (rt_edit_process_result(state) != BRLCAD_OK) {
	if (state->log_str && bu_vls_strlen(state->log_str)) {
	    if (gedp->ged_result_str)
		bu_vls_printf(gedp->ged_result_str, "%s",
		    bu_vls_cstr(state->log_str));
	    bu_vls_trunc(state->log_str, 0);
	}
	return GED_EDIT_REJECTED;
    }
    /* Selection changes are revisioned session state, but they do not modify
     * database geometry and must not turn a clean session into a pending
     * database write.  Presenters and command queries still receive the same
     * typed UPDATE event. */
    const bool geometry_changed = !command->category ||
	!BU_STR_EQUAL(command->category, "selection");
    _edit_buf_notify_update(gedp, &it->second.dfp, command->cmd_id,
	geometry_changed);
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_checkpoint(struct ged *gedp, ged_edit_session_ref session)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    if (rt_edit_checkpoint(it->second.s) != BRLCAD_OK)
	return GED_EDIT_ERROR;
    ged_edit_buf_notify_checkpoint(gedp, &it->second.dfp);
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_revert(struct ged *gedp, ged_edit_session_ref session)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    if (rt_edit_revert(it->second.s) != BRLCAD_OK)
	return GED_EDIT_ERROR;
    ged_edit_buf_notify_revert(gedp, &it->second.dfp);
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_commit(struct ged *gedp, ged_edit_session_ref session)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    if (!it->second.dirty) {
	/* Selection and other session-only state can advance revisions without
	 * changing the primitive.  Closing such a session is still a successful
	 * commit lifecycle event, but must not rewrite the database or emit a
	 * false object-modified notification. */
	const std::string key = it->first;
	const ged_scene_edit_scope_ref scope = it->second.edit_scope;
	const Ged_Internal::ged_edit_buf_entry finalEntry = it->second;
	rt_edit_destroy(it->second.s);
	db_free_full_path(&it->second.dfp);
	gi->edit_buf.erase(it);
	_edit_buf_release_scope(gedp, scope, GED_SCENE_EDIT_CANCEL);
	_edit_notify(gedp, GED_EDIT_SESSION_COMMIT, key, finalEntry);
	return GED_EDIT_OK;
    }
    struct db_full_path path;
    db_full_path_init(&path);
    db_dup_full_path(&path, &it->second.dfp);
    const int ret = ged_edit_buf_promote(gedp, &path);
    db_free_full_path(&path);
    return ret == BRLCAD_OK ? GED_EDIT_OK : GED_EDIT_ERROR;
}


enum ged_edit_status
ged_edit_session_cancel(struct ged *gedp, ged_edit_session_ref session)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    struct db_full_path path;
    db_full_path_init(&path);
    db_dup_full_path(&path, &it->second.dfp);
    ged_edit_buf_abandon(gedp, &path);
    db_free_full_path(&path);
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_internal_copy(struct ged *gedp,
	ged_edit_session_ref session, struct rt_db_internal *intern)
{
    if (!gedp || !intern || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *gi = gedp->i->i;
    auto it = _edit_find_ref(gi, session);
    if (it == gi->edit_buf.end())
	return session.owner == _edit_owner(gi) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;

    struct bu_external external;
    BU_EXTERNAL_INIT(&external);
    if (rt_obj_export(&external, &it->second.s->es_int, 1.0,
	    it->second.s->dbip) < 0)
	return GED_EDIT_ERROR;
    RT_DB_INTERNAL_INIT(intern);
    intern->idb_minor_type = it->second.s->es_int.idb_type;
    mat_t identity;
    MAT_IDN(identity);
    const int ret = rt_obj_import(intern, &external, identity,
	it->second.s->dbip);
    bu_free_external(&external);
    return ret < 0 ? GED_EDIT_ERROR : GED_EDIT_OK;
}


ged_edit_observer_token
ged_edit_observer_add(struct ged *gedp, ged_edit_observer_func_t callback,
	void *client_data)
{
    if (!gedp || !callback || !gedp->i || !gedp->i->i)
	return 0;
    Ged_Internal *gi = gedp->i->i;
    ged_edit_observer_token token = gi->edit_next_observer++;
    if (!token)
	token = gi->edit_next_observer++;
    gi->edit_observers[token] = std::make_pair(callback, client_data);
    return token;
}


int
ged_edit_observer_remove(struct ged *gedp, ged_edit_observer_token token)
{
    if (!gedp || !token || !gedp->i || !gedp->i->i)
	return 0;
    return gedp->i->i->edit_observers.erase(token) ? 1 : 0;
}


/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
