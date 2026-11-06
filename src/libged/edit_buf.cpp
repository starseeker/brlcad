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

#include <atomic>
#include <string>
#include <vector>

#include "raytrace.h"
#include "rt/db_fullpath.h"
#include "rt/edit.h"
#include "rt/func.h"
#include "rt/functab.h"

#include "ged/edit.h"

#include "./ged_private.h"
#include "./dbi.h"


/* ------------------------------------------------------------------ *
 * Ged_Internal destructor — abandons all pending edits on close.
 * ------------------------------------------------------------------ */

Ged_Internal::~Ged_Internal()
{
    for (auto &kv : edit_buf) {
	rt_edit_destroy(kv.second.s);
	db_free_full_path(&kv.second.dfp);
    }
    edit_buf.clear();
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

static uint64_t
_edit_owner(Ged_Internal *internal)
{
    static std::atomic<uint64_t> next_owner(1);
    if (internal && !internal->edit_owner) {
	internal->edit_owner = next_owner.fetch_add(1);
	if (!internal->edit_owner)
	    internal->edit_owner = next_owner.fetch_add(1);
    }
    return internal ? internal->edit_owner : 0;
}

static ged_edit_session_ref
_edit_ref(Ged_Internal *internal,
	const Ged_Internal::ged_edit_buf_entry &entry)
{
    ged_edit_session_ref ref = {_edit_owner(internal), entry.id,
	entry.generation};
    return ref;
}

static std::unordered_map<std::string,
    Ged_Internal::ged_edit_buf_entry>::iterator
_edit_find_ref(Ged_Internal *internal, ged_edit_session_ref ref)
{
    if (!internal || ref.owner != _edit_owner(internal) ||
	!ref.id || !ref.generation)
	return internal ? internal->edit_buf.end() :
	    std::unordered_map<std::string,
		Ged_Internal::ged_edit_buf_entry>::iterator();
    for (auto entry = internal->edit_buf.begin();
	entry != internal->edit_buf.end(); ++entry) {
	if (entry->second.id == ref.id &&
	    entry->second.generation == ref.generation)
	    return entry;
    }
    return internal->edit_buf.end();
}

static bool
_edit_resolve_key(std::string &key, struct ged *gedp, const char *path,
	struct db_full_path *resolved)
{
    if (!gedp || !gedp->dbip || !path || !path[0])
	return false;
    struct db_full_path local_path;
    db_full_path_init(&local_path);
    if (db_string_to_path(&local_path, gedp->dbip, path) != 0) {
	db_free_full_path(&local_path);
	return false;
    }
    key = _edit_buf_key(&local_path);
    if (resolved)
	db_dup_full_path(resolved, &local_path);
    db_free_full_path(&local_path);
    return true;
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
	rt_edit_destroy(it->second.s);
	db_free_full_path(&it->second.dfp);
	gi->edit_buf.erase(it);
    }

    Ged_Internal::ged_edit_buf_entry entry;
    db_full_path_init(&entry.dfp);
    db_dup_full_path(&entry.dfp, dfp);
    entry.s = s;
    entry.id = gi->edit_next_id++;
    if (!entry.id)
	entry.id = gi->edit_next_id++;
    gi->edit_buf[key] = entry;
}


void
ged_edit_buf_notify_update(struct ged *gedp, const struct db_full_path *dfp,
	int command_id)
{
    if (!gedp || !dfp || !gedp->i || !gedp->i->i)
	return;
    Ged_Internal *internal = gedp->i->i;
    auto entry = internal->edit_buf.find(_edit_buf_key(dfp));
    if (entry == internal->edit_buf.end())
	return;
    entry->second.revision++;
    entry->second.last_cmd = command_id;
    entry->second.dirty = true;
}


void
ged_edit_buf_notify_checkpoint(struct ged *gedp,
	const struct db_full_path *dfp)
{
    if (!gedp || !dfp || !gedp->i || !gedp->i->i)
	return;
    Ged_Internal *internal = gedp->i->i;
    auto entry = internal->edit_buf.find(_edit_buf_key(dfp));
    if (entry != internal->edit_buf.end())
	entry->second.checkpoint_dirty = entry->second.dirty;
}


void
ged_edit_buf_notify_revert(struct ged *gedp,
	const struct db_full_path *dfp)
{
    if (!gedp || !dfp || !gedp->i || !gedp->i->i)
	return;
    Ged_Internal *internal = gedp->i->i;
    auto entry = internal->edit_buf.find(_edit_buf_key(dfp));
    if (entry == internal->edit_buf.end())
	return;
    entry->second.revision++;
    entry->second.last_cmd = 0;
    entry->second.dirty = entry->second.checkpoint_dirty;
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

    /* Notify DbiState of the change */
    if (ret == BRLCAD_OK && gedp->dbi_state) {
	DbiState *dbis = (DbiState *)gedp->dbi_state;
	dbis->update();
    }

    return ret;
}


void
ged_edit_buf_abandon(struct ged *gedp, const struct db_full_path *dfp)
{
    if (!gedp || !dfp || !gedp->i || !gedp->i->i)
	return;

    Ged_Internal *gi = gedp->i->i;
    std::string key = _edit_buf_key(dfp);
    auto it = gi->edit_buf.find(key);
    if (it == gi->edit_buf.end())
	return;

    rt_edit_destroy(it->second.s);
    db_free_full_path(&it->second.dfp);
    gi->edit_buf.erase(it);
}


void
ged_edit_buf_flush(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return;

    Ged_Internal *gi = gedp->i->i;
    if (gi->edit_buf.empty())
	return;

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

	if (dp && gedp->dbip) {
	    rt_db_put_internal(dp, gedp->dbip, &s->es_int);
	    any_written = true;
	}

	rt_edit_destroy(s);
	db_free_full_path(&it->second.dfp);
	gi->edit_buf.erase(it);
    }

    if (any_written && gedp->dbi_state) {
	DbiState *dbis = (DbiState *)gedp->dbi_state;
	dbis->update();
    }
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

    std::string key;
    if (!_edit_resolve_key(key, gedp, path, NULL))
	return GED_EDIT_NOT_FOUND;
    Ged_Internal *internal = gedp->i->i;
    auto entry = internal->edit_buf.find(key);
    if (entry == internal->edit_buf.end())
	return GED_EDIT_NOT_FOUND;
    *session = _edit_ref(internal, entry->second);
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_begin(struct ged *gedp, const char *path,
	ged_edit_session_ref *session)
{
    if (session)
	*session = GED_EDIT_SESSION_REF_NULL;
    if (!gedp || !path || !path[0] || !session ||
	!gedp->i || !gedp->i->i || !gedp->dbip)
	return GED_EDIT_INVALID;

    struct db_full_path full_path;
    db_full_path_init(&full_path);
    std::string key;
    if (!_edit_resolve_key(key, gedp, path, &full_path))
	return GED_EDIT_NOT_FOUND;

    Ged_Internal *internal = gedp->i->i;
    auto existing = internal->edit_buf.find(key);
    if (existing != internal->edit_buf.end()) {
	*session = _edit_ref(internal, existing->second);
	db_free_full_path(&full_path);
	return GED_EDIT_OK;
    }

    struct directory *terminal = DB_FULL_PATH_CUR_DIR(&full_path);
    for (const auto &active : internal->edit_buf) {
	if (DB_FULL_PATH_CUR_DIR(&active.second.dfp) == terminal) {
	    db_free_full_path(&full_path);
	    return GED_EDIT_CONFLICT;
	}
    }

    struct bn_tol tolerance = BN_TOL_INIT_TOL;
    struct rt_edit *state = rt_edit_create(&full_path, gedp->dbip,
	&tolerance, NULL);
    if (!state) {
	db_free_full_path(&full_path);
	return GED_EDIT_ERROR;
    }
    ged_edit_buf_set(gedp, &full_path, state);
    auto added = internal->edit_buf.find(key);
    if (added == internal->edit_buf.end()) {
	db_free_full_path(&full_path);
	return GED_EDIT_ERROR;
    }
    *session = _edit_ref(internal, added->second);
    db_free_full_path(&full_path);
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_info_get(struct ged *gedp, ged_edit_session_ref session,
	struct ged_edit_session_info *info)
{
    if (!gedp || !info || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *internal = gedp->i->i;
    auto entry = _edit_find_ref(internal, session);
    if (entry == internal->edit_buf.end())
	return session.owner == _edit_owner(internal) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    info->session = session;
    info->revision = entry->second.revision;
    info->command_id = entry->second.last_cmd;
    info->primitive_type = entry->second.s ?
	entry->second.s->es_int.idb_type : -1;
    info->dirty = entry->second.dirty ? 1 : 0;
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_path_get(struct ged *gedp, ged_edit_session_ref session,
	struct bu_vls *path)
{
    if (!gedp || !path || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *internal = gedp->i->i;
    auto entry = _edit_find_ref(internal, session);
    if (entry == internal->edit_buf.end())
	return session.owner == _edit_owner(internal) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    bu_vls_strcat(path, entry->first.c_str());
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
    Ged_Internal *internal = gedp->i->i;
    auto entry = _edit_find_ref(internal, session);
    if (entry == internal->edit_buf.end())
	return session.owner == _edit_owner(internal) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    const int type = entry->second.s ? entry->second.s->es_int.idb_type : -1;
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
    Ged_Internal *internal = gedp->i->i;
    auto entry = _edit_find_ref(internal, session);
    if (entry == internal->edit_buf.end())
	return session.owner == _edit_owner(internal) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    const enum rt_edit_value_status status = rt_edit_cmd_values_get(
	entry->second.s, command_id, result);
    if (status == RT_EDIT_VALUE_ERROR)
	return GED_EDIT_ERROR;
    return status == RT_EDIT_VALUE_OK ? GED_EDIT_OK : GED_EDIT_UNSUPPORTED;
}


enum ged_edit_status
ged_edit_session_checkpoint(struct ged *gedp, ged_edit_session_ref session)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *internal = gedp->i->i;
    auto entry = _edit_find_ref(internal, session);
    if (entry == internal->edit_buf.end())
	return session.owner == _edit_owner(internal) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    if (rt_edit_checkpoint(entry->second.s) != BRLCAD_OK)
	return GED_EDIT_ERROR;
    ged_edit_buf_notify_checkpoint(gedp, &entry->second.dfp);
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_revert(struct ged *gedp, ged_edit_session_ref session)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *internal = gedp->i->i;
    auto entry = _edit_find_ref(internal, session);
    if (entry == internal->edit_buf.end())
	return session.owner == _edit_owner(internal) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    if (rt_edit_revert(entry->second.s) != BRLCAD_OK)
	return GED_EDIT_ERROR;
    ged_edit_buf_notify_revert(gedp, &entry->second.dfp);
    return GED_EDIT_OK;
}


enum ged_edit_status
ged_edit_session_commit(struct ged *gedp, ged_edit_session_ref session)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *internal = gedp->i->i;
    auto entry = _edit_find_ref(internal, session);
    if (entry == internal->edit_buf.end())
	return session.owner == _edit_owner(internal) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    if (!entry->second.dirty) {
	struct db_full_path path;
	db_full_path_init(&path);
	db_dup_full_path(&path, &entry->second.dfp);
	ged_edit_buf_abandon(gedp, &path);
	db_free_full_path(&path);
	return GED_EDIT_OK;
    }
    struct db_full_path path;
    db_full_path_init(&path);
    db_dup_full_path(&path, &entry->second.dfp);
    const int result = ged_edit_buf_promote(gedp, &path);
    db_free_full_path(&path);
    return result == BRLCAD_OK ? GED_EDIT_OK : GED_EDIT_ERROR;
}


enum ged_edit_status
ged_edit_session_cancel(struct ged *gedp, ged_edit_session_ref session)
{
    if (!gedp || !gedp->i || !gedp->i->i)
	return GED_EDIT_INVALID;
    Ged_Internal *internal = gedp->i->i;
    auto entry = _edit_find_ref(internal, session);
    if (entry == internal->edit_buf.end())
	return session.owner == _edit_owner(internal) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;
    struct db_full_path path;
    db_full_path_init(&path);
    db_dup_full_path(&path, &entry->second.dfp);
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
    Ged_Internal *internal = gedp->i->i;
    auto entry = _edit_find_ref(internal, session);
    if (entry == internal->edit_buf.end())
	return session.owner == _edit_owner(internal) ? GED_EDIT_STALE :
	    GED_EDIT_INVALID;

    struct bu_external external;
    BU_EXTERNAL_INIT(&external);
    if (rt_obj_export(&external, &entry->second.s->es_int, 1.0,
	entry->second.s->dbip) < 0)
	return GED_EDIT_ERROR;
    RT_DB_INTERNAL_INIT(intern);
    intern->idb_minor_type = entry->second.s->es_int.idb_type;
    mat_t identity;
    MAT_IDN(identity);
    const int result = rt_obj_import(intern, &external, identity,
	entry->second.s->dbip);
    bu_free_external(&external);
    return result < 0 ? GED_EDIT_ERROR : GED_EDIT_OK;
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
