/*                         E D I T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file ged/edit.h
 *
 * Path-scoped GED edit sessions.
 */

#ifndef GED_EDIT_H
#define GED_EDIT_H

#include <stdint.h>

#include "ged/defines.h"

struct bu_vls;
struct ged;
struct rt_db_internal;
struct rt_edit_cmd_values;
struct rt_edit_prim_desc;

/** Opaque, generation-checked reference to one active edit session. */
typedef struct ged_edit_session_ref {
    uint64_t owner;
    uint64_t id;
    uint64_t generation;
} ged_edit_session_ref;

#define GED_EDIT_SESSION_REF_NULL_INIT {0, 0, 0}
#ifdef __cplusplus
#  define GED_EDIT_SESSION_REF_NULL ged_edit_session_ref{0, 0, 0}
#else
#  define GED_EDIT_SESSION_REF_NULL ((ged_edit_session_ref){0, 0, 0})
#endif

enum ged_edit_status {
    GED_EDIT_OK = 0,
    GED_EDIT_INVALID,
    GED_EDIT_NOT_FOUND,
    GED_EDIT_STALE,
    GED_EDIT_UNSUPPORTED,
    GED_EDIT_REJECTED,
    GED_EDIT_CONFLICT,
    GED_EDIT_ERROR
};

struct ged_edit_session_info {
    ged_edit_session_ref session;
    uint64_t revision;
    int command_id;
    int primitive_type;
    int dirty;
};

__BEGIN_DECLS

/** Return non-zero when @p session is the null reference. */
GED_EXPORT extern int
ged_edit_session_ref_is_null(ged_edit_session_ref session);

/** Begin a path transaction, or join the existing transaction for that path. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_begin(struct ged *gedp, const char *path,
	ged_edit_session_ref *session);

/** Find an active path transaction without creating one. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_find(struct ged *gedp, const char *path,
	ged_edit_session_ref *session);

/** Return lightweight state for an active transaction. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_info_get(struct ged *gedp, ged_edit_session_ref session,
	struct ged_edit_session_info *info);

/** Append the canonical transaction path to caller-owned @p path. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_path_get(struct ged *gedp, ged_edit_session_ref session,
	struct bu_vls *path);

/** Return borrowed static edit metadata for the transaction primitive. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_descriptor_get(struct ged *gedp,
	ged_edit_session_ref session,
	const struct rt_edit_prim_desc **descriptor);

/** Read typed current values for one descriptor command. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_command_values_get(struct ged *gedp,
	ged_edit_session_ref session, int command_id,
	struct rt_edit_cmd_values *result);

/** Save a single-level restore point without closing the transaction. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_checkpoint(struct ged *gedp, ged_edit_session_ref session);

/** Restore the most recent checkpoint without closing the transaction. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_revert(struct ged *gedp, ged_edit_session_ref session);

/** Write changed intermediate geometry and close the transaction. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_commit(struct ged *gedp, ged_edit_session_ref session);

/** Discard intermediate geometry and close the transaction. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_cancel(struct ged *gedp, ged_edit_session_ref session);

/** Copy the intermediate primitive into initialized caller-owned storage. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_internal_copy(struct ged *gedp,
	ged_edit_session_ref session, struct rt_db_internal *intern);

__END_DECLS

#endif /* GED_EDIT_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
