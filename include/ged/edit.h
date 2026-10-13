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
 * GED edit-session queries and notifications.
 *
 * Each session owns one intermediate geometry variant.  Applications may
 * create multiple sessions for the same path and decide which clients, if any,
 * share a session.  References become stale when their session is committed,
 * canceled, or invalidated.
 */

#ifndef GED_EDIT_H
#define GED_EDIT_H

#include <stddef.h>
#include <stdint.h>

#include "vmath.h"
#include "ged/defines.h"

struct bu_vls;
struct ged;
struct rt_db_internal;
struct rt_edit_cmd_values;
struct rt_edit_param_bounds;
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

/** Result of an edit-session API request. */
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

/** Database entity changed by an edit session. */
enum ged_edit_session_target {
    GED_EDIT_TARGET_OBJECT = 0,
    GED_EDIT_TARGET_OCCURRENCE
};

/** Semantic transition reported to edit-session observers. */
enum ged_edit_session_event_kind {
    GED_EDIT_SESSION_BEGIN = 0,
    GED_EDIT_SESSION_UPDATE,
    GED_EDIT_SESSION_CHECKPOINT,
    GED_EDIT_SESSION_REVERT,
    GED_EDIT_SESSION_COMMIT,
    GED_EDIT_SESSION_CANCEL,
    /**
     * The database object or an ancestor of the edited path changed outside
     * this session.  The intermediate state has been discarded and the
     * session reference is stale.  A rename may supply replacement_path so a
     * persistent presenter can attach to the renamed object without retaining
     * the stale transaction.
     */
    GED_EDIT_SESSION_INVALIDATE
};

/** Why GED invalidated an edit session. */
enum ged_edit_session_invalidation_reason {
    GED_EDIT_INVALIDATION_NONE = 0,
    GED_EDIT_INVALIDATION_SOURCE_CHANGED,
    GED_EDIT_INVALIDATION_SOURCE_RENAMED,
    GED_EDIT_INVALIDATION_SOURCE_REMOVED,
    GED_EDIT_INVALIDATION_PATH_REMOVED,
    GED_EDIT_INVALIDATION_DATABASE_CLOSED
};

/** Callback-lifetime description of one edit-session transition. */
struct ged_edit_session_event {
    enum ged_edit_session_event_kind kind;
    ged_edit_session_ref session;
    const char *path;       /**< Borrowed canonical path. */
    const char *replacement_path; /**< Borrowed renamed path, or NULL. */
    enum ged_edit_session_invalidation_reason invalidation_reason;
    uint64_t revision;      /**< Monotonic within this session. */
    int command_id;         /**< Librt edit command, or zero if inapplicable. */
};

/** Lightweight current-state information for an active session. */
struct ged_edit_session_info {
    ged_edit_session_ref session;
    enum ged_edit_session_target target;
    uint64_t revision;
    int command_id;
    int primitive_type;
    int dirty;
};

/** Typed input for one librt descriptor command. */
struct ged_edit_command_input {
    int command_id;
    const fastf_t *values; /**< Indexed as described by rt_edit_param_desc. */
    size_t value_count;
    const char *const *strings; /**< Indexed string parameters. */
    size_t string_count;
};

typedef uint64_t ged_edit_observer_token;

typedef void (*ged_edit_observer_func_t)(
    struct ged *gedp,
    const struct ged_edit_session_event *event,
    void *client_data);

__BEGIN_DECLS

/** Return non-zero when @p session is the null reference. */
GED_EXPORT extern int
ged_edit_session_ref_is_null(ged_edit_session_ref session);

/**
 * Begin an independent edit session for @p path and @p target.
 *
 * Object sessions change terminal object geometry.  Occurrence sessions change
 * only the matrix on the final parent-to-child combination arc and require an
 * instance-qualified path when duplicate child names are present.
 */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_begin(struct ged *gedp, const char *path,
	enum ged_edit_session_target target, ged_edit_session_ref *session);

/**
 * List active sessions for @p path and @p target.
 *
 * @p count receives the total number of matching sessions.  Up to @p capacity
 * references are written to @p sessions.  Pass NULL and zero to query only the
 * count.  The result order is unspecified.
 */
GED_EXPORT extern enum ged_edit_status
ged_edit_sessions_get(struct ged *gedp, const char *path,
	enum ged_edit_session_target target, ged_edit_session_ref *sessions,
	size_t capacity, size_t *count);

/** Read the current final-arc matrix from an occurrence session. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_occurrence_matrix_get(struct ged *gedp,
	ged_edit_session_ref session, mat_t matrix);

/** Replace the final-arc matrix in an occurrence session. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_occurrence_matrix_set(struct ged *gedp,
	ged_edit_session_ref session, const mat_t matrix);

/** Return current scalar session information. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_info_get(struct ged *gedp, ged_edit_session_ref session,
	struct ged_edit_session_info *info);

/** Append the canonical session path to caller-owned @p path. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_path_get(struct ged *gedp, ged_edit_session_ref session,
	struct bu_vls *path);

/** Return the borrowed static librt descriptor for an active session. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_descriptor_get(struct ged *gedp,
	ged_edit_session_ref session,
	const struct rt_edit_prim_desc **descriptor);

/** Read current numeric and string values for @p command_id. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_command_values_get(struct ged *gedp,
	ged_edit_session_ref session, int command_id,
	struct rt_edit_cmd_values *result);

/** Read effective static/current-topology bounds for one command parameter. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_parameter_bounds_get(struct ged *gedp,
	ged_edit_session_ref session, int command_id, int parameter_index,
	struct rt_edit_param_bounds *bounds);

/** Apply one typed descriptor command to the intermediate session. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_apply(struct ged *gedp, ged_edit_session_ref session,
	const struct ged_edit_command_input *input);

/** Save a single-level restore point for the intermediate state. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_checkpoint(struct ged *gedp, ged_edit_session_ref session);

/** Restore the most recent checkpoint without closing the session. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_revert(struct ged *gedp, ged_edit_session_ref session);

/** Commit intermediate geometry to the database and close the session. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_commit(struct ged *gedp, ged_edit_session_ref session);

/** Discard intermediate geometry and close the session. */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_cancel(struct ged *gedp, ged_edit_session_ref session);

/**
 * Copy intermediate object geometry into caller-owned @p intern.
 * Initialize @p intern before calling and release a successful result with
 * rt_db_free_internal().  Occurrence sessions return GED_EDIT_UNSUPPORTED;
 * query their matrix instead.
 */
GED_EXPORT extern enum ged_edit_status
ged_edit_session_internal_copy(struct ged *gedp,
	ged_edit_session_ref session, struct rt_db_internal *intern);

/** Register a synchronous edit-session observer. */
GED_EXPORT extern ged_edit_observer_token
ged_edit_observer_add(struct ged *gedp, ged_edit_observer_func_t callback,
	void *client_data);

/** Remove an observer; removal from inside its callback is supported. */
GED_EXPORT extern int
ged_edit_observer_remove(struct ged *gedp, ged_edit_observer_token token);

struct bn_tol;
struct bu_list;
struct bview;
struct db_full_path;
struct rt_edit;

enum ged_edit_mode {
    GED_EDIT_VIEW,
    GED_EDIT_PRIMITIVE,
    GED_EDIT_INSTANCE
};

/** Begin an edit of path.  Instance edits affect the object at path_index.
 * A non-NULL edit is borrowed and must outlive the session.  Otherwise the
 * session owns a persistent rt_edit.  Tolerance is copied; view is borrowed.
 * Beginning an edit while another is active is an error.
 */
GED_EXPORT int ged_edit_interactive_begin(struct ged *gedp, struct rt_edit *edit,
    const struct db_full_path *path, enum ged_edit_mode mode, size_t path_index,
    const struct bn_tol *tol, struct bview *view);

/** Return the current editing mode, or GED_EDIT_VIEW when idle. */
GED_EXPORT enum ged_edit_mode ged_edit_interactive_mode(const struct ged *gedp);

/** Return the active working edit, or NULL.  The session retains ownership.
 * Owned edits follow ged_gvp; borrowed edits retain their host's view binding.
 */
GED_EXPORT struct rt_edit *ged_edit_interactive_edit(struct ged *gedp);

/** Return the active selection path, or NULL.  The path is borrowed. */
GED_EXPORT const struct db_full_path *ged_edit_interactive_path(const struct ged *gedp);

/** Return the selected object index in the full path, or SIZE_MAX when idle. */
GED_EXPORT size_t ged_edit_interactive_path_index(const struct ged *gedp);

/** Return nonzero when path names the active edit target in mode.  For an
 * instance edit this compares the selected ancestor, not its preview leaf.
 */
GED_EXPORT int ged_edit_interactive_matches(const struct ged *gedp,
    const char *path, enum ged_edit_mode mode);

/** Return a borrowed active menu label, or NULL when idle or past the end.
 * Row zero is the heading.  Other labels can be selected with press.
 */
GED_EXPORT const char *ged_edit_interactive_menu_label(const struct ged *gedp, size_t row);

/** Write the working edit.  A nonzero continue_editing retains the session
 * and establishes a new reset baseline.  Failure retains the working edit.
 */
GED_EXPORT int ged_edit_interactive_apply(struct ged *gedp, int continue_editing);

/** Reload the database geometry and retain the active editing session. */
GED_EXPORT int ged_edit_interactive_reset(struct ged *gedp);

/** Abandon working edits without writing the database.  Safe when idle. */
GED_EXPORT void ged_edit_interactive_end(struct ged *gedp);

/** Process the active rt_edit inputs and update its working geometry. */
GED_EXPORT int ged_edit_interactive_process(struct ged *gedp);

/** Apply a normalized view-space mouse position to the active operation. */
GED_EXPORT int ged_edit_interactive_mouse(struct ged *gedp, const vect_t position);

/** Dispatch contextual transforms to the active edit, otherwise execute the
 * ordinary GED command.  Hosts can call ged_exec directly for view actions.
 */
GED_EXPORT int ged_edit_interactive_exec(struct ged *gedp, int argc, const char *argv[]);

/** Advance active knob rates by steps nominal event-loop ticks.  Hosts choose
 * their tick interval; a zero or negative step does nothing.  Returns 1 when
 * geometry/view changed, 0 when idle, and -1 when an operation fails.
 */
GED_EXPORT int ged_edit_interactive_tick(struct ged *gedp, fastf_t steps);

/** Plot a displayed path using working geometry.  Returns 1 for a preview,
 * 0 for an unaffected path, and -1 on failure.  Output is appended to vlist.
 * Primitive edits preview every use; instance edits preview only the selected
 * occurrence and its descendants.  No database or display objects are changed.
 */
GED_EXPORT int ged_edit_interactive_plot(struct ged *gedp, const struct db_full_path *path,
    struct bview *view, struct bu_list *vlist);

/** Execute a classic editing command using the shared session. */
GED_EXPORT int ged_edit_interactive_command(struct ged *gedp, int argc, const char *argv[]);

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
