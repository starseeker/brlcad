/* BRL-CAD
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

#ifndef GED_EDIT_H
#define GED_EDIT_H

#include "common.h"
#include "ged/defines.h"
#include "rt/edit.h"

__BEGIN_DECLS

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
GED_EXPORT int ged_edit_session_begin(struct ged *gedp, struct rt_edit *edit,
    const struct db_full_path *path, enum ged_edit_mode mode, size_t path_index,
    const struct bn_tol *tol, struct bview *view);

/** Return the current editing mode, or GED_EDIT_VIEW when idle. */
GED_EXPORT enum ged_edit_mode ged_edit_session_mode(const struct ged *gedp);

/** Return the active working edit, or NULL.  The session retains ownership.
 * Owned edits follow ged_gvp; borrowed edits retain their host's view binding.
 */
GED_EXPORT struct rt_edit *ged_edit_session_edit(struct ged *gedp);

/** Return the active selection path, or NULL.  The path is borrowed. */
GED_EXPORT const struct db_full_path *ged_edit_session_path(const struct ged *gedp);

/** Return the selected object index in the full path, or SIZE_MAX when idle. */
GED_EXPORT size_t ged_edit_session_path_index(const struct ged *gedp);

/** Return nonzero when path names the active edit target in mode.  For an
 * instance edit this compares the selected ancestor, not its preview leaf.
 */
GED_EXPORT int ged_edit_session_matches(const struct ged *gedp,
    const char *path, enum ged_edit_mode mode);

/** Return a borrowed active menu label, or NULL when idle or past the end.
 * Row zero is the heading.  Other labels can be selected with press.
 */
GED_EXPORT const char *ged_edit_session_menu_label(const struct ged *gedp, size_t row);

/** Write the working edit.  A nonzero continue_editing retains the session
 * and establishes a new reset baseline.  Failure retains the working edit.
 */
GED_EXPORT int ged_edit_session_apply(struct ged *gedp, int continue_editing);

/** Reload the database geometry and retain the active editing session. */
GED_EXPORT int ged_edit_session_reset(struct ged *gedp);

/** Abandon working edits without writing the database.  Safe when idle. */
GED_EXPORT void ged_edit_session_end(struct ged *gedp);

/** Process the active rt_edit inputs and update its working geometry. */
GED_EXPORT int ged_edit_session_process(struct ged *gedp);

/** Apply a normalized view-space mouse position to the active operation. */
GED_EXPORT int ged_edit_session_mouse(struct ged *gedp, const vect_t position);

/** Dispatch contextual transforms to the active edit, otherwise execute the
 * ordinary GED command.  Hosts can call ged_exec directly for view actions.
 */
GED_EXPORT int ged_edit_session_exec(struct ged *gedp, int argc, const char *argv[]);

/** Advance active knob rates by steps nominal event-loop ticks.  Hosts choose
 * their tick interval; a zero or negative step does nothing.  Returns 1 when
 * geometry/view changed, 0 when idle, and -1 when an operation fails.
 */
GED_EXPORT int ged_edit_session_tick(struct ged *gedp, fastf_t steps);

/** Plot a displayed path using working geometry.  Returns 1 for a preview,
 * 0 for an unaffected path, and -1 on failure.  Output is appended to vlist.
 * Primitive edits preview every use; instance edits preview only the selected
 * occurrence and its descendants.  No database or display objects are changed.
 */
GED_EXPORT int ged_edit_session_plot(struct ged *gedp, const struct db_full_path *path,
    struct bview *view, struct bu_list *vlist);

/** Execute a classic editing command using the shared session. */
GED_EXPORT int ged_edit_session_command(struct ged *gedp, int argc, const char *argv[]);

__END_DECLS

#endif

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
