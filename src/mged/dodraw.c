/*                        D O D R A W . C
 * BRL-CAD
 *
 * Copyright (c) 1985-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */

#include "common.h"

#include <string.h>

#include "vmath.h"
#include "bn.h"
#include "nmg.h"
#include "rt/geom.h"		/* for ID_POLY special support */
#include "raytrace.h"
#include "rt/db4.h"
#include "ged/view.h"

#include "./mged.h"
#include "./mged_dm.h"
#include "./cmd.h"

#define MGED_EDIT_PREVIEW_PREFIX "_mged_edit_preview::"

static int
mged_check_shape_ref(struct mged_state *s, ged_draw_shape_ref ref, const char *caller)
{
    if (!s || !s->gedp || ged_draw_shape_ref_is_null(ref)) {
	if (s && s->interp && caller)
	    Tcl_AppendResult(s->interp, caller, "() ref is NULL\n", (char *)NULL);
	return 0;
    }

    struct ged_draw_shape_record rec;
    if (!ged_draw_shape_record_get(s->gedp, ref, &rec)) {
	if (s->interp && caller)
	    Tcl_AppendResult(s->interp, caller, "() stale draw ref\n", (char *)NULL);
	return 0;
    }
    return 1;
}


static void
mged_edit_preview_name(struct bu_vls *name, const char *source_path)
{
    if (!name)
	return;
    bu_vls_sprintf(name, "%s%s", MGED_EDIT_PREVIEW_PREFIX,
	    source_path ? source_path : "");
}


static int
mged_edit_preview_view_exists(void *view_ctx, const char *name)
{
    return (view_ctx && name && name[0] &&
	    ged_draw_view_context_feature_exists(view_ctx, name)) ? 1 : 0;
}


static int
mged_edit_preview_exists_any(struct mged_state *s, const char *name)
{
    int checked = 0;

    if (!s || !name || !name[0])
	return 0;

    for (size_t di = 0; di < BU_PTBL_LEN(&active_dm_set); di++) {
	struct mged_dm *m_dmp = (struct mged_dm *)BU_PTBL_GET(&active_dm_set, di);
	void *view_ctx = (m_dmp && m_dmp->dm_view_state) ?
	    m_dmp->dm_view_state->vs_gvp : NULL;
	if (!view_ctx)
	    continue;
	checked = 1;
	if (mged_edit_preview_view_exists(view_ctx, name))
	    return 1;
    }

    if (!checked && view_state && view_state->vs_gvp)
	return mged_edit_preview_view_exists(view_state->vs_gvp, name);

    return 0;
}


static int
mged_edit_preview_publish_view(
    struct mged_state *s,
    struct _view_state *vsp,
    const char *name,
    const char *source_path,
    struct rt_db_internal *ip,
    const mat_t mat)
{
    void *view_ctx = vsp ? vsp->vs_gvp : NULL;
    if (!s || !view_ctx || !name || !name[0] || !ip)
	return 0;

    struct ged_draw_view_edit_transaction transaction =
	GED_DRAW_VIEW_EDIT_TRANSACTION_INIT;
    transaction.event = GED_DRAW_VIEW_EDIT_PREVIEW_UPDATE;
    transaction.feature_name = name;
    transaction.owner = (const void *)s;
    transaction.source_path = source_path;
    transaction.edit_intent_id = source_path;
    transaction.edit_intent_role = "primitive-edit";
    transaction.dbip = s->dbip;
    transaction.internal = ip;
    transaction.matrix = mat;
    transaction.ttol = &s->tol.ttol;
    transaction.tol = &s->tol.tol;
    if (color_scheme) {
	transaction.color_valid = 1;
	transaction.color[0] = color_scheme->cs_edit_info[0];
	transaction.color[1] = color_scheme->cs_edit_info[1];
	transaction.color[2] = color_scheme->cs_edit_info[2];
    }

    if (!ged_draw_view_context_edit_transaction_apply(view_ctx,
	    &transaction, NULL))
	return -1;

    mged_refresh_request_view(s, vsp, GED_VIEW_REFRESH_VIEW);
    return 1;
}


static int
mged_edit_preview_publish_all(
    struct mged_state *s,
    const char *name,
    const char *source_path,
    struct rt_db_internal *ip,
    const mat_t mat)
{
    int published = 0;
    int attempted = 0;
    int failed = 0;

    if (!s || !name || !name[0] || !ip)
	return 0;

    for (size_t di = 0; di < BU_PTBL_LEN(&active_dm_set); di++) {
	struct mged_dm *m_dmp = (struct mged_dm *)BU_PTBL_GET(&active_dm_set, di);
	if (!m_dmp || !m_dmp->dm_view_state || !m_dmp->dm_view_state->vs_gvp)
	    continue;
	attempted = 1;
	int ret = mged_edit_preview_publish_view(s, m_dmp->dm_view_state,
	    name, source_path, ip, mat);
	if (ret > 0)
	    published = 1;
	else if (ret < 0)
	    failed = 1;
    }

    if (!attempted && view_state && view_state->vs_gvp) {
	int ret = mged_edit_preview_publish_view(s, view_state, name,
	    source_path, ip, mat);
	if (ret > 0)
	    published = 1;
	else if (ret < 0)
	    failed = 1;
    }

    return failed ? -1 : published;
}


static int
mged_edit_preview_clear_view(struct mged_state *s,
			     struct _view_state *vsp,
			     const char *name)
{
    void *view_ctx = vsp ? vsp->vs_gvp : NULL;
    int removed = 0;

    if (!view_ctx || !name || !name[0])
	return 0;

    struct ged_draw_view_edit_transaction transaction =
	GED_DRAW_VIEW_EDIT_TRANSACTION_INIT;
    transaction.event = GED_DRAW_VIEW_EDIT_PREVIEW_CANCEL;
    transaction.feature_name = name;
    removed = ged_draw_view_context_edit_transaction_apply(view_ctx,
	    &transaction, NULL);
    if (removed)
	mged_refresh_request_view(s, vsp, GED_VIEW_REFRESH_VIEW);
    return removed;
}


static int
mged_edit_preview_clear_all(struct mged_state *s, const char *name)
{
    int removed = 0;
    int attempted = 0;

    if (!s || !name || !name[0])
	return 0;

    for (size_t di = 0; di < BU_PTBL_LEN(&active_dm_set); di++) {
	struct mged_dm *m_dmp = (struct mged_dm *)BU_PTBL_GET(&active_dm_set, di);
	if (!m_dmp || !m_dmp->dm_view_state || !m_dmp->dm_view_state->vs_gvp)
	    continue;
	attempted = 1;
	if (mged_edit_preview_clear_view(s, m_dmp->dm_view_state, name))
	    removed = 1;
    }

    if (!attempted && view_state && view_state->vs_gvp)
	removed = mged_edit_preview_clear_view(s, view_state, name);

    return removed;
}


/*
 * Given an existing solid structure that may have been subjected to
 * solid editing, recompute the vector list, etc., to make the solid
 * the same as it originally was.
 *
 * Returns -
 * -1 error
 * 0 OK
 */
int
replot_original_solid(struct mged_state *s, ged_draw_shape_ref ref)
{
    if (s->dbip == DBI_NULL)
	return 0;

    struct ged_draw_shape_record rec;
    if (!ged_draw_shape_record_get(s->gedp, ref, &rec))
	return 0;
    if (!rec.fullpath || rec.fullpath->fp_len <= 0)
	return 0;
    char *source_path = db_path_to_string(rec.fullpath);
    if (!source_path)
	return -1;

    struct bu_vls feature_name = BU_VLS_INIT_ZERO;
    mged_edit_preview_name(&feature_name, source_path);
    if (mged_edit_preview_clear_all(s, bu_vls_cstr(&feature_name)))
	(void)ged_draw_shape_ref_set_visible(s->gedp, ref, 1);

    bu_vls_free(&feature_name);
    bu_free(source_path, "MGED edit-preview source path");
    return 0;
}


/*
 * Given the solid structure of a solid that has already been drawn,
 * and a new database record and transform matrix,
 * create a new vector list for that solid, and substitute.
 * Used for solid editing mode.
 *
 * Returns -
 * -1 error
 * 0 OK
 */
int
replot_modified_solid(
	struct mged_state *s,
	ged_draw_shape_ref ref,
	struct rt_db_internal *ip,
	const mat_t mat)
{
    if (!mged_check_shape_ref(s, ref, "replot_modified_solid")) {
	Tcl_AppendResult(s->interp, "replot_modified_solid() sp==NULL?\n", (char *)NULL);
	return -1;
    }

    struct ged_draw_shape_record rec;
    if (!ged_draw_shape_record_get(s->gedp, ref, &rec) ||
	    !rec.fullpath || rec.fullpath->fp_len <= 0) {
	Tcl_AppendResult(s->interp, "replot_modified_solid() stale draw path\n", (char *)NULL);
	return -1;
    }
    char *source_path = db_path_to_string(rec.fullpath);
    if (!source_path) {
	Tcl_AppendResult(s->interp, "replot_modified_solid() unable to resolve draw path\n", (char *)NULL);
	return -1;
    }

    struct bu_vls feature_name = BU_VLS_INIT_ZERO;
    mged_edit_preview_name(&feature_name, source_path);
    const int existing_preview =
	mged_edit_preview_exists_any(s, bu_vls_cstr(&feature_name));
    if (!rec.visible && !existing_preview) {
	bu_vls_free(&feature_name);
	bu_free(source_path, "MGED edit-preview source path");
	return 0;
    }

    RT_CK_DB_INTERNAL(ip);

    s->tol.ttol.magic = BG_TESS_TOL_MAGIC;
    s->tol.ttol.abs = s->tol.abs_tol;
    s->tol.ttol.rel = s->tol.rel_tol;
    s->tol.ttol.norm = s->tol.nrm_tol;

    int preview_status = mged_edit_preview_publish_all(s,
	bu_vls_cstr(&feature_name), source_path, ip, mat);
    if (preview_status < 0) {
	if (rec.leaf_name)
	    Tcl_AppendResult(s->interp, rec.leaf_name, ": edit-preview plot failure\n",
		    (char *)NULL);
	bu_vls_free(&feature_name);
	bu_free(source_path, "MGED edit-preview source path");
	return -1;
    }

    if (preview_status > 0 && rec.visible)
	(void)ged_draw_shape_ref_set_visible(s->gedp, ref, 0);
    if (preview_status > 0)
	ged_draw_shape_set_highlighted(s->gedp, ref, 1);

    bu_vls_free(&feature_name);
    bu_free(source_path, "MGED edit-preview source path");
    return 0;
}

void
add_solid_record_path_to_result(
    Tcl_Interp *interp,
    const struct ged_draw_shape_record *rec)
{
    struct bu_vls str = BU_VLS_INIT_ZERO;
    if (!rec || !rec->fullpath)
	return;
    db_path_to_vls(&str, rec->fullpath);
    Tcl_AppendResult(interp, bu_vls_addr(&str), " ", NULL);
    bu_vls_free(&str);
}

int
redraw_visible_objects(struct mged_state *s)
{
    const char *av[1] = {"redraw"};
    int ret = ged_exec_redraw(s->gedp, 1, av);

    if (ret & BRLCAD_ERROR)
	return TCL_ERROR;

    return TCL_OK;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
