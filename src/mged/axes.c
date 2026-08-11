/*                          A X E S . C
 * BRL-CAD
 *
 * Copyright (c) 1998-2026 United States Government as represented by
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
/** @file mged/axes.c
 *
 */

#include "common.h"

#include <math.h>
#include <string.h>

#include "vmath.h"
#include "bn.h"
#include "ged.h"
#include "ged/view.h"
#include "ged/view_feature_batch.h"
#include "rt/view.h"

#include "./mged.h"
#include "./mged_display.h"

/* local sp_hook function */
static void ax_set_dirty_flag(const struct bu_structparse *, const char *, void *, const char *, void *);

static int
mged_hud_axes_replace(struct ged_view_context *view_ctx, const char *name,
	const struct bv_axes_state *axes, const mat_t rotation)
{
    struct ged_view_feature_batch_desc desc = GED_VIEW_FEATURE_BATCH_DESC_INIT;
    desc.owner_id = "mged-edit-axes";
    desc.owner_role = "faceplate";
    desc.local = 1;
    struct ged_view_feature_batch *batch =
	ged_view_feature_batch_begin(view_ctx, &desc);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_hud_axes_replace(batch, name, axes,
	    rotation)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

struct _axes_state default_axes_state = {
    /* ax_rc */			1,
    /* ax_model_draw */    	0,
    /* ax_model_size */		500,
    /* ax_model_linewidth */	1,
    /* ax_model_pos */		VINIT_ZERO,
    /* ax_view_draw */    	0,
    /* ax_view_size */		500,
    /* ax_view_linewidth */	1,
    /* ax_view_pos */		{ 0, 0 },
    /* ax_edit_draw */		0,
    /* ax_edit_size1 */		500,
    /* ax_edit_size2 */		500,
    /* ax_edit_linewidth1 */	1,
    /* ax_edit_linewidth2 */	1
};


#define AX_O(_m) bu_offsetof(struct _axes_state, _m)
struct bu_structparse axes_vparse[] = {
    {"%d", 1, "model_draw",	AX_O(ax_model_draw),		ax_set_dirty_flag, NULL, NULL },
    {"%d", 1, "model_size",	AX_O(ax_model_size),		ax_set_dirty_flag, NULL, NULL },
    {"%d", 1, "model_linewidth",AX_O(ax_model_linewidth),	ax_set_dirty_flag, NULL, NULL },
    {"%f", 3, "model_pos",	AX_O(ax_model_pos),		ax_set_dirty_flag, NULL, NULL },
    {"%d", 1, "view_draw",	AX_O(ax_view_draw),		ax_set_dirty_flag, NULL, NULL },
    {"%d", 1, "view_size",	AX_O(ax_view_size),		ax_set_dirty_flag, NULL, NULL },
    {"%d", 1, "view_linewidth",	AX_O(ax_view_linewidth),	ax_set_dirty_flag, NULL, NULL },
    {"%d", 2, "view_pos",	AX_O(ax_view_pos),		ax_set_dirty_flag, NULL, NULL },
    {"%d", 1, "edit_draw",	AX_O(ax_edit_draw),		ax_set_dirty_flag, NULL, NULL },
    {"%d", 1, "edit_size1",	AX_O(ax_edit_size1),		ax_set_dirty_flag, NULL, NULL },
    {"%d", 1, "edit_size2",	AX_O(ax_edit_size2),		ax_set_dirty_flag, NULL, NULL },
    {"%d", 1, "edit_linewidth1",AX_O(ax_edit_linewidth1),	ax_set_dirty_flag, NULL, NULL },
    {"%d", 1, "edit_linewidth2",AX_O(ax_edit_linewidth2),	ax_set_dirty_flag, NULL, NULL },
    {"",   0, (char *)0,	0,			 BU_STRUCTPARSE_FUNC_NULL, NULL, NULL }
};


static void
ax_set_dirty_flag(const struct bu_structparse *UNUSED(sdp),
		  const char *UNUSED(name),
		  void *base,
		  const char *UNUSED(value),
		  void *data)
{
    struct mged_state *s = (struct mged_state *)data;
    struct _axes_state *changed_state = (struct _axes_state *)base;
    MGED_CK_STATE(s);
    for (size_t i = 0; i < BU_PTBL_LEN(&active_display_set); i++) {
	struct mged_display *m_dmp = (struct mged_display *)BU_PTBL_GET(&active_display_set, i);
	if (m_dmp->display_axes_state == changed_state) {
	    m_dmp->display_axes_state_dirty = 1;
	    mged_display_repaint_request(m_dmp, MGED_REPAINT_DEVICE_SETTING);
	}
    }
}


void
mged_edit_axes_state_sync(struct mged_state *s)
{
    point_t v_ap1;                 /* axes position in view coordinates */
    point_t v_ap2;                 /* axes position in view coordinates */
    mat_t model2view;
    mat_t rot_mat;
    mat_t view_rotation;
    struct bv_axes_state gas;
    struct ged_view_context *view_ctx = view_state->vs_gvp;

    bv_model2view_get(model2view, mged_view_context_view_const(view_ctx));
    bv_rotation_get(view_rotation, mged_view_context_view_const(view_ctx));

    if (s->global_editing_state == ST_S_EDIT) {
	MAT4X3PNT(v_ap1, model2view, MEDIT(s)->e_axes_pos);
	MAT4X3PNT(v_ap2, model2view, MEDIT(s)->curr_e_axes_pos);
    } else if (s->global_editing_state == ST_O_EDIT) {
	point_t m_ap2;

	MAT4X3PNT(v_ap1, model2view, MEDIT(s)->e_keypoint);
	MAT4X3PNT(m_ap2, MEDIT(s)->model_changes, MEDIT(s)->e_keypoint);
	MAT4X3PNT(v_ap2, model2view, m_ap2);
    } else {
	(void)mged_hud_axes_replace(view_ctx,
		"_faceplate/edit_axes/initial", NULL, NULL);
	(void)mged_hud_axes_replace(view_ctx,
		"_faceplate/edit_axes/current", NULL, NULL);
	return;
    }

    memset(&gas, 0, sizeof(gas));
    gas.draw = axes_state->ax_edit_draw;
    gas.label_flag = 1;
    VMOVE(gas.axes_pos, v_ap1);
    gas.axes_size = axes_state->ax_edit_size1 * RT_INV_VIEW;
    VMOVE(gas.axes_color, color_scheme->cs_edit_axes1);
    VMOVE(gas.label_color, color_scheme->cs_edit_axes_label1);
    gas.line_width = axes_state->ax_edit_linewidth1;

    (void)mged_hud_axes_replace(view_ctx,
	    "_faceplate/edit_axes/initial", &gas, view_rotation);

    memset(&gas, 0, sizeof(gas));
    gas.draw = axes_state->ax_edit_draw;
    gas.label_flag = 1;
    VMOVE(gas.axes_pos, v_ap2);
    gas.axes_size = axes_state->ax_edit_size2 * RT_INV_VIEW;
    VMOVE(gas.axes_color, color_scheme->cs_edit_axes2);
    VMOVE(gas.label_color, color_scheme->cs_edit_axes_label2);
    gas.line_width = axes_state->ax_edit_linewidth2;

    bn_mat_mul(rot_mat, view_rotation, MEDIT(s)->acc_rot_sol);
    (void)mged_hud_axes_replace(view_ctx,
	    "_faceplate/edit_axes/current", &gas, rot_mat);
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
