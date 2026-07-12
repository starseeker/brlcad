/*                        D O Z O O M . C
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
/** @file mged/scene_refresh.c
 *
 * MGED retained scene refresh.
 */

#include "common.h"

#include "vmath.h"
#include "bn.h"
#include "dm/view.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "ged/view.h"

#include "./mged.h"
#include "./sedit.h"
#include "./mged_dm.h"

mat_t perspective_mat;

static void
_mged_obol_refresh(struct mged_state *s, void *view_ctx)
{
    if (!s || !s->gedp || !view_ctx ||
	    !ged_draw_obol_controller_opaque_for_view(view_ctx))
	return;

    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    txn.view = view_ctx;
    (void)ged_draw_apply_transaction(s->gedp, &txn, &result);
    ged_draw_transaction_result_free(&result);
}

void
mged_obol_scene_refresh(struct mged_state *s)
{
    void *view_ctx = view_state->vs_gvp;
    struct bv *view = mged_view_context_view(view_ctx);
    fastf_t view_perspective;
    point_t view_eye_pos;

    bv_refresh_drawn_count_set(view, 0);

    view_perspective = bv_perspective_get(view);
    bv_eye_pos_get(view_eye_pos, view);

    if (view_perspective >= SMALL_FASTF) {
	if (!EQUAL(view_eye_pos[Z], 1.0)) {
	    point_t l, h;
	    VSET(l, -1.0, -1.0, -1.0);
	    VSET(h, 1.0, 1.0, 200.0);
	    deering_persp_mat(perspective_mat, l, h, view_eye_pos);
	} else {
	    persp_mat(perspective_mat, view_perspective, 1.0, 0.01,
		1.0e10, 1.0);
	}
	bv_pmat_set(view, perspective_mat);
    }

    _mged_obol_refresh(s, view_ctx);
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
