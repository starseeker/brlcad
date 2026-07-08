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
/** @file mged/dozoom.c
 *
 * MGED scene refresh for normal and stereo eyes.
 *
 * Obol display managers refresh through the high-level GED draw transaction
 * path and let the outer MGED frame's dm_draw_end render the retained
 * Obol/Coin scene.  Older non-Obol display managers still use the legacy
 * immediate draw fallback until those display paths are retired.
 */

#include "common.h"

#include <math.h>
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
mat_t identity;

/* Assumed physical screen width in mm (stereo eye-separation calculation) */
#ifndef SCR_WIDTH_PHYS
#  define SCR_WIDTH_PHYS 330
#endif

/* This is a holding place for the current display managers default wireframe color */
unsigned char geometry_default_color[] = { 255, 0, 0 };

/* Count draw records whose drawn revision matches the view's current
 * current retained-view frame revision, i.e. shapes that were actually
 * painted in the frame just rendered.  Stored in the view refresh record for
 * mouse-pick sequencing. */
struct _mged_count_drawn_ctx {
    int *np;
    uint64_t frame_rev;
};
static int
_mged_count_drawn_cb(const struct ged_draw_shape_record *rec, void *userdata)
{
    struct _mged_count_drawn_ctx *ctx =
	(struct _mged_count_drawn_ctx *)userdata;
    if (rec && rec->drawn_revision == ctx->frame_rev)
	(*ctx->np)++;
    return 1; /* continue traversal */
}

static int
_mged_count_visible_cb(const struct ged_draw_shape_record *rec, void *userdata)
{
    int *np = (int *)userdata;
    if (rec && rec->visible)
	(*np)++;
    return 1;
}

static int
_mged_high_level_refresh(struct mged_state *s, void *view_ctx)
{
    if (!s || !s->gedp || !view_ctx || !DMP ||
	    !ged_draw_obol_controller_opaque_for_view(view_ctx))
	return 0;

    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    txn.view = view_ctx;
    int ret = ged_draw_apply_transaction(s->gedp, &txn, &result);
    ged_draw_transaction_result_free(&result);
    return ret >= 0 ? 1 : 0;
}

/*
 * Paint one eye of the scene.
 *
 * which_eye == 0  Normal (non-stereo) view.
 * which_eye == 1  Stereo right eye.
 * which_eye == 2  Stereo left eye.
 *
 * In the Obol case rendering is prepared by a GED redraw transaction; the
 * legacy fallback still calls dm_draw_objs().  The stereo case differs from
 * the non-stereo case only in that:
 *   - v->gv_pmat is overridden with a Deering eye-offset perspective for
 *     the duration of the call;
 *   - dm_loadmatrix(DMP, gv_model2view, which_eye) is called once before
 *     refresh so legacy GL display managers can select the correct stereo
 *     viewport/scissor region.
 */
void
dozoom(struct mged_state *s, int which_eye)
{
    void *view_ctx = view_state->vs_gvp;
    struct bv *view = mged_view_context_view(view_ctx);
    fastf_t view_perspective;
    mat_t model2view;
    point_t view_eye_pos;

    /*
     * The vectorThreshold stuff in libdm may turn the
     * Tcl-crank causing s->mged_curr_dm to change.
     */
    struct mged_dm *save_dm_list = s->mged_curr_dm;

    bv_refresh_drawn_count_set(view, 0);

    /* Keep the retained view's display manager in sync so that
     * refresh code can find the DM.  This must be done every frame
     * because set_curr_dm() (called from refresh()) updates
     * s->mged_curr_dm without updating the view's display pointer. */
    ged_view_context_display_manager_set(view_ctx, (void *)DMP);

    /* gv_pmat may be replaced for the stereo path; remember the original
     * so we can restore it before returning. */
    mat_t saved_pmat;
    bv_pmat_get(saved_pmat, view);

    view_perspective = bv_perspective_get(view);
    bv_eye_pos_get(view_eye_pos, view);
    bv_model2view_get(model2view, view);

    if (which_eye == 0) {
	/* ----- Non-stereo: keep gv_pmat in sync with the perspective state.
	 * The GED "perspective" command normally maintains this, but the
	 * shear-perspective (gv_eye_pos[Z] != 1.0) mode needs an explicit
	 * rebuild here. */
	if (view_perspective >= SMALL_FASTF) {
	    if (!EQUAL(view_eye_pos[Z], 1.0)) {
		point_t l, h;
		VSET(l, -1.0, -1.0, -1.0);
		VSET(h,  1.0,  1.0, 200.0);
		deering_persp_mat(perspective_mat, l, h, view_eye_pos);
		bv_pmat_set(view, perspective_mat);
	    } else {
		persp_mat(perspective_mat, view_perspective,
			  (fastf_t)1.0f, (fastf_t)0.01f,
			  (fastf_t)1.0e10f, (fastf_t)1.0f);
		bv_pmat_set(view, perspective_mat);
	    }
	}
    } else {
	/* ----- Stereo: install a Deering eye-offset perspective into
	 * v->gv_pmat so the display host can load it via dm_loadpmatrix().
	 *
	 * Stereo requires a non-zero gv_perspective: the eye-distance
	 * to_eye_scr below derives from it.  When mv_perspective_mode is
	 * enabled but gv_perspective happens to be 0, fall back to a
	 * sensible default so we don't divide by zero. */
	fastf_t persp = view_perspective;
	if (persp < SMALL_FASTF)
	    persp = 90.0;
	fastf_t to_eye_scr = 1 / tan(persp * DEG2RAD * 0.5);
	fastf_t eye_delta_scr = mged_variables->mv_eye_sep_dist * 0.5 / SCR_WIDTH_PHYS;
	point_t l, h, eye;
	VSET(l, -1.0, -1.0, -1.0);
	VSET(h,  1.0,  1.0, 200.0);
	VSET(eye, 0.0, 0.0, to_eye_scr);

	if (which_eye == 1) {
	    eye[X] = eye_delta_scr;
	    printf("d=%gscr, d=%gmm, delta=%gscr\n",
		   to_eye_scr, to_eye_scr * SCR_WIDTH_PHYS, eye_delta_scr);
	    VPRINT("l", l); VPRINT("h", h);
	} else {
	    eye[X] = -eye_delta_scr;
	}
	deering_persp_mat(perspective_mat, l, h, eye);
	bv_pmat_set(view, perspective_mat);

	/* Force the display host to apply the perspective matrix even when the
	 * retained view perspective was 0; legacy drawing gates the projection
	 * load on a non-zero perspective angle. */
	if (view_perspective < SMALL_FASTF)
	    bv_perspective_set(view, persp);

	/* Stereo viewport / scissor selection.  gl_loadMatrix() inspects
	 * which_eye (1 = right, 2 = left) and adjusts glViewport+glScissor
	 * accordingly. */
	dm_loadmatrix(DMP, model2view, which_eye);
    }

    int high_level_refresh = _mged_high_level_refresh(s, view_ctx);
    if (!high_level_refresh) {
	/* Legacy fallback for non-Obol display managers. */
	dm_draw_objs(view_ctx);
    }

    /* Restore gv_pmat (no-op for which_eye == 0). */
    bv_pmat_set(view, saved_pmat);

    /* Count drawn objects for usepen.c zone-based picking. */
    if (s->gedp && ged_draw_scene_available(s->gedp)) {
	int ndrawn = 0;
	if (high_level_refresh) {
	    ged_draw_foreach_shape_record(s->gedp, _mged_count_visible_cb,
		    &ndrawn);
	} else {
	    struct _mged_count_drawn_ctx ctx;
	    ctx.np = &ndrawn;
	    ctx.frame_rev = bv_frame_revision_get(view);
	    ged_draw_foreach_shape_record(s->gedp, _mged_count_drawn_cb, &ctx);
	}
	bv_refresh_drawn_count_set(view, ndrawn);
    }

    if (s->mged_curr_dm != save_dm_list) set_curr_dm(s, save_dm_list);
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
