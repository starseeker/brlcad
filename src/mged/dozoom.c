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
 * MGED scene rendering via the BSG draw path, for both non-stereo and stereo
 * eyes.
 *
 * Both paths delegate to dm_draw_objs(), which resolves the retained BSG
 * scene to render batches.  Highlighted edit-mode objects are rendered at
 * their edited position by staging the edit matrix on the retained view
 * before the draw call; render-item drawing picks this up and swaps the
 * modelview matrix for the duration of each such object.
 *
 * The stereo path (which_eye == 1 or 2) builds a Deering perspective matrix
 * with a left/right eye offset and installs it into v->gv_pmat; dm_draw_objs()
 * then loads it via dm_loadpmatrix() in the standard separated
 * projection/modelview pipeline.  Stereo viewport selection (which controls
 * the split-screen left/right scissor regions in gl_loadMatrix) is performed
 * by an explicit dm_loadmatrix() with which_eye != 0 before dm_draw_objs()
 * runs.  v->gv_pmat is saved and restored around the stereo render so the
 * non-stereo perspective state is not perturbed.
 */

#include "common.h"

#include <math.h>
#include "vmath.h"
#include "bn.h"
#include "dm/view.h"
#include "ged/bsg_ged_draw.h"
#include "rt/view_legacy_bsg.h"

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

/*
 * Paint one eye of the scene.
 *
 * which_eye == 0  Normal (non-stereo) view.
 * which_eye == 1  Stereo right eye.
 * which_eye == 2  Stereo left eye.
 *
 * In all cases rendering is performed by dm_draw_objs(), which uses the
 * retained scene when the GED draw scene is available.  The stereo case differs from
 * the non-stereo case only in that:
 *   - v->gv_pmat is overridden with a Deering eye-offset perspective for
 *     the duration of the call (saved/restored around dm_draw_objs);
 *   - dm_loadmatrix(DMP, gv_model2view, which_eye) is called once before
 *     dm_draw_objs() so that gl_loadMatrix() can select the correct stereo
 *     viewport/scissor region (the matrix upload there is harmless because
 *     dm_draw_objs() will re-upload model2view immediately after).
 */
void
dozoom(struct mged_state *s, int which_eye)
{
    struct bsg_view *v = view_state->vs_gvp;
    fastf_t view_perspective;
    mat_t model2view;
    point_t view_eye_pos;

    /*
     * The vectorThreshold stuff in libdm may turn the
     * Tcl-crank causing s->mged_curr_dm to change.
     */
    struct mged_dm *save_dm_list = s->mged_curr_dm;

    rt_view_context_refresh_drawn_count_set_bsg(v, 0);

    /* Keep the retained view's display manager in sync so that
     * dm_draw_objs() can find the DM.  This must be done every frame
     * because set_curr_dm() (called from refresh()) updates
     * s->mged_curr_dm without updating the view's display pointer. */
    rt_view_context_display_manager_set_bsg(v, (void *)DMP);

    /* gv_pmat may be replaced for the stereo path; remember the original
     * so we can restore it before returning. */
    mat_t saved_pmat;
    rt_view_context_pmat_from_bsg(saved_pmat, v);

    view_perspective = rt_view_context_perspective_from_bsg(v);
    rt_view_context_eye_pos_from_bsg(view_eye_pos, v);
    rt_view_context_model2view_from_bsg(model2view, v);

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
		rt_view_context_pmat_set_bsg(v, perspective_mat);
	    } else {
		persp_mat(perspective_mat, view_perspective,
			  (fastf_t)1.0f, (fastf_t)0.01f,
			  (fastf_t)1.0e10f, (fastf_t)1.0f);
		rt_view_context_pmat_set_bsg(v, perspective_mat);
	    }
	}
    } else {
	/* ----- Stereo: install a Deering eye-offset perspective into
	 * v->gv_pmat so dm_draw_objs() will load it via dm_loadpmatrix().
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
	rt_view_context_pmat_set_bsg(v, perspective_mat);

	/* Force dm_draw_objs() to apply the perspective matrix even when the
	 * retained view perspective was 0; dm_draw_objs() gates the projection
	 * load on a non-zero perspective angle. */
	if (view_perspective < SMALL_FASTF)
	    rt_view_context_perspective_set_bsg(v, persp);

	/* Stereo viewport / scissor selection.  gl_loadMatrix() inspects
	 * which_eye (1 = right, 2 = left) and adjusts glViewport+glScissor
	 * accordingly; the matrix upload itself is then redone by
	 * dm_draw_objs() with which_eye=0 (which is a no-op for viewport). */
	dm_loadmatrix(DMP, model2view, which_eye);
    }

    /* Expose the edit-mode matrix on the view so render-item drawing can use
     * it for highlighted edit objects without a second pass. */
    if (s->global_editing_state != ST_VIEW)
	rt_view_context_edit_matrix_set_bsg(v, view_state->vs_model2objview);
    else
	rt_view_context_edit_matrix_clear_bsg(v);

    /* dm_draw_objs() handles:
     *   - framebuffer overlay/underlay
     *   - dm_loadmatrix(gv_model2view)
     *   - dm_loadpmatrix(gv_pmat) for perspective (incl. stereo eye-offset)
     *   - retained-scene render via render batches without exposing the
     *     view's scene attachment
     *   - per-object edit matrix swap for highlighted edit objects
    */
    dm_draw_objs(v);

    /* Clear edit-mat pointer now that the frame is done. */
    rt_view_context_edit_matrix_clear_bsg(v);

    /* Restore gv_pmat (no-op for which_eye == 0). */
    rt_view_context_pmat_set_bsg(v, saved_pmat);

    /* Count drawn objects for usepen.c zone-based picking.  Each rendered
     * shape records the bsg_view frame revision when painted; comparing the
     * semantic draw record against the retained frame revision gives a
     * frame-by-frame "what got drawn" count without application graph
     * traversal. */
    if (s->gedp && ged_draw_scene_available(s->gedp)) {
	int ndrawn = 0;
	struct _mged_count_drawn_ctx ctx;
	ctx.np = &ndrawn;
	ctx.frame_rev = rt_view_context_frame_revision_from_bsg(v);
	ged_draw_foreach_shape_record(s->gedp, _mged_count_drawn_cb, &ctx);
	rt_view_context_refresh_drawn_count_set_bsg(v, ndrawn);
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
