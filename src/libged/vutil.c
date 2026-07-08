/*                         V U T I L . C
 * BRL-CAD
 *
 * Copyright (c) 2008-2026 United States Government as represented by
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
/** @file libged/vutil.c
 *
 * This file contains view related utility functions.
 *
 */

#include "common.h"

#include <string.h>

#include "bv.h"
#include "bg/line_layer.h"
#include "nmg/display.h"
#include "./ged_private.h"
#include "ged/view.h"
#include "ged/draw.h"

int
_ged_do_rot(struct ged *gedp,
	    char coord,
	    mat_t rmat,
	    int (*func)(struct ged *, char, char, mat_t))
{
    mat_t temp1, temp2;
    void *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx)
	return BRLCAD_ERROR;
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    char rotate_about = bv_rotate_about_get(view);

    if (func != (int (*)(struct ged *, char, char, mat_t))0)
	return (*func)(gedp, coord, rotate_about, rmat);

    switch (coord) {
	case 'm': {
	    /* transform model rotations into view rotations */
	    mat_t view_rotation;
	    bv_rotation_get(view_rotation, view);
	    bn_mat_inv(temp1, view_rotation);
	    bn_mat_mul(temp2, view_rotation, rmat);
	    bn_mat_mul(rmat, temp2, temp1);
	    break;
	}
	case 'v':
	default:
	    break;
    }

    /* Calculate new view center */
    if (rotate_about != 'v') {
	point_t rot_pt;
	point_t new_origin;
	mat_t viewchg, viewchginv;
	mat_t model2view;
	mat_t view2model;
	point_t new_cent_view;
	point_t new_cent_model;

	bv_model2view_get(model2view, view);

	switch (rotate_about) {
	    case 'e':
		VSET(rot_pt, 0.0, 0.0, 1.0);
		break;
	    case 'k': {
		point_t keypoint;
		bv_keypoint_get(keypoint, view);
		MAT4X3PNT(rot_pt, model2view, keypoint);
		break;
	    }
	    case 'm':
		/* rotate around model center (0, 0, 0) */
		VSET(new_origin, 0.0, 0.0, 0.0);
		MAT4X3PNT(rot_pt, model2view, new_origin);
		break;
	    default:
		return BRLCAD_ERROR;
	}

	bn_mat_xform_about_pnt(viewchg, rmat, rot_pt);
	bn_mat_inv(viewchginv, viewchg);

	/* Convert origin in new (viewchg) coords back to old view coords */
	VSET(new_origin, 0.0, 0.0, 0.0);
	MAT4X3PNT(new_cent_view, viewchginv, new_origin);
	bv_view2model_get(view2model, view);
	MAT4X3PNT(new_cent_model, view2model, new_cent_view);
	bv_center_set(view, new_cent_model);
    }

    /* pure rotation */
    {
	mat_t view_rotation;
	bv_rotation_get(view_rotation, view);
	bn_mat_mul2(rmat, view_rotation);
	bv_rotation_set(view, view_rotation);
    }
    ged_view_context_update(view_ctx);

    return BRLCAD_OK;
}


int
_ged_do_slew(struct ged *gedp, vect_t svec)
{
    point_t model_center;
    mat_t view2model;
    void *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx)
	return BRLCAD_ERROR;
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);

    bv_view2model_get(view2model, view);
    MAT4X3PNT(model_center, view2model, svec);
    bv_center_set(view, model_center);
    ged_view_context_update(view_ctx);

    return BRLCAD_OK;
}


int
_ged_do_tra(struct ged *gedp,
	    char coord,
	    vect_t tvec,
	    int (*func)(struct ged *, char, vect_t))
{
    point_t delta;
    point_t work;
    point_t vc, nvc;
    void *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx)
	return BRLCAD_ERROR;
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);

    if (func != (int (*)(struct ged *, char, vect_t))0)
	return (*func)(gedp, coord, tvec);

    switch (coord) {
	case 'm':
	    VSCALE(delta, tvec, -gedp->dbip->dbi_base2local);
	    {
		bv_center_get(vc, view);
	    }
	    break;
	case 'v':
	default:
	    VSCALE(tvec, tvec, -2.0 * gedp->dbip->dbi_base2local *
		    (1.0 / bv_size_get(view)));
	    {
		mat_t view2model;
		bv_view2model_get(view2model, view);
		MAT4X3PNT(work, view2model, tvec);
		bv_center_get(vc, view);
	    }
	    VSUB2(delta, work, vc);
	    break;
    }

    VSUB2(nvc, vc, delta);
    bv_center_set(view, nvc);
    ged_view_context_update(view_ctx);

    return BRLCAD_OK;
}

void
nmg_plot_eu(struct ged *gedp, struct edgeuse *es_eu, const struct bn_tol *tol)
{
    if (!gedp || !es_eu || !tol)
	return;

    if (*es_eu->g.magic_p != NMG_EDGE_G_LSEG_MAGIC)
	return;

    struct model *m = nmg_find_model(&es_eu->l.magic);
    NMG_CK_MODEL(m);

    /* get space for list of items processed */
    long *tab = (long *)bu_calloc(m->maxindex+1, sizeof(long), "nmg_ed tab[]");
    struct bg_line_layer_builder *plot = bg_line_layer_builder_create();

    nmg_line_layer_around_eu(plot, es_eu, tab, 1, tol);
    (void)_ged_line_layer_builder_publish_command_scene_feature(gedp,
	    "nmg::_EU_", plot, "nmg", "command-result", "nmg::_EU_",
	    "nmg-edgeuse", 0);

    bg_line_layer_builder_free(plot);
    bu_free((void *)tab, "nmg_ed tab[]");
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
