/*                      C H G M O D E L . C
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
/** @file mged/chgmodel.c
 *
 * This module contains functions which change particulars of the
 * model, generally on a single solid or combination.
 * Changes to the tree structure of the model are done in chgtree.c
 *
 */

#include "common.h"
#include "ged/edit.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "vmath.h"
#include "bu/units.h"
#include "bn.h"
#include "raytrace.h"
#include "nmg.h"
#include "rt/geom.h"
#include "ged.h"
#include "wdb.h"

#include "./mged.h"
#include "./mged_dm.h"
#include "./sedit.h"
#include "./cmd.h"
#include "./f_cmd.h"


/* defined in chgview.c */
extern int edit_com(struct mged_state *s, int argc, const char *argv[]);

/* defined in buttons.c */
extern int be_s_trans(ClientData, Tcl_Interp *, int, char **);

/*
 * Create a new solid of a given type
 * (Generic, or explicit)
 */
int
f_make(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;
    Tcl_DString ds;
    int ret;

    CHECK_DBI_NULL;
    CHECK_READ_ONLY;

    if (argc == 3) {
	const char *av[8];
	char center[512];
	char scale[128];

	sprintf(center, "%.17f %.17f %.17f",
		(ZERO(view_state->vs_gvp->gv_center[MDX])) ? 0.0 : -view_state->vs_gvp->gv_center[MDX],
		(ZERO(view_state->vs_gvp->gv_center[MDY])) ? 0.0 : -view_state->vs_gvp->gv_center[MDY],
		(ZERO(view_state->vs_gvp->gv_center[MDZ])) ? 0.0 : -view_state->vs_gvp->gv_center[MDZ]);
	sprintf(scale, "%.17f", view_state->vs_gvp->gv_scale * 2.0);

	av[0] = argv[0];
	av[1] = "-o";
	av[2] = center;
	av[3] = "-s";
	av[4] = scale;
	av[5] = argv[1];
	av[6] = argv[2];
	av[7] = (char *)0;

	ret = ged_exec(s->gedp, 7, (const char **)av);
    } else
	ret = ged_exec(s->gedp, argc, (const char **)argv);

    Tcl_DStringInit(&ds);
    Tcl_DStringAppend(&ds, bu_vls_addr(s->gedp->ged_result_str), -1);
    Tcl_DStringResult(interp, &ds);

    if (ret == BRLCAD_OK) {
	const char *av[4];

	av[0] = "draw";
	av[1] = "-R";
	av[2] = argv[argc-2];
	av[3] = NULL;
	edit_com(s, 3, av);
    } else if (ret != GED_HELP) {
	/* not OK and not help */
	return TCL_ERROR;
    }

    return TCL_OK;
}


int
mged_rot_obj(struct mged_state *s, int iflag, fastf_t *argvect)
{
    point_t model_pt;
    point_t point;
    point_t s_point;
    mat_t temp;
    vect_t v_work;

    s->update_views = 1;
    dm_set_dirty(DMP, 1);

    if (movedir != ROTARROW) {
	/* NOT in object rotate mode - put it in obj rot */
	movedir = ROTARROW;
    }

    /* find point for rotation to take place wrt */
    VMOVE(model_pt, MEDIT(s)->e_keypoint);

    MAT4X3PNT(point, MEDIT(s)->model_changes, model_pt);

    /* Find absolute translation vector to go from "model_pt" to
     * "point" without any of the rotations in "MEDIT(s)->model_changes"
     */
    VSCALE(s_point, point, MEDIT(s)->model_changes[15]);
    VSUB2(v_work, s_point, model_pt);

    /* REDO "MEDIT(s)->model_changes" such that:
     * 1. NO rotations (identity)
     * 2. trans == v_work
     * 3. same scale factor
     */
    MAT_IDN(temp);
    MAT_DELTAS_VEC(temp, v_work);
    temp[15] = MEDIT(s)->model_changes[15];
    MAT_COPY(MEDIT(s)->model_changes, temp);

    /* build new rotation matrix */
    MAT_IDN(temp);
    bn_mat_angles(temp, argvect[0], argvect[1], argvect[2]);

    if (iflag) {
	/* apply accumulated rotations */
	bn_mat_mul2(MEDIT(s)->acc_rot_sol, temp);
    }

    /*XXX*/ MAT_COPY(MEDIT(s)->acc_rot_sol, temp); /* used to rotate solid/object axis */

    /* Record the new rotation matrix into the revised
     * MEDIT(s)->model_changes matrix wrt "point"
     */
    wrt_point(MEDIT(s)->model_changes, temp, MEDIT(s)->model_changes, point);

    new_edit_mats(s);

    return TCL_OK;
}


static int
mged_object_transform(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;
    bu_vls_trunc(s->gedp->ged_result_str, 0);
    int ret = ged_edit_session_command(s->gedp, argc, argv);
    Tcl_AppendResult(interp, bu_vls_cstr(s->gedp->ged_result_str), (char *)NULL);
    if (ret != BRLCAD_OK)
        return TCL_ERROR;
    new_edit_mats(s);
    return TCL_OK;
}

int
f_rot_obj(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    return mged_object_transform(clientData, interp, argc, argv);
}

int
f_sc_obj(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    return mged_object_transform(clientData, interp, argc, argv);
}

int
f_tr_obj(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    return mged_object_transform(clientData, interp, argc, argv);
}

int
f_qorot(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    return mged_object_transform(clientData, interp, argc, argv);
}


void
set_localunit_TclVar(struct mged_state *s)
{
    struct bu_vls vls = BU_VLS_INIT_ZERO;
    struct bu_vls units_vls = BU_VLS_INIT_ZERO;
    const char *str = NULL;

    if (s->dbip == DBI_NULL)
	return;

    str = bu_units_string(s->dbip->dbi_local2base);
    if (str)
	bu_vls_strcpy(&units_vls, str);
    else
	bu_vls_printf(&units_vls, "%gmm", s->dbip->dbi_local2base);

    bu_vls_strcpy(&vls, "localunit");
    Tcl_SetVar(s->interp, bu_vls_addr(&vls), bu_vls_addr(&units_vls), TCL_GLOBAL_ONLY);

    bu_vls_free(&vls);
    bu_vls_free(&units_vls);
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
