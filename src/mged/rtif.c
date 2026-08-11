/*                          R T I F . C
 * BRL-CAD
 *
 * Copyright (c) 1988-2026 United States Government as represented by
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
/** @file mged/rtif.c
 *
 * Routines to interface to RT, and RT-style command files
 *
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>		/* For struct timeval */
#endif
#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif
#include <sys/stat.h>
#include "bresource.h"

#include "tcl.h"

#include "vmath.h"
#include "raytrace.h"
#include "ged/scene.h"
#include "ged/view.h"

#include "./sedit.h"
#include "./mged.h"
#include "./mged_display.h"
#include "./cmd.h"


/* Callback: find the displayed shape for database solid "EYE" for rmats. */
struct _rtif_eye_data {
    struct directory *dp;
    Tcl_Interp *interp;
    vect_t sav_center;
    struct db_full_path *path;
    int found;
};

static int
_rtif_eye_shape_cb(const struct ged_scene_occurrence_info *rec, void *ud)
{
    struct _rtif_eye_data *d = (struct _rtif_eye_data *)ud;
    if (!rec || !rec->fullpath || rec->fullpath->fp_len <= 0) return 1;
    if (DB_FULL_PATH_CUR_DIR(rec->fullpath) != d->dp) return 1;
    VMOVE(d->sav_center, rec->center);
    db_dup_full_path(d->path, rec->fullpath);
    Tcl_AppendResult(d->interp, "animating EYE solid\n", (char *)NULL);
    d->found = 1;
    return 0; /* stop visiting */
}

static void
_rtif_use_current_view(struct mged_state *s)
{
    if (!s || !s->gedp || !view_state || !view_state->vs_gvp)
	return;

    ged_view_active_ctx_set(s->gedp, view_state->vs_gvp);
}

static void
_rtif_request_current_view_refresh(struct mged_state *s)
{
    mged_refresh_request_view(s, view_state, GED_VIEW_REFRESH_VIEW);
    if (s)
	mged_display_repaint_request(s->mged_curr_display, MGED_REPAINT_INTERACTION);
}

/**
 * rt, rtarea, rtweight, rtcheck, and rtedge all use this.
 */
int
cmd_rt(ClientData clientData,
       Tcl_Interp *interp,
       int argc,
       const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;

    int ret;
    Tcl_DString ds;

    CHECK_DBI_NULL;

    /* skip past _mged_ */
    if (argv[0][0] == '_' && argv[0][1] == 'm' &&
	bu_strncmp(argv[0], "_mged_", 6) == 0)
	argv[0] += 6;

    Tcl_DStringInit(&ds);

    ret = ged_exec(s->gedp, argc, (const char **)argv);

    Tcl_DStringAppend(&ds, bu_vls_addr(s->gedp->ged_result_str), -1);
    Tcl_DStringResult(interp, &ds);

    if (ret == BRLCAD_OK)
	return TCL_OK;

    return TCL_ERROR;
}


/**
 * Invoke any program with the current view & stuff, just like
 * an "rt" command (above).
 * Typically used to invoke a remote RT (hence the name).
 */
int
cmd_rrt(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;
    int ret;
    Tcl_DString ds;

    CHECK_DBI_NULL;

    if (argc < 2) {
	struct bu_vls vls = BU_VLS_INIT_ZERO;

	bu_vls_printf(&vls, "help rrt");
	Tcl_Eval(interp, bu_vls_addr(&vls));
	bu_vls_free(&vls);
	return TCL_ERROR;
    }

    if (not_state(s, ST_VIEW, "Ray-trace of current view"))
	return TCL_ERROR;

    Tcl_DStringInit(&ds);

    ret = ged_exec(s->gedp, argc, (const char **)argv);
    Tcl_DStringAppend(&ds, bu_vls_addr(s->gedp->ged_result_str), -1);
    Tcl_DStringResult(interp, &ds);

    if (ret == BRLCAD_OK)
	return TCL_OK;

    return TCL_ERROR;
}


/**
 * Read in one view in the old RT format.
 */
static int
rt_read(FILE *fp, fastf_t *scale, fastf_t *eye, fastf_t *mat)
{
    int i;
    double d;

    if (fscanf(fp, "%lf", &d) != 1) return -1;
    *scale = d*0.5;
    if (fscanf(fp, "%lf", &d) != 1) return -1;
    eye[X] = d;
    if (fscanf(fp, "%lf", &d) != 1) return -1;
    eye[Y] = d;
    if (fscanf(fp, "%lf", &d) != 1) return -1;
    eye[Z] = d;
    for (i=0; i < 16; i++) {
	if (fscanf(fp, "%lf", &d) != 1)
	    return -1;
	mat[i] = d;
    }
    return 0;
}


/**
 * Load view matrices from a file.  rmats filename [mode]
 *
 * Modes:
 * -1 put eye in viewcenter (default)
 * 0 put eye in viewcenter, don't rotate.
 * 1 leave view alone, animate solid named "EYE"
 */
int
f_rmats(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;
    FILE *fp = NULL;
    fastf_t scale = 0.0;
    mat_t rot;
    struct directory *dp = NULL;
    vect_t eye_model = VINIT_ZERO;
    vect_t sav_center = VINIT_ZERO;
    vect_t xlate = VINIT_ZERO;
    struct db_full_path eye_path;
    struct rt_db_internal eye_internal;
    mat_t eye_path_mat;
    char *eye_path_name = NULL;
    int eye_preview_active = 0;

    /* static due to setjmp */
    static int mode = 0;

    db_full_path_init(&eye_path);
    RT_DB_INTERNAL_INIT(&eye_internal);
    MAT_IDN(eye_path_mat);

    CHECK_DBI_NULL;

    if (argc < 2 || 3 < argc) {
	struct bu_vls vls = BU_VLS_INIT_ZERO;

	bu_vls_printf(&vls, "help rmats");
	Tcl_Eval(interp, bu_vls_addr(&vls));
	bu_vls_free(&vls);
	return TCL_ERROR;
    }

    if (not_state(s, ST_VIEW, "animate from matrix file"))
	return TCL_ERROR;

    if ((fp = fopen(argv[1], "r")) == NULL) {
	perror(argv[1]);
	return TCL_ERROR;
    }

    mode = -1;
    if (argc > 2)
	mode = atoi(argv[2]);
    switch (mode) {
	case 1:
	    if ((dp = db_lookup(s->dbip, "EYE", LOOKUP_NOISY)) == RT_DIR_NULL) {
		mode = -1;
		break;
	    }

	    {
		struct _rtif_eye_data d;
		d.dp = dp;
		d.interp = interp;
		VSETALL(d.sav_center, 0.0);
		d.path = &eye_path;
		d.found = 0;
		ged_scene_occurrences_visit(s->gedp, _rtif_eye_shape_cb, &d);
		if (d.found) {
		    VMOVE(sav_center, d.sav_center);
		    if (rt_db_get_internal(&eye_internal, dp, s->dbip,
			    NULL) < 0) {
			db_free_full_path(&eye_path);
			fclose(fp);
			return TCL_ERROR;
		    }
		    (void)db_path_to_mat(s->dbip, &eye_path, eye_path_mat,
			eye_path.fp_len - 1);
		    eye_path_name = db_path_to_string(&eye_path);
		    if (!eye_path_name) {
			rt_db_free_internal(&eye_internal);
			db_free_full_path(&eye_path);
			fclose(fp);
			return TCL_ERROR;
		    }
		    struct ged_scene_path_request request;
		    ged_scene_path_request_init(&request);
		    request.path = eye_path_name;
		    request.match = GED_SCENE_PATH_MATCH_EXACT;
		    (void)ged_scene_visibility_set(s->gedp, &request, 0, NULL);
		    eye_preview_active = 1;
		    goto work;
		}
	    }
	    /* Fall through */
	default:
	case -1:
	    mode = -1;
	    Tcl_AppendResult(interp, "default mode:  eyepoint at (0, 0, 1) viewspace\n", (char *)NULL);
	    break;
	case 0:
	    Tcl_AppendResult(interp, "rotation suppressed, center is eyepoint\n", (char *)NULL);
	    break;
    }
work:
    /* FIXME: this isn't portable or seem well thought-out */
    if (setjmp(jmp_env) == 0)
	(void)signal(SIGINT, sig3);  /* allow interrupts */
    else
	return TCL_OK;

    struct bv *view = mged_view_context_view(view_state->vs_gvp);
    while (!feof(fp) &&
	   rt_read(fp, &scale, eye_model, rot) >= 0) {
	switch (mode) {
	    case -1:
		/* First step:  put eye in center */
		bv_scale_set(view, scale);
		bv_rotation_set(view, rot);
		bv_center_set(view, eye_model);
		new_mats(s);
		/* Second step:  put eye in front */
		VSET(xlate, 0.0, 0.0, -1.0);	/* correction factor */
		mat_t view2model;
		bv_view2model_get(view2model, view);
		MAT4X3PNT(eye_model, view2model, xlate);
		bv_center_set(view, eye_model);
		new_mats(s);
		break;
	    case 0: {
		mat_t top_view;
		MAT_IDN(top_view);
		bv_scale_set(view, scale);
		bv_rotation_set(view, top_view);	/* top view */
		bv_center_set(view, eye_model);
		new_mats(s);
		break;
	    }
	    case 1:
		if (eye_preview_active) {
		    mat_t translation;
		    mat_t preview_mat;
		    MAT_IDN(translation);
		    VSUB2(xlate, eye_model, sav_center);
		    MAT_DELTAS_VEC(translation, xlate);
		    bn_mat_mul(preview_mat, translation, eye_path_mat);
		    struct ged_view_edit_transaction transaction =
			GED_VIEW_EDIT_TRANSACTION_INIT;
		    transaction.event = GED_VIEW_EDIT_PREVIEW_UPDATE;
		    transaction.feature_name = "_mged_rmats_eye";
		    transaction.owner = (const void *)s;
		    transaction.source_path = eye_path_name;
		    transaction.edit_intent_id = "rmats-eye";
		    transaction.edit_intent_role = "animation";
		    transaction.dbip = s->dbip;
		    transaction.internal = &eye_internal;
		    transaction.matrix = preview_mat;
		    transaction.ttol = &s->tol.ttol;
		    transaction.tol = &s->tol.tol;
		    (void)ged_view_edit_transaction_apply_all(s->gedp,
			&transaction);
		}
		break;
	}
	mged_refresh_request_view(s, view_state, GED_VIEW_REFRESH_VIEW);
	refresh(s);	/* Draw new display */
    }

    if (eye_preview_active) {
	struct ged_view_edit_transaction transaction =
	    GED_VIEW_EDIT_TRANSACTION_INIT;
	transaction.event = GED_VIEW_EDIT_PREVIEW_CANCEL;
	transaction.feature_name = "_mged_rmats_eye";
	transaction.owner = (const void *)s;
	(void)ged_view_edit_transaction_apply_all(s->gedp, &transaction);
	struct ged_scene_path_request request;
	ged_scene_path_request_init(&request);
	request.path = eye_path_name;
	request.match = GED_SCENE_PATH_MATCH_EXACT;
	(void)ged_scene_visibility_set(s->gedp, &request, 1, NULL);
    }

    if (eye_path_name)
	bu_free(eye_path_name, "rmats EYE path");
    rt_db_free_internal(&eye_internal);
    db_free_full_path(&eye_path);

    fclose(fp);
    (void)mged_svbase(s);

    (void)signal(SIGINT, SIG_IGN);
    return TCL_OK;
}


/**
 * Invoke nirt with the current view & stuff
 */
int
f_nirt(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;
    int ret;
    Tcl_DString ds;

    CHECK_DBI_NULL;

    /* skip past _mged_ */
    if (argv[0][0] == '_' && argv[0][1] == 'm' &&
	bu_strncmp(argv[0], "_mged_", 6) == 0)
	argv[0] += 6;

    Tcl_DStringInit(&ds);
    _rtif_use_current_view(s);

    if (mged_variables->mv_use_air) {
	int insertArgc = 2;
	char *insertArgv[3];
	int newArgc;
	char **newArgv;

	insertArgv[0] = "-u";
	insertArgv[1] = "1";
	insertArgv[2] = (char *)0;
	newArgv = bu_argv_dupinsert(1, insertArgc, (const char **)insertArgv, argc, (const char **)argv);
	newArgc = argc + insertArgc;
	ret = ged_exec(s->gedp, newArgc, (const char **)newArgv);
	bu_argv_free(newArgc, newArgv);
    } else {
	ret = ged_exec(s->gedp, argc, (const char **)argv);
    }

    Tcl_DStringAppend(&ds, bu_vls_addr(s->gedp->ged_result_str), -1);
    Tcl_DStringResult(interp, &ds);

    if (ret == BRLCAD_OK) {
	_rtif_request_current_view_refresh(s);
	return TCL_OK;
    }

    return TCL_ERROR;
}


int
f_vnirt(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;
    int ret;
    Tcl_DString ds;

    CHECK_DBI_NULL;

    /* skip past _mged_ */
    if (argv[0][0] == '_' && argv[0][1] == 'm' &&
	bu_strncmp(argv[0], "_mged_", 6) == 0)
	argv[0] += 6;

    Tcl_DStringInit(&ds);
    _rtif_use_current_view(s);

    ret = ged_exec(s->gedp, argc, (const char **)argv);

    Tcl_DStringAppend(&ds, bu_vls_addr(s->gedp->ged_result_str), -1);
    Tcl_DStringResult(interp, &ds);

    if (ret == BRLCAD_OK) {
	_rtif_request_current_view_refresh(s);
	return TCL_OK;
    }

    return TCL_ERROR;
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
