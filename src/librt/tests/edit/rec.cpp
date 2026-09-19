/*                           R E C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file rec.cpp
 *
 * Runtime edit coverage for REC's constraint-preserving operations.
 */

#include "common.h"

#include "bu/app.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "raytrace.h"
#include "edit_test_view.h"
#include "rt/rt_ecmds.h"


static void
check_rec_constraints(const struct rt_tgc_internal *rec)
{
    if (!VNEAR_EQUAL(rec->a, rec->c, VUNITIZE_TOL) ||
	!VNEAR_EQUAL(rec->b, rec->d, VUNITIZE_TOL) ||
	!NEAR_ZERO(VDOT(rec->h, rec->a), VUNITIZE_TOL) ||
	!NEAR_ZERO(VDOT(rec->h, rec->b), VUNITIZE_TOL))
	bu_exit(1, "ERROR: REC constraints were not preserved\n");
}


int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);
    if (argc != 1)
	return BRLCAD_ERROR;

    struct db_i *dbip = db_open_inmem();
    if (dbip == DBI_NULL)
	bu_exit(1, "ERROR: unable to create database instance\n");
    struct rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_INMEM);
    struct rt_tgc_internal *rec;
    BU_ALLOC(rec, struct rt_tgc_internal);
    rec->magic = RT_TGC_INTERNAL_MAGIC;
    VSET(rec->v, 1.0, 2.0, 3.0);
    VSET(rec->h, 0.0, 0.0, 8.0);
    VSET(rec->a, 3.0, 0.0, 0.0);
    VSET(rec->b, 0.0, 2.0, 0.0);
    VMOVE(rec->c, rec->a);
    VMOVE(rec->d, rec->b);
    wdb_export(wdbp, "rec", rec, ID_REC, 1.0);

    struct directory *dp = db_lookup(dbip, "rec", LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
	bu_exit(1, "ERROR: unable to create REC\n");
    struct db_full_path path;
    db_full_path_init(&path);
    db_add_node_to_full_path(&path, dp);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct rt_edit_view view;
    rt_edit_test_view_init(&view);
    struct rt_edit *edit = rt_edit_create(&path, dbip, &tol, &view);
    if (!edit)
	bu_exit(1, "ERROR: unable to create REC editor\n");

    point_t origin = {1.0, 2.0, 3.0};
    vect_t height = {0.0, 0.0, 8.0};
    if (!rt_edit_test_point_value(edit, ECMD_REC_SET_V, origin) ||
	!rt_edit_test_point_value(edit, ECMD_REC_SET_H, height) ||
	!rt_edit_test_scalar_value(edit, ECMD_REC_SCALE_R1, 3.0) ||
	!rt_edit_test_scalar_value(edit, ECMD_REC_SCALE_R2, 2.0) ||
	!rt_edit_test_scalar_value(edit, ECMD_REC_SCALE_R, 2.5))
	bu_exit(1, "ERROR: REC descriptor readback failed\n");

    EDOBJ[dp->d_minor_type].ft_set_edit_mode(edit, ECMD_REC_SET_V);
    edit->e_inpara = 3;
    VSET(edit->e_para, -1.0, 4.0, 6.0);
    if (rt_edit_process_result(edit) != BRLCAD_OK)
	bu_exit(1, "ERROR: REC base edit failed\n");

    EDOBJ[dp->d_minor_type].ft_set_edit_mode(edit, ECMD_REC_SET_H);
    edit->e_inpara = 3;
    VSET(edit->e_para, 0.0, 6.0, 0.0);
    if (rt_edit_process_result(edit) != BRLCAD_OK)
	bu_exit(1, "ERROR: REC height edit failed\n");
    struct rt_tgc_internal *edited =
	(struct rt_tgc_internal *)edit->es_int.idb_ptr;
    check_rec_constraints(edited);

    const int radius_commands[] = {
	ECMD_REC_SCALE_R1, ECMD_REC_SCALE_R2, ECMD_REC_SCALE_R
    };
    const fastf_t radii[] = {4.0, 5.0, 6.0};
    for (size_t i = 0; i < sizeof(radius_commands) / sizeof(radius_commands[0]); i++) {
	EDOBJ[dp->d_minor_type].ft_set_edit_mode(edit, radius_commands[i]);
	edit->e_inpara = 1;
	edit->e_para[0] = radii[i];
	if (rt_edit_process_result(edit) != BRLCAD_OK)
	    bu_exit(1, "ERROR: REC radius edit failed\n");
	check_rec_constraints(edited);
    }

    const vect_t a_before = {V3ARGS(edited->a)};
	EDOBJ[dp->d_minor_type].ft_set_edit_mode(edit, ECMD_REC_SCALE_R1);
	edit->e_inpara = 1;
	edit->e_para[0] = 0.0;
	if (rt_edit_process_result(edit) != BRLCAD_ERROR ||
	!VNEAR_EQUAL(a_before, edited->a, VUNITIZE_TOL))
	bu_exit(1, "ERROR: invalid REC radius was not rejected without mutation\n");

    rt_edit_destroy(edit);
    db_free_full_path(&path);
    db_close(dbip);
    return BRLCAD_OK;
}
