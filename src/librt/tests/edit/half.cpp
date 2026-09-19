/*                         H A L F . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file half.cpp
 *
 * Runtime edit and readback coverage for the halfspace descriptor.
 */

#include "common.h"

#include "bu/app.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "raytrace.h"
#include "edit_test_view.h"
#include "rt/rt_ecmds.h"


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
    struct rt_half_internal *half;
    BU_ALLOC(half, struct rt_half_internal);
    half->magic = RT_HALF_INTERNAL_MAGIC;
    VSET(half->eqn, 0.0, 0.0, 1.0);
    half->eqn[W] = -4.0;
    wdb_export(wdbp, "half", half, ID_HALF, 1.0);

    struct directory *dp = db_lookup(dbip, "half", LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
	bu_exit(1, "ERROR: unable to create halfspace\n");
    struct db_full_path path;
    db_full_path_init(&path);
    db_add_node_to_full_path(&path, dp);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct rt_edit_view view;
    rt_edit_test_view_init(&view);
    struct rt_edit *edit = rt_edit_create(&path, dbip, &tol, &view);
    if (!edit)
	bu_exit(1, "ERROR: unable to create halfspace editor\n");

    if (!rt_edit_test_scalar_value(edit, ECMD_HALF_SET_D, -4.0))
	bu_exit(1, "ERROR: halfspace descriptor readback failed\n");

    EDOBJ[dp->d_minor_type].ft_set_edit_mode(edit, ECMD_HALF_SET_D);
    edit->e_inpara = 1;
    edit->e_para[0] = -9.0;
    if (rt_edit_process_result(edit) != BRLCAD_OK)
	bu_exit(1, "ERROR: halfspace edit failed\n");
    struct rt_half_internal *edited =
	(struct rt_half_internal *)edit->es_int.idb_ptr;
    if (!NEAR_EQUAL(edited->eqn[W], -9.0, VUNITIZE_TOL))
	bu_exit(1, "ERROR: expected D=-9, got %g\n", edited->eqn[W]);
    if (!rt_edit_test_scalar_value(edit, ECMD_HALF_SET_D, -9.0))
	bu_exit(1, "ERROR: halfspace updated readback failed\n");

    rt_edit_destroy(edit);
    db_free_full_path(&path);
    db_close(dbip);
    return BRLCAD_OK;
}
