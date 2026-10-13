/*                       F A C E D E F . C
 * BRL-CAD
 *
 * Copyright (c) 1986-2026 United States Government as represented by
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
/** @file libged/facedef.c
 *
 */

#include "common.h"
#include "rt/edit.h"
#include "rt/primitives/arb8.h"
#include "../ged_private.h"
#include "./ged_edit.h"

int
edarb_facedef(void *data, int argc, const char *argv[])
{
    struct ged *gedp = (struct ged *)data;
    if (argc < 5) {
        bu_vls_printf(gedp->ged_result_str, "Usage: edarb facedef arb face form parameters");
        return argc == 2 ? GED_HELP : BRLCAD_ERROR;
    }
    struct db_full_path path = DB_FULL_PATH_INIT_ZERO;
    if (db_string_to_path(&path, gedp->dbip, argv[2]) < 0 || !path.fp_len) {
        db_free_full_path(&path);
        return BRLCAD_ERROR;
    }
    struct bn_tol tolerance = BN_TOL_INIT_TOL;
    struct rt_edit *edit = rt_edit_create(&path, gedp->dbip, &tolerance, gedp->ged_gvp);
    if (!edit) {
        db_free_full_path(&path);
        return BRLCAD_ERROR;
    }
    const char **parameters = (const char **)bu_calloc(argc - 1, sizeof(char *), "facedef arguments");
    parameters[0] = "facedef";
    for (int i = 3; i < argc; i++)
        parameters[i - 2] = argv[i];
    int ret = rt_arb_facedef(edit, argc - 2, parameters);
    bu_vls_vlscat(gedp->ged_result_str, edit->log_str);
    bu_free(parameters, "facedef arguments");
    if (ret == BRLCAD_OK && rt_db_put_internal(DB_FULL_PATH_CUR_DIR(&path), gedp->dbip, &edit->es_int) < 0)
        ret = BRLCAD_ERROR;
    rt_edit_destroy(edit);
    db_free_full_path(&path);
    return ret;
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
