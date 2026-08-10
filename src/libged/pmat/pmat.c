/*                         P M A T . C
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
/** @file libged/pmat.c
 *
 * The pmat command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "bv.h"

#include "../ged_private.h"


int
ged_pmat_core(struct ged *gedp, int argc, const char *argv[])
{
    mat_t pmat;

    GED_CHECK_DATABASE_OPEN(gedp, BRLCAD_ERROR);
    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* get the perspective matrix */
    if (argc == 1) {
	const struct bv *view = bv_context_view_const((const struct bv_context *)view_ctx);
	bv_pmat_get(pmat, view);
	bn_encode_mat(gedp->ged_result_str, pmat, 1);
	return BRLCAD_OK;
    } else if (argc == 2) {
	/* set perspective matrix */
	if (bn_decode_mat(pmat, argv[1]) != 16)
	    return BRLCAD_ERROR;

	struct bv *view = bv_context_view((struct bv_context *)view_ctx);
	bv_pmat_set(view, pmat);
	ged_view_context_update(view_ctx);

	return BRLCAD_OK;
    }

    bu_vls_printf(gedp->ged_result_str, "Usage: %s", argv[0]);
    return BRLCAD_ERROR;
}


#include "../include/plugin.h"

#define GED_PMAT_COMMANDS(X, XID) \
    X(pmat, ged_pmat_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_PMAT_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_pmat", 1, GED_PMAT_COMMANDS)

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
