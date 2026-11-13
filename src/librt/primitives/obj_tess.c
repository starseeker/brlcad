/*                    O B J _ T E S S . C
 * BRL-CAD
 *
 * Copyright (c) 2010-2026 United States Government as represented by
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

#include "common.h"


#include "bn.h"
#include "raytrace.h"


int
rt_obj_tess(struct nmgregion **r, struct model *m, struct rt_db_internal *ip, const struct bg_tess_tol *ttol, const struct bn_tol *tol)
{
    const struct rt_functab *ft;

    if (!r || !ip)
	return -1;

    if (*r) NMG_CK_REGION(*r);
    if (m) NMG_CK_MODEL(m);
    RT_CK_DB_INTERNAL(ip);
    if (ttol) BG_CK_TESS_TOL(ttol);
    if (tol) BN_CK_TOL(tol);

    if (ip->idb_minor_type < 0)
	return -2;

    /* idb_minor_type is a primitive ID only for BRL-CAD geometry.  Binary
     * uniform objects retain their storage subtype there, and those numeric
     * values overlap ordinary primitive IDs.  The importer-provided method
     * table is the authoritative runtime type and prevents dispatching a
     * binunif payload to an unrelated geometry tessellator. */
    ft = ip->idb_meth;
    if (!ft)
	return -3;
    if (!ft->ft_tessellate)
	return -4;

    return ft->ft_tessellate(r, m, ip, ttol, tol);
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
