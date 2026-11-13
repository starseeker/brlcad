/*                    G E D _ D R A W _ M A T E R I A L . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
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
/** @file libged/ged_draw_material.c
 *
 * GED material/color revision bookkeeping.
 */

#include "common.h"

#include "ged.h"
#include "ged/draw.h"
#include "./ged_private.h"
#include "./ged_draw_private.h"

/**
 * Return the current mater-revision counter.  The counter is bumped by
 * ged_draw_bump_material_revision() whenever the material/color table changes.
 * The backend adapter uses this value to stamp one committed recolor batch.
 */
uint64_t
ged_draw_material_revision(const struct ged *gedp)
{
    if (!gedp)
	return 0;
    return gedp->i->ged_gdp->gd_mater_rev;
}


/**
 * Bump the mater-revision counter (B4 activation).
 *
 * Called by the semantic material-change reducer before it emits its
 * committed recolor delta.
 */
void
ged_draw_bump_material_revision(struct ged *gedp)
{
    if (!gedp)
	return;
    gedp->i->ged_gdp->gd_mater_rev++;
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
