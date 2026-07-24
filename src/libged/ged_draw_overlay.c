/*                  G E D _ D R A W _ O V E R L A Y . C
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
/** @file ged_draw_overlay.c
 *
 * Typed GED overlay insertion helpers.
 */

#include "common.h"

#include "bu/log.h"

#include "ged.h"
#include "ged/draw.h"
#include "./ged_private.h"
#include "./ged_draw_private.h"


void
ged_draw_overlay_erase_name(struct ged *gedp, const char *name)
{
    struct ged_view_context *view_ctx = ged_draw_active_view_ctx(gedp);
    ged_draw_overlay_erase_name_context(gedp, view_ctx, name);
}


int
ged_draw_overlay_geometry_insert(struct ged *gedp, const char *name,
	const struct ged_draw_overlay_geometry *geometry,
	long int rgb, fastf_t transparency, int draw_mode, int csoltab,
	ged_draw_shape_ref *out)
{
    if (out)
	*out = GED_DRAW_SHAPE_REF_NULL;
    if (!gedp || !name || !geometry)
	return -1;
    struct ged_view_context *view_ctx = ged_draw_active_view_ctx(gedp);
    if (!view_ctx)
	return 0;

    struct db_i *dbip = gedp->dbip;
    if (dbip == DBI_NULL)
	return 0;

    if (db_lookup(dbip, name, LOOKUP_QUIET) != RT_DIR_NULL) {
	bu_log("ged_draw_overlay_geometry_insert(%s) would clobber existing database entry, "
		"ignored\n", name);
	return -1;
    }

    ged_draw_overlay_erase_name(gedp, name);

    unsigned char basecolor[3] = {
	(unsigned char)((rgb >> 16) & 0xFF),
	(unsigned char)((rgb >>  8) & 0xFF),
	(unsigned char)(rgb & 0xFF)
    };

    return ged_draw_overlay_geometry_insert_context(gedp, view_ctx, name,
	    geometry, basecolor, transparency, draw_mode, csoltab, out);
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
