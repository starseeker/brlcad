/*                B S G _ G E D _ D R A W _ M A T E R I A L . C
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
/** @file libged/b_s_g___g_e_d___d_r_a_w___m_a_t_e_r_i_a_l_._c.c
 *
 * GED material/color revision bookkeeping for Obol database sources.
 */

#include "common.h"

#include "ged.h"
#include "ged/draw.h"
#include "./ged_private.h"
#include "./bsg_ged_draw_private.h"

struct _color_ctx {
    struct ged  *gedp;
    struct db_i *dbip;
    uint64_t     mater_rev;
};


static int
_color_shape_record_cb(const struct ged_draw_shape_record *rec, void *ud)
{
    struct _color_ctx *ctx = (struct _color_ctx *)ud;
    if (!ctx || !rec)
	return 1;

    /* Color-table interpretation lives in libbrlobol.  GED only requests
     * refresh for drawn Obol database sources at the current material
     * revision. */
    (void)ged_draw_shape_ref_refresh_material_color(ctx->gedp, rec->ref,
	    ctx->dbip, ctx->mater_rev);
    return 1;
}


static void
_sg_color_soltab(struct ged *gedp)
{
    if (ged_draw_obol_database_sources_refresh_material_colors(gedp,
	    gedp->dbip, gedp->i->ged_gdp->gd_mater_rev))
	return;

    struct _color_ctx ctx;
    ctx.gedp      = gedp;
    ctx.dbip      = gedp->dbip;
    ctx.mater_rev = gedp->i->ged_gdp->gd_mater_rev;
    ged_draw_foreach_shape_record(gedp, _color_shape_record_cb, &ctx);
    /* B4 activation: do NOT bump gd_mater_rev here.  The counter is
     * event-driven — only ged_draw_bump_material_revision() (called by
     * material-change commands) advances it. */
}


void
ged_draw_refresh_material_colors(struct ged *gedp)
{
    if (!gedp)
	return;
    _sg_color_soltab(gedp);
}


/**
 * Return the current mater-revision counter.  The counter is bumped by
 * ged_draw_bump_material_revision() whenever the material/color table changes.
 * color_from_soltab() does NOT bump the counter; it only stamps per-shape
 * s_color_rev fields to match the current counter value.
 *
 * Consumers that cache per-solid colors can store a snapshot of this
 * value and skip re-querying when the counter is unchanged.  For example:
 *
 *   if (saved_mater_rev != ged_draw_material_revision(gedp)) {
 *       ged_draw_refresh_material_colors(gedp);
 *       saved_mater_rev = ged_draw_material_revision(gedp);
 *   }
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
 * Must be called after any operation that changes the effective material
 * or color table so that the next ged_draw_refresh_material_colors() call
 * recolors shapes whose s_color_rev is now stale.
 *
 * Typical callers: 'color', 'mater', 'rmater', 'edmater' commands and
 * any other code path that mutates dbip->dbi_mater or per-combination
 * shader/rgb attributes.
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
