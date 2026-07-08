/*                          V I E W . C
 * BRL-CAD
 *
 * Copyright (c) 2007-2026 United States Government as represented by
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
#include "bv.h"
#include "dm.h"

int
dm_view_context_width_get(const void *view_ctx)
{
    return bv_context_width_get((const struct bv_context *)view_ctx);
}

int
dm_view_context_height_get(const void *view_ctx)
{
    return bv_context_height_get((const struct bv_context *)view_ctx);
}

int
dm_view_context_dimensions_set(void *view_ctx, int width, int height)
{
    return bv_context_dimensions_set((struct bv_context *)view_ctx,
	    width, height);
}

int
dm_view_context_refresh_request(void *view_ctx, uint32_t flags)
{
    return bv_context_refresh_request((struct bv_context *)view_ctx,
	    flags);
}

void
dm_draw_arrow(struct dm *dmp, point_t A, point_t B, fastf_t tip_length, fastf_t tip_width, fastf_t sf)
{
    point_t points[16];
    point_t BmA;
    point_t offset;
    point_t perp1, perp2;
    point_t a_base;
    point_t a_pt1, a_pt2, a_pt3, a_pt4;

    VSUB2(BmA, B, A);

    VUNITIZE(BmA);
    VSCALE(offset, BmA, -tip_length * sf);

    bn_vec_perp(perp1, BmA);
    VUNITIZE(perp1);

    VCROSS(perp2, BmA, perp1);
    VUNITIZE(perp2);

    VSCALE(perp1, perp1, tip_width * sf);
    VSCALE(perp2, perp2, tip_width * sf);

    VADD2(a_base, B, offset);
    VADD2(a_pt1, a_base, perp1);
    VADD2(a_pt2, a_base, perp2);
    VSUB2(a_pt3, a_base, perp1);
    VSUB2(a_pt4, a_base, perp2);

    VMOVE(points[0], B);
    VMOVE(points[1], a_pt1);
    VMOVE(points[2], B);
    VMOVE(points[3], a_pt2);
    VMOVE(points[4], B);
    VMOVE(points[5], a_pt3);
    VMOVE(points[6], B);
    VMOVE(points[7], a_pt4);
    VMOVE(points[8], a_pt1);
    VMOVE(points[9], a_pt2);
    VMOVE(points[10], a_pt2);
    VMOVE(points[11], a_pt3);
    VMOVE(points[12], a_pt3);
    VMOVE(points[13], a_pt4);
    VMOVE(points[14], a_pt4);
    VMOVE(points[15], a_pt1);

    (void)dm_draw_lines_3d(dmp, 16, points, 0);
}

void
dm_draw_faceplate(void *UNUSED(view_ctx))
{
}

void
dm_draw_objs(void *UNUSED(view_ctx))
{
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
