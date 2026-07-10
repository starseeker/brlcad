/*                        V L I S T . C
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
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

#include <stdio.h>

#include "vmath.h"
#include "bu/list.h"
#include "bu/log.h"
#include "bg/plot3.h"
#include "bsg/vlist.h"


size_t
bsg_vlist_cmd_cnt(bsg_vlist *vlist)
{
    return bg_vlist_cmd_cnt(vlist);
}


int
bsg_vlist_bbox(struct bu_list *vlistp, point_t *bmin, point_t *bmax, size_t *length, int *dispmode)
{
    return bg_vlist_bbox(vlistp, bmin, bmax, length, dispmode);
}


const char *
bsg_vlist_get_cmd_description(int cmd)
{
    return bg_vlist_get_cmd_description(cmd);
}


size_t
bsg_ck_vlist(const struct bu_list *vhead)
{
    return bg_ck_vlist(vhead);
}


void
bsg_vlist_copy(struct bu_list *vlists, struct bu_list *dest, const struct bu_list *src)
{
    bg_vlist_copy(vlists, dest, src);
}


void
bsg_vlist_cleanup(struct bu_list *hd)
{
    bg_vlist_cleanup(hd);
}


void
bsg_vlist_export(struct bu_vls *vls, struct bu_list *hp, const char *name)
{
    bg_vlist_export(vls, hp, name);
}


void
bsg_vlist_import(struct bu_list *vlists, struct bu_list *hp, struct bu_vls *namevls, const unsigned char *buf)
{
    bg_vlist_import(vlists, hp, namevls, buf);
}


struct bsg_vlblock *
bsg_vlblock_init(struct bu_list *free_vlist_hd, int max_ent)
{
    return bg_vlblock_init(free_vlist_hd, max_ent);
}


void
bsg_vlblock_free(struct bsg_vlblock *vbp)
{
    bg_vlblock_free(vbp);
}


struct bu_list *
bsg_vlblock_find(struct bsg_vlblock *vbp, int r, int g, int b)
{
    return bg_vlblock_find(vbp, r, g, b);
}


void
bsg_vlist_rpp(struct bu_list *vlists, struct bu_list *hd, const point_t minn, const point_t maxx)
{
    bg_vlist_rpp(vlists, hd, minn, maxx);
}


void
bsg_vlist_arb8(struct bu_list *vhead, struct bu_list *vlfree, point_t pts[8])
{
    bg_vlist_arb8(vhead, vlfree, pts);
}


void
bsg_plot_vlblock(FILE *fp, const struct bsg_vlblock *vbp)
{
    bg_plot_vlblock(fp, vbp);
}


void
bsg_vlist_to_uplot(FILE *fp, const struct bu_list *vhead)
{
    bg_vlist_to_uplot(fp, vhead);
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
