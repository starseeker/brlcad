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

#include <math.h>
#include <string.h>
#include "bnetwork.h"

#include "vmath.h"
#include "bn/mat.h"
#include "bg/plot3.h"
#include "bg/vlist.h"
#include "bu/cv.h"
#include "bu/list.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "bu/str.h"


size_t
bg_vlist_cmd_cnt(bg_vlist *vlist)
{
    size_t num_commands;
    bg_vlist *vp;

    if (UNLIKELY(vlist == NULL)) {
	return 0;
    }

    num_commands = 0;
    for (BU_LIST_FOR(vp, bg_vlist, &(vlist->l))) {
	num_commands += vp->nused;
    }

    return num_commands;
}


static int
bg_vlist_bbox_internal(bg_vlist *vp, point_t *bmin, point_t *bmax, int *disp_mode, int *dispmode_used)
{
    size_t i;
    size_t nused = vp->nused;
    int *cmd = vp->cmd;
    point_t *pt = vp->pt;

    for (i = 0; i < nused; i++, cmd++, pt++) {
	if (*disp_mode == 1 && *cmd != BG_VLIST_MODEL_MAT)
	    continue;
	*disp_mode = 0;
	switch (*cmd) {
	    case BG_VLIST_POLY_START:
	    case BG_VLIST_POLY_VERTNORM:
	    case BG_VLIST_TRI_START:
	    case BG_VLIST_TRI_VERTNORM:
	    case BG_VLIST_POINT_SIZE:
	    case BG_VLIST_LINE_WIDTH:
	    case BG_VLIST_MODEL_MAT:
		break;
	    case BG_VLIST_LINE_MOVE:
	    case BG_VLIST_LINE_DRAW:
	    case BG_VLIST_POLY_MOVE:
	    case BG_VLIST_POLY_DRAW:
	    case BG_VLIST_POLY_END:
	    case BG_VLIST_TRI_MOVE:
	    case BG_VLIST_TRI_DRAW:
	    case BG_VLIST_TRI_END:
		V_MIN((*bmin)[X], (*pt)[X]);
		V_MAX((*bmax)[X], (*pt)[X]);
		V_MIN((*bmin)[Y], (*pt)[Y]);
		V_MAX((*bmax)[Y], (*pt)[Y]);
		V_MIN((*bmin)[Z], (*pt)[Z]);
		V_MAX((*bmax)[Z], (*pt)[Z]);
		break;
	    case BG_VLIST_DISPLAY_MAT:
		*disp_mode = 1;
		*dispmode_used = 1;
		/* fall through */
	    case BG_VLIST_POINT_DRAW:
		V_MIN((*bmin)[X], (*pt)[X]-1.0);
		V_MAX((*bmax)[X], (*pt)[X]+1.0);
		V_MIN((*bmin)[Y], (*pt)[Y]-1.0);
		V_MAX((*bmax)[Y], (*pt)[Y]+1.0);
		V_MIN((*bmin)[Z], (*pt)[Z]-1.0);
		V_MAX((*bmax)[Z], (*pt)[Z]+1.0);
		break;
	    default:
		return *cmd;
	}
    }

    return 0;
}


int
bg_vlist_bbox(struct bu_list *vlistp, point_t *bmin, point_t *bmax, size_t *length, int *dispmode)
{
    bg_vlist *vp;
    int cmd = 0;
    int disp_mode = 0;
    int dispmode_used = 0;
    size_t len = 0;

    for (BU_LIST_FOR(vp, bg_vlist, vlistp)) {
	cmd = bg_vlist_bbox_internal(vp, bmin, bmax, &disp_mode, &dispmode_used);
	if (cmd)
	    break;
	len += vp->nused;
    }

    if (length)
	*length = len;
    if (dispmode)
	*dispmode = dispmode_used;
    return cmd;
}


const char *
bg_vlist_get_cmd_description(int cmd)
{
    static const char *bg_vlist_cmd_descriptions[] = {
	"line move",
	"line draw",
	"poly start",
	"poly move",
	"poly draw",
	"poly end",
	"poly vnorm",
	"tri start",
	"tri move",
	"tri draw",
	"tri end",
	"tri vnorm",
	"point draw",
	"point size",
	"line width",
	"display mat",
	"model mat"
    };

    if (cmd >= 0 && cmd <= BG_VLIST_CMD_MAX)
	return bg_vlist_cmd_descriptions[cmd];

    return "**unknown*";
}


size_t
bg_ck_vlist(const struct bu_list *vhead)
{
    bg_vlist *vp;
    size_t npts = 0;

    for (BU_LIST_FOR(vp, bg_vlist, vhead)) {
	size_t i;
	size_t nused = vp->nused;
	int *cmd = vp->cmd;
	point_t *pt = vp->pt;

	BG_CK_VLIST(vp);
	npts += nused;

	for (i = 0; i < nused; i++, cmd++, pt++) {
	    int j;

	    for (j = 0; j < 3; j++) {
		if ((*pt)[j] > -INFINITY && (*pt)[j] < INFINITY) {
		    /* Number is good. */
		} else {
		    bu_log("  %s (%g, %g, %g)\n",
			   bg_vlist_get_cmd_description(*cmd),
			   V3ARGS(*pt));
		    bu_bomb("bg_ck_vlist() bad coordinate value\n");
		}
		if (*cmd < 0 || *cmd > BG_VLIST_CMD_MAX) {
		    bu_log("cmd = x%x (%d.)\n", *cmd, *cmd);
		    bu_bomb("bg_ck_vlist() bad vlist command\n");
		}
	    }
	}
    }
    return npts;
}


void
bg_vlist_copy(struct bu_list *vlists, struct bu_list *dest, const struct bu_list *src)
{
    bg_vlist *vp;

    for (BU_LIST_FOR(vp, bg_vlist, src)) {
	size_t i;
	size_t nused = vp->nused;
	int *cmd = vp->cmd;
	point_t *pt = vp->pt;
	for (i = 0; i < nused; i++, cmd++, pt++) {
	    BG_ADD_VLIST(vlists, dest, *pt, *cmd);
	}
    }
}


void
bg_vlist_cleanup(struct bu_list *hd)
{
    bg_vlist *vp;

    if (!BU_LIST_IS_INITIALIZED(hd)) {
	BU_LIST_INIT(hd);
	return;
    }

    while (BU_LIST_WHILE(vp, bg_vlist, hd)) {
	BG_CK_VLIST(vp);
	BU_LIST_DEQUEUE(&(vp->l));
	bu_free((char *)vp, "bg_vlist");
    }
}


void
bg_vlist_export(struct bu_vls *vls, struct bu_list *hp, const char *name)
{
    bg_vlist *vp;
    size_t nelem;
    size_t namelen;
    size_t nbytes;
    unsigned char *buf;
    unsigned char *bp;

    BU_CK_VLS(vls);

    nelem = 0;
    for (BU_LIST_FOR(vp, bg_vlist, hp)) {
	nelem += vp->nused;
    }

    namelen = strlen(name) + 1;
    nbytes = namelen + 4 + nelem * (1 + ELEMENTS_PER_VECT * SIZEOF_NETWORK_DOUBLE) + 2;

    bu_vls_extend(vls, nbytes + 1);
    bu_vls_setlen(vls, (int)nbytes);
    buf = (unsigned char *)bu_vls_addr(vls);
    *(uint32_t *)buf = htonl((uint32_t)nelem);
    bp = buf + sizeof(uint32_t);
    bu_strlcpy((char *)bp, name, namelen);
    bp += namelen;

    for (BU_LIST_FOR(vp, bg_vlist, hp)) {
	size_t i;
	size_t nused = vp->nused;
	int *cmd = vp->cmd;
	for (i = 0; i < nused; i++) {
	    *bp++ = *cmd++;
	}
    }

    for (BU_LIST_FOR(vp, bg_vlist, hp)) {
	size_t i;
	size_t nused = vp->nused;
	point_t *pt = vp->pt;
	double point[ELEMENTS_PER_POINT];

	for (i = 0; i < nused; i++) {
	    VMOVE(point, pt[i]);
	    bu_cv_htond(bp, (unsigned char *)point, ELEMENTS_PER_VECT);
	    bp += ELEMENTS_PER_VECT * SIZEOF_NETWORK_DOUBLE;
	}
    }
    *bp = '\0';
}


void
bg_vlist_import(struct bu_list *vlists, struct bu_list *hp, struct bu_vls *namevls, const unsigned char *buf)
{
    const unsigned char *bp;
    const unsigned char *pp;
    size_t nelem;
    size_t namelen;
    size_t i;
    double point[ELEMENTS_PER_POINT];

    BU_CK_VLS(namevls);

    nelem = ntohl(*(uint32_t *)buf);
    bp = buf + 4;

    namelen = strlen((char *)bp) + 1;
    bu_vls_strncpy(namevls, (char *)bp, namelen);
    bp += namelen;

    pp = bp + nelem;

    for (i = 0; i < nelem; i++) {
	int cmd = *bp++;
	bu_cv_ntohd((unsigned char *)point, pp, ELEMENTS_PER_POINT);
	pp += ELEMENTS_PER_POINT * SIZEOF_NETWORK_DOUBLE;
	BG_ADD_VLIST(vlists, hp, point, cmd);
    }
}


struct bg_vlblock *
bg_vlblock_init(struct bu_list *free_vlist_hd, int max_ent)
{
    struct bg_vlblock *vbp;
    size_t i;

    if (!BU_LIST_IS_INITIALIZED(free_vlist_hd))
	BU_LIST_INIT(free_vlist_hd);

    BU_ALLOC(vbp, struct bg_vlblock);
    vbp->magic = BG_VLBLOCK_MAGIC;
    vbp->free_vlist_hd = free_vlist_hd;
    vbp->max = max_ent;
    vbp->head = (struct bu_list *)bu_calloc(vbp->max,
					    sizeof(struct bu_list), "head[]");
    vbp->rgb = (long *)bu_calloc(vbp->max,
				 sizeof(long), "rgb[]");

    for (i = 0; i < vbp->max; i++) {
	vbp->rgb[i] = 0;
	BU_LIST_INIT(&(vbp->head[i]));
    }

    vbp->rgb[0] = 0xFFFF00L;
    vbp->rgb[1] = 0xFFFFFFL;
    vbp->nused = 2;

    return vbp;
}


void
bg_vlblock_free(struct bg_vlblock *vbp)
{
    size_t i;

    BG_CK_VLBLOCK(vbp);
    for (i = 0; i < vbp->nused; i++) {
	if (vbp->rgb[i] == 0)
	    continue;
	if (BU_LIST_IS_EMPTY(&(vbp->head[i])))
	    continue;
	BG_FREE_VLIST(vbp->free_vlist_hd, &(vbp->head[i]));
    }

    bu_free((char *)(vbp->head), "head[]");
    bu_free((char *)(vbp->rgb), "rgb[]");
    bu_free((char *)vbp, "bg_vlblock");
}


struct bu_list *
bg_vlblock_find(struct bg_vlblock *vbp, int r, int g, int b)
{
    long newrgb;
    size_t n;
    size_t omax;

    BG_CK_VLBLOCK(vbp);

    newrgb = ((r & 0xFF) << 16) | ((g & 0xFF) << 8) | (b & 0xFF);

    for (n = 0; n < vbp->nused; n++) {
	if (vbp->rgb[n] == newrgb)
	    return &(vbp->head[n]);
    }
    if (vbp->nused < vbp->max) {
	n = vbp->nused++;
	vbp->rgb[n] = newrgb;
	return &(vbp->head[n]);
    }

    omax = vbp->max;
    vbp->max *= 2;

    for (n = 0; n < omax; n++)
	if (BU_LIST_IS_EMPTY(&vbp->head[n]))
	    vbp->head[n].forw = BU_LIST_NULL;

    vbp->head = (struct bu_list *)bu_realloc((void *)vbp->head,
					     vbp->max * sizeof(struct bu_list),
					     "head[]");
    vbp->rgb = (long *)bu_realloc((void *)vbp->rgb,
				  vbp->max * sizeof(long),
				  "rgb[]");

    for (n = 0; n < omax; n++) {
	if (vbp->head[n].forw == BU_LIST_NULL) {
	    vbp->head[n].forw = &vbp->head[n];
	    vbp->head[n].back = &vbp->head[n];
	} else {
	    vbp->head[n].forw->back = &vbp->head[n];
	    vbp->head[n].back->forw = &vbp->head[n];
	}
    }

    for (n = omax; n < vbp->max; n++) {
	vbp->rgb[n] = 0;
	BU_LIST_INIT(&vbp->head[n]);
    }

    return bg_vlblock_find(vbp, r, g, b);
}


void
bg_vlist_rpp(struct bu_list *vlists, struct bu_list *hd, const point_t minn, const point_t maxx)
{
    point_t p;

    VSET(p, minn[X], minn[Y], minn[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_MOVE);

    VSET(p, minn[X], maxx[Y], minn[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);
    VSET(p, minn[X], maxx[Y], maxx[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);
    VSET(p, minn[X], minn[Y], maxx[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);
    VSET(p, minn[X], minn[Y], minn[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);

    VSET(p, maxx[X], minn[Y], minn[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);

    VSET(p, maxx[X], maxx[Y], minn[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);
    VSET(p, maxx[X], maxx[Y], maxx[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);
    VSET(p, maxx[X], minn[Y], maxx[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);
    VSET(p, maxx[X], minn[Y], minn[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);

    VSET(p, minn[X], maxx[Y], minn[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_MOVE);
    VSET(p, maxx[X], maxx[Y], minn[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);

    VSET(p, minn[X], minn[Y], maxx[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_MOVE);
    VSET(p, maxx[X], minn[Y], maxx[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);

    VSET(p, minn[X], maxx[Y], maxx[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_MOVE);
    VSET(p, maxx[X], maxx[Y], maxx[Z]);
    BG_ADD_VLIST(vlists, hd, p, BG_VLIST_LINE_DRAW);
}


void
bg_vlist_arb8(struct bu_list *vhead, struct bu_list *vlfree, point_t pts[8])
{
    if (!vhead || !vlfree || !pts)
	return;

    /* Bottom face loop: MOVE(0) DRAW(1) DRAW(2) DRAW(3) DRAW(0)  - 5 cmds */
    BG_ADD_VLIST(vlfree, vhead, pts[0], BG_VLIST_LINE_MOVE);
    BG_ADD_VLIST(vlfree, vhead, pts[1], BG_VLIST_LINE_DRAW);
    BG_ADD_VLIST(vlfree, vhead, pts[2], BG_VLIST_LINE_DRAW);
    BG_ADD_VLIST(vlfree, vhead, pts[3], BG_VLIST_LINE_DRAW);
    BG_ADD_VLIST(vlfree, vhead, pts[0], BG_VLIST_LINE_DRAW);

    /* Top face loop: MOVE(4) DRAW(5) DRAW(6) DRAW(7) DRAW(4)     - 5 cmds */
    BG_ADD_VLIST(vlfree, vhead, pts[4], BG_VLIST_LINE_MOVE);
    BG_ADD_VLIST(vlfree, vhead, pts[5], BG_VLIST_LINE_DRAW);
    BG_ADD_VLIST(vlfree, vhead, pts[6], BG_VLIST_LINE_DRAW);
    BG_ADD_VLIST(vlfree, vhead, pts[7], BG_VLIST_LINE_DRAW);
    BG_ADD_VLIST(vlfree, vhead, pts[4], BG_VLIST_LINE_DRAW);

    /* Four vertical edges: 4 x (MOVE + DRAW)                     - 8 cmds */
    BG_ADD_VLIST(vlfree, vhead, pts[0], BG_VLIST_LINE_MOVE);
    BG_ADD_VLIST(vlfree, vhead, pts[4], BG_VLIST_LINE_DRAW);
    BG_ADD_VLIST(vlfree, vhead, pts[1], BG_VLIST_LINE_MOVE);
    BG_ADD_VLIST(vlfree, vhead, pts[5], BG_VLIST_LINE_DRAW);
    BG_ADD_VLIST(vlfree, vhead, pts[2], BG_VLIST_LINE_MOVE);
    BG_ADD_VLIST(vlfree, vhead, pts[6], BG_VLIST_LINE_DRAW);
    BG_ADD_VLIST(vlfree, vhead, pts[3], BG_VLIST_LINE_MOVE);
    BG_ADD_VLIST(vlfree, vhead, pts[7], BG_VLIST_LINE_DRAW);

    /* Total: 5 + 5 + 8 = 18 vlist commands */
}


void
bg_vlist_3string(struct bu_list *vhead,
		 struct bu_list *free_hd,
		 const char *string,
		 const point_t origin,
		 const mat_t rot,
		 double scale)
{
    register unsigned char *cp;
    double offset;
    int ysign;
    vect_t temp;
    vect_t loc;
    mat_t xlate_to_origin;
    mat_t mat;

    if (!string || !*string)
	return;

    MAT_IDN(xlate_to_origin);
    MAT_DELTAS_VEC(xlate_to_origin, origin);
    bn_mat_mul(mat, xlate_to_origin, rot);

    offset = 0;
    for (cp = (unsigned char *)string; *cp; cp++, offset += scale) {
	register int *p;
	register int stroke;

	VSET(temp, offset, 0, 0);
	MAT4X3PNT(loc, mat, temp);
	BG_ADD_VLIST(free_hd, vhead, loc, BG_VLIST_LINE_MOVE);

	for (p = plot3_font_getchar(cp); ((stroke = *p)) != PLOT3_FONT_LAST; p++) {
	    int draw;

	    if (stroke == PLOT3_FONT_NEGY) {
		ysign = -1;
		stroke = *++p;
	    } else {
		ysign = 1;
	    }

	    if (stroke < 0) {
		stroke = -stroke;
		draw = 0;
	    } else {
		draw = 1;
	    }

	    VSET(temp,
		 (stroke / 11) * 0.1 * scale + offset,
		 (ysign * (stroke % 11)) * 0.1 * scale,
		 0);
	    MAT4X3PNT(loc, mat, temp);
	    BG_ADD_VLIST(free_hd, vhead, loc,
		    draw ? BG_VLIST_LINE_DRAW : BG_VLIST_LINE_MOVE);
	}
    }
}


void
bg_vlist_2string(struct bu_list *vhead,
		 struct bu_list *free_hd,
		 const char *string,
		 double x,
		 double y,
		 double scale,
		 double theta)
{
    mat_t mat;
    vect_t p;

    bn_mat_angles(mat, 0.0, 0.0, theta);
    VSET(p, x, y, 0);
    bg_vlist_3string(vhead, free_hd, string, p, mat, scale);
}


void
bg_vlist_to_uplot(FILE *fp, const struct bu_list *vhead)
{
    bg_vlist *vp;

    for (BU_LIST_FOR(vp, bg_vlist, vhead)) {
	size_t i;
	size_t nused = vp->nused;
	const int *cmd = vp->cmd;
	point_t *pt = vp->pt;

	for (i = 0; i < nused; i++, cmd++, pt++) {
	    switch (*cmd) {
		case BG_VLIST_POLY_START:
		case BG_VLIST_TRI_START:
		    break;
		case BG_VLIST_POLY_MOVE:
		case BG_VLIST_LINE_MOVE:
		case BG_VLIST_TRI_MOVE:
		    pdv_3move(fp, *pt);
		    break;
		case BG_VLIST_POLY_DRAW:
		case BG_VLIST_POLY_END:
		case BG_VLIST_LINE_DRAW:
		case BG_VLIST_TRI_DRAW:
		case BG_VLIST_TRI_END:
		    pdv_3cont(fp, *pt);
		    break;
		default:
		    bu_log("bg_vlist_to_uplot: unknown vlist cmd x%x\n", *cmd);
	    }
	}
    }
}


void
bg_plot_vlblock(FILE *fp, const struct bg_vlblock *vbp)
{
    size_t i;

    BG_CK_VLBLOCK(vbp);

    for (i = 0; i < vbp->nused; i++) {
	if (vbp->rgb[i] == 0)
	    continue;
	if (BU_LIST_IS_EMPTY(&(vbp->head[i])))
	    continue;
	pl_color(fp,
		 (vbp->rgb[i] >> 16) & 0xFF,
		 (vbp->rgb[i] >> 8) & 0xFF,
		 (vbp->rgb[i]) & 0xFF);
	bg_vlist_to_uplot(fp, &(vbp->head[i]));
    }
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
