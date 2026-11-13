/*                        V L I S T . H
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
/** @file bg/vlist.h
 *
 * Low-level geometric vertex command list storage.
 */

#ifndef BG_VLIST_H
#define BG_VLIST_H

#include "common.h"

#include <stdio.h> /* for FILE */

#include "vmath.h"
#include "bu/list.h"
#include "bu/magic.h"
#include "bu/malloc.h"
#include "bu/vls.h"
#include "bg/defines.h"

__BEGIN_DECLS

#define BG_VLIST_CHUNK 35   /**< @brief 32-bit mach => just less than 1k */

struct bg_vlist {
    struct bu_list l;               /**< @brief magic, forw, back */
    size_t nused;                   /**< @brief elements 0..nused active */
    int cmd[BG_VLIST_CHUNK];        /**< @brief BG_VLIST_* command codes */
    point_t pt[BG_VLIST_CHUNK];     /**< @brief associated 3-point/vect */
};
typedef struct bg_vlist bg_vlist;

#define BG_VLIST_NULL  ((bg_vlist *)0)
#define BG_VLIST_MAGIC 0x98237474
#define BG_CK_VLIST(_p) BU_CKMAG((_p), BG_VLIST_MAGIC, "bg_vlist")

/* Geometric command constants. */
#define BG_VLIST_LINE_MOVE         0   /**< @brief specify new line */
#define BG_VLIST_LINE_DRAW         1   /**< @brief subsequent line vertex */
#define BG_VLIST_POLY_START        2   /**< @brief pt[] has surface normal */
#define BG_VLIST_POLY_MOVE         3   /**< @brief move to first poly vertex */
#define BG_VLIST_POLY_DRAW         4   /**< @brief subsequent poly vertex */
#define BG_VLIST_POLY_END          5   /**< @brief last vert (repeats 1st), draw poly */
#define BG_VLIST_POLY_VERTNORM     6   /**< @brief per-vertex normal */
#define BG_VLIST_TRI_START         7   /**< @brief pt[] has surface normal */
#define BG_VLIST_TRI_MOVE          8   /**< @brief move to first triangle vertex */
#define BG_VLIST_TRI_DRAW          9   /**< @brief subsequent triangle vertex */
#define BG_VLIST_TRI_END           10  /**< @brief last vert (repeats 1st), draw poly */
#define BG_VLIST_TRI_VERTNORM      11  /**< @brief per-vertex normal */
#define BG_VLIST_POINT_DRAW        12  /**< @brief draw a single point */
#define BG_VLIST_POINT_SIZE        13  /**< @brief specify point pixel size */
#define BG_VLIST_LINE_WIDTH        14  /**< @brief specify line pixel width */
#define BG_VLIST_DISPLAY_MAT       15  /**< @brief specify the model matrix */
#define BG_VLIST_MODEL_MAT         16  /**< @brief specify the display matrix */
#define BG_VLIST_CMD_MAX           16  /**< @brief max command number */

/* Applications must call BU_LIST_INIT on their free-list head before use.
 * These macros are non-PARALLEL. */
#define BG_GET_VLIST(_free_hd, p) do { \
(p) = BU_LIST_FIRST(bg_vlist, (_free_hd)); \
if (BU_LIST_IS_HEAD((p), (_free_hd))) { \
    BU_ALLOC((p), bg_vlist); \
    (p)->l.magic = BG_VLIST_MAGIC; \
} else { \
    BU_LIST_DEQUEUE(&((p)->l)); \
} \
(p)->nused = 0; \
    } while (0)

#define BG_FREE_VLIST(_free_hd, hd) do { \
BU_CK_LIST_HEAD((hd)); \
BU_LIST_APPEND_LIST((_free_hd), (hd)); \
    } while (0)

#define BG_ADD_VLIST(_free_hd, _dest_hd, pnt, draw) do { \
bg_vlist *_vp; \
BU_CK_LIST_HEAD(_dest_hd); \
_vp = BU_LIST_LAST(bg_vlist, (_dest_hd)); \
if (BU_LIST_IS_HEAD(_vp, (_dest_hd)) || _vp->nused >= BG_VLIST_CHUNK) { \
    BG_GET_VLIST(_free_hd, _vp); \
    BU_LIST_INSERT((_dest_hd), &(_vp->l)); \
} \
VMOVE(_vp->pt[_vp->nused], (pnt)); \
_vp->cmd[_vp->nused++] = (draw); \
    } while (0)

#define BG_VLIST_SET_DISP_MAT(_free_hd, _dest_hd, _ref_pt) do { \
bg_vlist *_vp; \
BU_CK_LIST_HEAD(_dest_hd); \
_vp = BU_LIST_LAST(bg_vlist, (_dest_hd)); \
if (BU_LIST_IS_HEAD(_vp, (_dest_hd)) || _vp->nused >= BG_VLIST_CHUNK) { \
    BG_GET_VLIST(_free_hd, _vp); \
    BU_LIST_INSERT((_dest_hd), &(_vp->l)); \
} \
VMOVE(_vp->pt[_vp->nused], (_ref_pt)); \
_vp->cmd[_vp->nused++] = BG_VLIST_DISPLAY_MAT; \
    } while (0)

#define BG_VLIST_SET_MODEL_MAT(_free_hd, _dest_hd) do { \
bg_vlist *_vp; \
BU_CK_LIST_HEAD(_dest_hd); \
_vp = BU_LIST_LAST(bg_vlist, (_dest_hd)); \
if (BU_LIST_IS_HEAD(_vp, (_dest_hd)) || _vp->nused >= BG_VLIST_CHUNK) { \
    BG_GET_VLIST(_free_hd, _vp); \
    BU_LIST_INSERT((_dest_hd), &(_vp->l)); \
} \
_vp->cmd[_vp->nused++] = BG_VLIST_MODEL_MAT; \
    } while (0)

#define BG_VLIST_SET_POINT_SIZE(_free_hd, _dest_hd, _size) do { \
bg_vlist *_vp; \
BU_CK_LIST_HEAD(_dest_hd); \
_vp = BU_LIST_LAST(bg_vlist, (_dest_hd)); \
if (BU_LIST_IS_HEAD(_vp, (_dest_hd)) || _vp->nused >= BG_VLIST_CHUNK) { \
    BG_GET_VLIST(_free_hd, _vp); \
    BU_LIST_INSERT((_dest_hd), &(_vp->l)); \
} \
_vp->pt[_vp->nused][0] = (_size); \
_vp->cmd[_vp->nused++] = BG_VLIST_POINT_SIZE; \
    } while (0)

#define BG_VLIST_SET_LINE_WIDTH(_free_hd, _dest_hd, _width) do { \
bg_vlist *_vp; \
BU_CK_LIST_HEAD(_dest_hd); \
_vp = BU_LIST_LAST(bg_vlist, (_dest_hd)); \
if (BU_LIST_IS_HEAD(_vp, (_dest_hd)) || _vp->nused >= BG_VLIST_CHUNK) { \
    BG_GET_VLIST(_free_hd, _vp); \
    BU_LIST_INSERT((_dest_hd), &(_vp->l)); \
} \
_vp->pt[_vp->nused][0] = (_width); \
_vp->cmd[_vp->nused++] = BG_VLIST_LINE_WIDTH; \
    } while (0)

struct bg_vlblock {
    uint32_t magic;
    size_t nused;
    size_t max;
    long *rgb;                      /**< @brief rgb[max] variable size array */
    struct bu_list *head;           /**< @brief head[max] variable size array */
    struct bu_list *free_vlist_hd;  /**< @brief where to get/put free vlists */
};
#define BG_VLBLOCK_MAGIC 0x981bd112
#define BG_CK_VLBLOCK(_p) BU_CKMAG((_p), BG_VLBLOCK_MAGIC, "bg_vlblock")

BG_EXPORT extern size_t bg_vlist_cmd_cnt(bg_vlist *vlist);
BG_EXPORT extern int bg_vlist_bbox(struct bu_list *vlistp, point_t *bmin, point_t *bmax, size_t *length, int *dispmode);
BG_EXPORT extern const char *bg_vlist_get_cmd_description(int cmd);
BG_EXPORT extern size_t bg_ck_vlist(const struct bu_list *vhead);
BG_EXPORT extern void bg_vlist_copy(struct bu_list *vlists, struct bu_list *dest, const struct bu_list *src);
BG_EXPORT extern void bg_vlist_export(struct bu_vls *vls, struct bu_list *hp, const char *name);
BG_EXPORT extern void bg_vlist_import(struct bu_list *vlists, struct bu_list *hp, struct bu_vls *namevls, const unsigned char *buf);
BG_EXPORT extern void bg_vlist_cleanup(struct bu_list *hd);
BG_EXPORT extern struct bg_vlblock *bg_vlblock_init(struct bu_list *free_vlist_hd, int max_ent);
BG_EXPORT extern void bg_vlblock_free(struct bg_vlblock *vbp);
BG_EXPORT extern struct bu_list *bg_vlblock_find(struct bg_vlblock *vbp, int r, int g, int b);
BG_EXPORT extern void bg_vlist_rpp(struct bu_list *vlists, struct bu_list *hd, const point_t minn, const point_t maxx);
BG_EXPORT extern void bg_vlist_arb8(struct bu_list *vhead, struct bu_list *vlfree, point_t pts[8]);
BG_EXPORT extern void bg_vlist_3string(struct bu_list *vhead, struct bu_list *free_hd, const char *string, const point_t origin, const mat_t rot, double scale);
BG_EXPORT extern void bg_vlist_2string(struct bu_list *vhead, struct bu_list *free_hd, const char *string, double x, double y, double scale, double theta);
BG_EXPORT extern void bg_vlist_to_uplot(FILE *fp, const struct bu_list *vhead);
BG_EXPORT extern void bg_plot_vlblock(FILE *fp, const struct bg_vlblock *vbp);

__END_DECLS

#endif /* BG_VLIST_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
