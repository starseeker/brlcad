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
/** @addtogroup libbsg
 *
 * @brief
 * Legacy BSG vlist compatibility names.
 *
 * The low-level vlist container is now defined by bg/vlist.h.  This header
 * preserves the historical BSG names and declares the existing legacy helper
 * functions until their users migrate to the BG/RT-facing API.
 */
/** @{ */
/* @file bsg/vlist.h */

#ifndef BSG_VLIST_H
#define BSG_VLIST_H

#include "common.h"
#include <stdio.h>
#include "bg/vlist.h"
#include "vmath.h"
#include "bu/magic.h"
#include "bu/list.h"
#include "bu/malloc.h"
#include "bu/vls.h"
#include "bsg/defines.h"

__BEGIN_DECLS

#define bsg_vlist bg_vlist
#define bsg_vlblock bg_vlblock

#define BSG_VLIST_CHUNK BG_VLIST_CHUNK
#define BSG_VLIST_NULL BG_VLIST_NULL
#define BSG_VLIST_MAGIC BG_VLIST_MAGIC
#define BSG_CK_VLIST(_p) BG_CK_VLIST(_p)

#define BSG_VLIST_LINE_MOVE BG_VLIST_LINE_MOVE
#define BSG_VLIST_LINE_DRAW BG_VLIST_LINE_DRAW
#define BSG_VLIST_POLY_START BG_VLIST_POLY_START
#define BSG_VLIST_POLY_MOVE BG_VLIST_POLY_MOVE
#define BSG_VLIST_POLY_DRAW BG_VLIST_POLY_DRAW
#define BSG_VLIST_POLY_END BG_VLIST_POLY_END
#define BSG_VLIST_POLY_VERTNORM BG_VLIST_POLY_VERTNORM
#define BSG_VLIST_TRI_START BG_VLIST_TRI_START
#define BSG_VLIST_TRI_MOVE BG_VLIST_TRI_MOVE
#define BSG_VLIST_TRI_DRAW BG_VLIST_TRI_DRAW
#define BSG_VLIST_TRI_END BG_VLIST_TRI_END
#define BSG_VLIST_TRI_VERTNORM BG_VLIST_TRI_VERTNORM
#define BSG_VLIST_POINT_DRAW BG_VLIST_POINT_DRAW
#define BSG_VLIST_POINT_SIZE BG_VLIST_POINT_SIZE
#define BSG_VLIST_LINE_WIDTH BG_VLIST_LINE_WIDTH
#define BSG_VLIST_DISPLAY_MAT BG_VLIST_DISPLAY_MAT
#define BSG_VLIST_MODEL_MAT BG_VLIST_MODEL_MAT
#define BSG_VLIST_CMD_MAX BG_VLIST_CMD_MAX

#define BSG_GET_VLIST(_free_hd, p) BG_GET_VLIST((_free_hd), (p))
#define BSG_FREE_VLIST(_free_hd, hd) BG_FREE_VLIST((_free_hd), (hd))
#define BSG_ADD_VLIST(_free_hd, _dest_hd, pnt, draw) \
    BG_ADD_VLIST((_free_hd), (_dest_hd), (pnt), (draw))
#define BSG_VLIST_SET_DISP_MAT(_free_hd, _dest_hd, _ref_pt) \
    BG_VLIST_SET_DISP_MAT((_free_hd), (_dest_hd), (_ref_pt))
#define BSG_VLIST_SET_MODEL_MAT(_free_hd, _dest_hd) \
    BG_VLIST_SET_MODEL_MAT((_free_hd), (_dest_hd))
#define BSG_VLIST_SET_POINT_SIZE(_free_hd, _dest_hd, _size) \
    BG_VLIST_SET_POINT_SIZE((_free_hd), (_dest_hd), (_size))
#define BSG_VLIST_SET_LINE_WIDTH(_free_hd, _dest_hd, _width) \
    BG_VLIST_SET_LINE_WIDTH((_free_hd), (_dest_hd), (_width))

#define BSG_VLBLOCK_MAGIC BG_VLBLOCK_MAGIC
#define BSG_CK_VLBLOCK(_p) BG_CK_VLBLOCK(_p)

/* bsg_vlist_* and bsg_vlblock_* function declarations until the rename slice. */
BSG_EXPORT extern size_t bsg_vlist_cmd_cnt(bsg_vlist *vlist);
BSG_EXPORT extern int bsg_vlist_bbox(struct bu_list *vlistp, point_t *bmin, point_t *bmax, size_t *length, int *dispmode);
BSG_EXPORT extern void bsg_vlist_3string(struct bu_list *vhead, struct bu_list *free_hd, const char *string, const point_t origin, const mat_t rot, double scale);
BSG_EXPORT extern void bsg_vlist_2string(struct bu_list *vhead, struct bu_list *free_hd, const char *string, double x, double y, double scale, double theta);
BSG_EXPORT extern void bsg_vlist_arb8(struct bu_list *vhead, struct bu_list *vlfree, point_t pts[8]);
BSG_EXPORT extern const char *bsg_vlist_get_cmd_description(int cmd);
BSG_EXPORT extern size_t bsg_ck_vlist(const struct bu_list *vhead);
BSG_EXPORT extern void bsg_vlist_copy(struct bu_list *vlists, struct bu_list *dest, const struct bu_list *src);
BSG_EXPORT extern void bsg_vlist_export(struct bu_vls *vls, struct bu_list *hp, const char *name);
BSG_EXPORT extern void bsg_vlist_import(struct bu_list *vlists, struct bu_list *hp, struct bu_vls *namevls, const unsigned char *buf);
BSG_EXPORT extern void bsg_vlist_cleanup(struct bu_list *hd);
BSG_EXPORT extern struct bsg_vlblock *bsg_vlblock_init(struct bu_list *free_vlist_hd, int max_ent);
BSG_EXPORT extern void bsg_vlblock_free(struct bsg_vlblock *vbp);
BSG_EXPORT extern struct bu_list *bsg_vlblock_find(struct bsg_vlblock *vbp, int r, int g, int b);
BSG_EXPORT void bsg_vlist_rpp(struct bu_list *vlists, struct bu_list *hd, const point_t minn, const point_t maxx);
BSG_EXPORT extern void bsg_plot_vlblock(FILE *fp, const struct bsg_vlblock *vbp);
BSG_EXPORT extern void bsg_vlist_to_uplot(FILE *fp, const struct bu_list *vhead);

__END_DECLS

#endif /* BSG_VLIST_H */

/** @} */
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
