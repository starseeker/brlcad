/*                      V L I S T . H
 * BRL-CAD
 *
 * Copyright (c) 1993-2026 United States Government as represented by
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
/** @file rt/vlist.h
 *
 */

#ifndef RT_VLIST_H
#define RT_VLIST_H

#include "common.h"

/* system headers */
#include <stdio.h> /* for FILE */

/* interface headers */
#include "vmath.h"
#include "bg/vlist.h"
#include "bu/vls.h"
#include "rt/defines.h"

__BEGIN_DECLS

/*
 * Neutral librt-facing aliases for low-level geometric command-list storage.
 * The container and pure helper implementation are owned by libbg; librt owns
 * RT-named transitional text/uplot bridges while legacy diagnostic plot output
 * and flattened text geometry are still supported.
 */
typedef bg_vlist rt_vlist;
typedef struct bg_vlblock rt_vlblock;

#define RT_VLIST_NULL BG_VLIST_NULL
#define RT_VLIST_MAGIC BG_VLIST_MAGIC
#define RT_CK_VLIST(_p) BG_CK_VLIST(_p)
#define RT_VLBLOCK_MAGIC BG_VLBLOCK_MAGIC
#define RT_CK_VLBLOCK(_p) BG_CK_VLBLOCK(_p)

#define RT_VLIST_LINE_MOVE BG_VLIST_LINE_MOVE
#define RT_VLIST_LINE_DRAW BG_VLIST_LINE_DRAW
#define RT_VLIST_POLY_START BG_VLIST_POLY_START
#define RT_VLIST_POLY_MOVE BG_VLIST_POLY_MOVE
#define RT_VLIST_POLY_DRAW BG_VLIST_POLY_DRAW
#define RT_VLIST_POLY_END BG_VLIST_POLY_END
#define RT_VLIST_POLY_VERTNORM BG_VLIST_POLY_VERTNORM
#define RT_VLIST_TRI_START BG_VLIST_TRI_START
#define RT_VLIST_TRI_MOVE BG_VLIST_TRI_MOVE
#define RT_VLIST_TRI_DRAW BG_VLIST_TRI_DRAW
#define RT_VLIST_TRI_END BG_VLIST_TRI_END
#define RT_VLIST_TRI_VERTNORM BG_VLIST_TRI_VERTNORM
#define RT_VLIST_POINT_DRAW BG_VLIST_POINT_DRAW
#define RT_VLIST_POINT_SIZE BG_VLIST_POINT_SIZE
#define RT_VLIST_LINE_WIDTH BG_VLIST_LINE_WIDTH
#define RT_VLIST_DISPLAY_MAT BG_VLIST_DISPLAY_MAT
#define RT_VLIST_MODEL_MAT BG_VLIST_MODEL_MAT
#define RT_VLIST_CMD_MAX BG_VLIST_CMD_MAX

#define RT_GET_VLIST(_free_hd, p) BG_GET_VLIST((_free_hd), (p))
#define RT_FREE_VLIST(_free_hd, hd) BG_FREE_VLIST((_free_hd), (hd))
#define RT_ADD_VLIST(_free_hd, _dest_hd, pnt, draw) \
    BG_ADD_VLIST((_free_hd), (_dest_hd), (pnt), (draw))
#define RT_VLIST_SET_DISP_MAT(_free_hd, _dest_hd, _ref_pt) \
    BG_VLIST_SET_DISP_MAT((_free_hd), (_dest_hd), (_ref_pt))
#define RT_VLIST_SET_MODEL_MAT(_free_hd, _dest_hd) \
    BG_VLIST_SET_MODEL_MAT((_free_hd), (_dest_hd))
#define RT_VLIST_SET_POINT_SIZE(_free_hd, _dest_hd, _size) \
    BG_VLIST_SET_POINT_SIZE((_free_hd), (_dest_hd), (_size))
#define RT_VLIST_SET_LINE_WIDTH(_free_hd, _dest_hd, _width) \
    BG_VLIST_SET_LINE_WIDTH((_free_hd), (_dest_hd), (_width))

/**
 * TODO - replace these with the appropriate libbn calls specifically
 * passing &rt_vlfree
 */
RT_EXPORT extern void rt_vlist_copy(struct bu_list *dest,
				    const struct bu_list *src);
RT_EXPORT extern void rt_vlist_cleanup(void);
RT_EXPORT extern void rt_vlist_import(struct bu_list *hp,
				      struct bu_vls *namevls,
				      const unsigned char *buf);
RT_EXPORT extern rt_vlblock *    rt_vlblock_init(void);

RT_EXPORT extern void rt_vlist_3string(struct bu_list *vhead,
				       struct bu_list *free_hd,
				       const char *string,
				       const point_t origin,
				       const mat_t rot,
				       double scale);
RT_EXPORT extern void rt_vlist_2string(struct bu_list *vhead,
				       struct bu_list *free_hd,
				       const char *string,
				       double x,
				       double y,
				       double scale,
				       double theta);
RT_EXPORT extern void rt_vlist_to_uplot(FILE *fp,
					const struct bu_list *vhead);
RT_EXPORT extern void rt_plot_vlblock(FILE *fp,
				      const rt_vlblock *vbp);



/************************************************************************
 *                                                                      *
 *                      UNIX-Plot VLIST import/export routines          *
 *                                                                      *
 ************************************************************************/

RT_EXPORT extern int rt_process_uplot_value(struct bu_list **vhead,
					    rt_vlblock *vbp,
					    FILE *fp,
					    int c,
					    double char_size,
					    int mode);


/**
 * Read a BRL-style 3-D UNIX-plot file into a vector list.  For now,
 * discard color information, only extract vectors.  This might be
 * more naturally located in mged/plot.c
 */
RT_EXPORT extern int rt_uplot_to_vlist(rt_vlblock *vbp,
				       FILE *fp,
				       double char_size,
				       int mode);


/**
 * Used by MGED's "labelvert" command.
 */
RT_EXPORT extern void rt_label_vlist_verts(rt_vlblock *vbp,
					   struct bu_list *src,
					   mat_t mat,
					   double sz,
					   double mm2local);

/**
 * Used by MGED's "labelface" command.
 */
RT_EXPORT extern void rt_label_vlist_faces(rt_vlblock *vbp,
					   struct bu_list* f_list,
					   mat_t mat,
					   double sz,
					   double mm2local);


__END_DECLS

#endif /* RT_VLIST_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
