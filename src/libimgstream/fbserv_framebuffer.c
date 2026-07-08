/*              F B S E R V _ F R A M E B U F F E R . C
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
/** @file libimgstream/fbserv_framebuffer.c */

#include "common.h"

#include "imgstream/fbserv.h"


static int
fbserv_has_backend_op(const struct fbserv_obj *fbsp)
{
    return fbsp && fbsp->fbs_fb_ops && fbsp->fbs_fb_ctx;
}


int
fbserv_backend_op_installed(const struct fbserv_obj *fbsp, enum fbserv_backend_op op)
{
    const struct fbserv_fb_ops *ops;

    if (!fbserv_has_backend_op(fbsp))
	return 0;

    ops = fbsp->fbs_fb_ops;
    switch (op) {
	case FBSERV_BACKEND_OP_INFO:
	    return ops->info != NULL;
	case FBSERV_BACKEND_OP_CLEAR:
	    return ops->clear != NULL;
	case FBSERV_BACKEND_OP_READ:
	    return ops->read != NULL;
	case FBSERV_BACKEND_OP_WRITE:
	    return ops->write != NULL;
	case FBSERV_BACKEND_OP_READRECT:
	    return ops->readrect != NULL;
	case FBSERV_BACKEND_OP_WRITERECT:
	    return ops->writerect != NULL;
	case FBSERV_BACKEND_OP_BWREADRECT:
	    return ops->bwreadrect != NULL;
	case FBSERV_BACKEND_OP_BWWRITERECT:
	    return ops->bwwriterect != NULL;
	case FBSERV_BACKEND_OP_CURSOR:
	    return ops->cursor != NULL;
	case FBSERV_BACKEND_OP_GETCURSOR:
	    return ops->getcursor != NULL;
	case FBSERV_BACKEND_OP_SETCURSOR:
	    return ops->setcursor != NULL;
	case FBSERV_BACKEND_OP_SCURSOR:
	    return ops->scursor != NULL;
	case FBSERV_BACKEND_OP_WINDOW:
	    return ops->window != NULL;
	case FBSERV_BACKEND_OP_ZOOM:
	    return ops->zoom != NULL;
	case FBSERV_BACKEND_OP_VIEW:
	    return ops->view != NULL;
	case FBSERV_BACKEND_OP_GETVIEW:
	    return ops->getview != NULL;
	case FBSERV_BACKEND_OP_RMAP:
	    return ops->rmap != NULL;
	case FBSERV_BACKEND_OP_WMAP:
	    return ops->wmap != NULL;
	case FBSERV_BACKEND_OP_FLUSH:
	    return ops->flush != NULL;
	case FBSERV_BACKEND_OP_POLL:
	    return ops->poll != NULL;
	case FBSERV_BACKEND_OP_HELP:
	    return ops->help != NULL;
    }

    return 0;
}


int
fbserv_backend_info(struct fbserv_obj *fbsp, struct fbserv_fb_info *info)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->info)
	return fbsp->fbs_fb_ops->info(fbsp->fbs_fb_ctx, info);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_clear(struct fbserv_obj *fbsp, const unsigned char rgb[3])
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->clear)
	return fbsp->fbs_fb_ops->clear(fbsp->fbs_fb_ctx, rgb);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


ssize_t
fbserv_backend_read(struct fbserv_obj *fbsp,
	int x,
	int y,
	unsigned char *rgb,
	size_t count)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->read)
	return fbsp->fbs_fb_ops->read(fbsp->fbs_fb_ctx, x, y, rgb, count);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


ssize_t
fbserv_backend_write(struct fbserv_obj *fbsp,
	int x,
	int y,
	const unsigned char *rgb,
	size_t count)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->write)
	return fbsp->fbs_fb_ops->write(fbsp->fbs_fb_ctx, x, y, rgb, count);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_readrect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	unsigned char *rgb)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->readrect)
	return fbsp->fbs_fb_ops->readrect(fbsp->fbs_fb_ctx, xmin, ymin,
	    width, height, rgb);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_writerect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	const unsigned char *rgb)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->writerect)
	return fbsp->fbs_fb_ops->writerect(fbsp->fbs_fb_ctx, xmin, ymin,
	    width, height, rgb);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_bwreadrect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	unsigned char *bw)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->bwreadrect)
	return fbsp->fbs_fb_ops->bwreadrect(fbsp->fbs_fb_ctx, xmin, ymin,
	    width, height, bw);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_bwwriterect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	const unsigned char *bw)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->bwwriterect)
	return fbsp->fbs_fb_ops->bwwriterect(fbsp->fbs_fb_ctx, xmin, ymin,
	    width, height, bw);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_cursor(struct fbserv_obj *fbsp, int mode, int x, int y)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->cursor)
	return fbsp->fbs_fb_ops->cursor(fbsp->fbs_fb_ctx, mode, x, y);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_getcursor(struct fbserv_obj *fbsp, int *mode, int *x, int *y)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->getcursor)
	return fbsp->fbs_fb_ops->getcursor(fbsp->fbs_fb_ctx, mode, x, y);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_setcursor(struct fbserv_obj *fbsp,
	const unsigned char *bits,
	int xbits,
	int ybits,
	int xorig,
	int yorig)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->setcursor)
	return fbsp->fbs_fb_ops->setcursor(fbsp->fbs_fb_ctx, bits, xbits,
	    ybits, xorig, yorig);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_scursor(struct fbserv_obj *fbsp, int mode, int x, int y)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->scursor)
	return fbsp->fbs_fb_ops->scursor(fbsp->fbs_fb_ctx, mode, x, y);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_window(struct fbserv_obj *fbsp, int xcenter, int ycenter)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->window)
	return fbsp->fbs_fb_ops->window(fbsp->fbs_fb_ctx, xcenter, ycenter);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_zoom(struct fbserv_obj *fbsp, int xzoom, int yzoom)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->zoom)
	return fbsp->fbs_fb_ops->zoom(fbsp->fbs_fb_ctx, xzoom, yzoom);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_view(struct fbserv_obj *fbsp,
	int xcenter,
	int ycenter,
	int xzoom,
	int yzoom)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->view)
	return fbsp->fbs_fb_ops->view(fbsp->fbs_fb_ctx, xcenter, ycenter,
	    xzoom, yzoom);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_getview(struct fbserv_obj *fbsp,
	int *xcenter,
	int *ycenter,
	int *xzoom,
	int *yzoom)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->getview)
	return fbsp->fbs_fb_ops->getview(fbsp->fbs_fb_ctx, xcenter, ycenter,
	    xzoom, yzoom);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_rmap(struct fbserv_obj *fbsp, struct fbserv_colormap *cmap)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->rmap)
	return fbsp->fbs_fb_ops->rmap(fbsp->fbs_fb_ctx, cmap);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_wmap(struct fbserv_obj *fbsp,
	const struct fbserv_colormap *cmap)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->wmap)
	return fbsp->fbs_fb_ops->wmap(fbsp->fbs_fb_ctx, cmap);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_flush(struct fbserv_obj *fbsp)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->flush)
	return fbsp->fbs_fb_ops->flush(fbsp->fbs_fb_ctx);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_poll(struct fbserv_obj *fbsp)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->poll)
	return fbsp->fbs_fb_ops->poll(fbsp->fbs_fb_ctx);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
}


int
fbserv_backend_help(struct fbserv_obj *fbsp)
{
    if (fbserv_has_backend_op(fbsp) && fbsp->fbs_fb_ops->help)
	return fbsp->fbs_fb_ops->help(fbsp->fbs_fb_ctx);
    return FBSERV_FRAMEBUFFER_UNHANDLED;
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
