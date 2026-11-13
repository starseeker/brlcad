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

#include <limits.h>

#include "imgstream/fb_compat.h"
#include "imgstream/fbserv.h"


static int
imgstream_backend_info(void *ctx, struct fbserv_fb_info *info)
{
    imgstream_fb_t *fb = (imgstream_fb_t *)ctx;
    if (!fb || !info)
	return -1;
    size_t width = imgstream_fb_width(fb);
    size_t height = imgstream_fb_height(fb);
    if (width > INT_MAX || height > INT_MAX)
	return -1;
    info->max_width = info->width = (int)width;
    info->max_height = info->height = (int)height;
    return 0;
}

static int imgstream_backend_clear(void *ctx, const unsigned char rgb[3]) { return imgstream_fb_clear((imgstream_fb_t *)ctx, rgb); }
static ssize_t imgstream_backend_read(void *ctx, int x, int y, unsigned char *rgb, size_t count) { return imgstream_fb_read((imgstream_fb_t *)ctx, x, y, rgb, count); }
static ssize_t imgstream_backend_write(void *ctx, int x, int y, const unsigned char *rgb, size_t count) { return imgstream_fb_write((imgstream_fb_t *)ctx, x, y, rgb, count); }
static int imgstream_backend_readrect(void *ctx, int x, int y, int width, int height, unsigned char *rgb) { return imgstream_fb_readrect((imgstream_fb_t *)ctx, x, y, width, height, rgb); }
static int imgstream_backend_writerect(void *ctx, int x, int y, int width, int height, const unsigned char *rgb) { return imgstream_fb_writerect((imgstream_fb_t *)ctx, x, y, width, height, rgb); }
static int imgstream_backend_bwreadrect(void *ctx, int x, int y, int width, int height, unsigned char *bw) { return imgstream_fb_bwreadrect((imgstream_fb_t *)ctx, x, y, width, height, bw); }
static int imgstream_backend_bwwriterect(void *ctx, int x, int y, int width, int height, const unsigned char *bw) { return imgstream_fb_bwwriterect((imgstream_fb_t *)ctx, x, y, width, height, bw); }
static int imgstream_backend_cursor(void *ctx, int mode, int x, int y) { return imgstream_fb_cursor((imgstream_fb_t *)ctx, mode, x, y); }
static int imgstream_backend_getcursor(void *ctx, int *mode, int *x, int *y) { return imgstream_fb_getcursor((imgstream_fb_t *)ctx, mode, x, y); }
static int imgstream_backend_setcursor(void *ctx, const unsigned char *bits, int xbits, int ybits, int xorig, int yorig) { return imgstream_fb_setcursor((imgstream_fb_t *)ctx, bits, xbits, ybits, xorig, yorig); }
static int imgstream_backend_scursor(void *ctx, int mode, int x, int y) { return imgstream_fb_scursor((imgstream_fb_t *)ctx, mode, x, y); }
static int imgstream_backend_window(void *ctx, int xcenter, int ycenter) { return imgstream_fb_window((imgstream_fb_t *)ctx, xcenter, ycenter); }
static int imgstream_backend_zoom(void *ctx, int xzoom, int yzoom) { return imgstream_fb_zoom((imgstream_fb_t *)ctx, xzoom, yzoom); }
static int imgstream_backend_view(void *ctx, int xcenter, int ycenter, int xzoom, int yzoom) { return imgstream_fb_view((imgstream_fb_t *)ctx, xcenter, ycenter, xzoom, yzoom); }
static int imgstream_backend_getview(void *ctx, int *xcenter, int *ycenter, int *xzoom, int *yzoom) { return imgstream_fb_getview((imgstream_fb_t *)ctx, xcenter, ycenter, xzoom, yzoom); }

static int
imgstream_backend_rmap(void *ctx, struct fbserv_colormap *cmap)
{
    struct imgstream_fb_colormap map;
    int ret = imgstream_fb_rmap((imgstream_fb_t *)ctx, &map);
    if (ret != 0 || !cmap)
	return ret;
    for (int i = 0; i < 256; i++) {
	cmap->red[i] = map.red[i];
	cmap->green[i] = map.green[i];
	cmap->blue[i] = map.blue[i];
    }
    return 0;
}

static int
imgstream_backend_wmap(void *ctx, const struct fbserv_colormap *cmap)
{
    if (!cmap)
	return imgstream_fb_wmap((imgstream_fb_t *)ctx, NULL);
    struct imgstream_fb_colormap map;
    for (int i = 0; i < 256; i++) {
	map.red[i] = cmap->red[i];
	map.green[i] = cmap->green[i];
	map.blue[i] = cmap->blue[i];
    }
    return imgstream_fb_wmap((imgstream_fb_t *)ctx, &map);
}

static int imgstream_backend_flush(void *ctx) { return imgstream_fb_flush((imgstream_fb_t *)ctx); }
static int imgstream_backend_poll(void *ctx) { return imgstream_fb_poll((imgstream_fb_t *)ctx); }
static int imgstream_backend_help(void *UNUSED(ctx)) { return 0; }

static const struct fbserv_fb_ops imgstream_backend_ops = {
    imgstream_backend_info,
    imgstream_backend_clear,
    imgstream_backend_read,
    imgstream_backend_write,
    imgstream_backend_readrect,
    imgstream_backend_writerect,
    imgstream_backend_bwreadrect,
    imgstream_backend_bwwriterect,
    imgstream_backend_cursor,
    imgstream_backend_getcursor,
    imgstream_backend_setcursor,
    imgstream_backend_scursor,
    imgstream_backend_window,
    imgstream_backend_zoom,
    imgstream_backend_view,
    imgstream_backend_getview,
    imgstream_backend_rmap,
    imgstream_backend_wmap,
    imgstream_backend_flush,
    imgstream_backend_poll,
    imgstream_backend_help
};


int
imgstream_fbserv_set_framebuffer(struct fbserv_obj *fbsp, imgstream_fb_t *fb)
{
    if (!fbsp)
	return BRLCAD_ERROR;
    if (!fb) {
	if (fbsp->fbs_fb_ops == &imgstream_backend_ops)
	    fbserv_clear_backend(fbsp);
	return BRLCAD_OK;
    }
    return fbserv_set_backend(fbsp, &imgstream_backend_ops, fb);
}


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
