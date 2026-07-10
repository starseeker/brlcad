/*                 F B S E R V _ L E G A C Y . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "dm.h"
#include "dm/fbserv_legacy.h"

static int
legacy_info(void *ctx, struct fbserv_fb_info *info)
{
    struct fb *fbp = (struct fb *)ctx;
    if (!fbp || !info)
	return -1;
    info->max_width = fb_get_max_width(fbp);
    info->max_height = fb_get_max_height(fbp);
    info->width = fb_getwidth(fbp);
    info->height = fb_getheight(fbp);
    return 0;
}

static int
legacy_clear(void *ctx, const unsigned char rgb[3])
{
    return fb_clear((struct fb *)ctx, (unsigned char *)rgb);
}

static ssize_t
legacy_read(void *ctx, int x, int y, unsigned char *rgb, size_t count)
{
    return fb_read((struct fb *)ctx, x, y, rgb, count);
}

static ssize_t
legacy_write(void *ctx, int x, int y, const unsigned char *rgb, size_t count)
{
    return fb_write((struct fb *)ctx, x, y, rgb, count);
}

static int
legacy_readrect(void *ctx, int x, int y, int width, int height,
	unsigned char *rgb)
{
    return fb_readrect((struct fb *)ctx, x, y, width, height, rgb);
}

static int
legacy_writerect(void *ctx, int x, int y, int width, int height,
	const unsigned char *rgb)
{
    return fb_writerect((struct fb *)ctx, x, y, width, height, rgb);
}

static int
legacy_bwreadrect(void *ctx, int x, int y, int width, int height,
	unsigned char *bw)
{
    return fb_bwreadrect((struct fb *)ctx, x, y, width, height, bw);
}

static int
legacy_bwwriterect(void *ctx, int x, int y, int width, int height,
	const unsigned char *bw)
{
    return fb_bwwriterect((struct fb *)ctx, x, y, width, height, bw);
}

static int
legacy_cursor(void *ctx, int mode, int x, int y)
{
    return fb_cursor((struct fb *)ctx, mode, x, y);
}

static int
legacy_getcursor(void *ctx, int *mode, int *x, int *y)
{
    return fb_getcursor((struct fb *)ctx, mode, x, y);
}

static int
legacy_setcursor(void *ctx, const unsigned char *bits, int xbits, int ybits,
	int xorig, int yorig)
{
    return fb_setcursor((struct fb *)ctx, bits, xbits, ybits, xorig, yorig);
}

static int
legacy_scursor(void *ctx, int mode, int x, int y)
{
    return fb_scursor((struct fb *)ctx, mode, x, y);
}

static int
legacy_window(void *ctx, int xcenter, int ycenter)
{
    return fb_window((struct fb *)ctx, xcenter, ycenter);
}

static int
legacy_zoom(void *ctx, int xzoom, int yzoom)
{
    return fb_zoom((struct fb *)ctx, xzoom, yzoom);
}

static int
legacy_view(void *ctx, int xcenter, int ycenter, int xzoom, int yzoom)
{
    return fb_view((struct fb *)ctx, xcenter, ycenter, xzoom, yzoom);
}

static int
legacy_getview(void *ctx, int *xcenter, int *ycenter, int *xzoom, int *yzoom)
{
    return fb_getview((struct fb *)ctx, xcenter, ycenter, xzoom, yzoom);
}

static int
legacy_rmap(void *ctx, struct fbserv_colormap *cmap)
{
    ColorMap map;
    int ret = fb_rmap((struct fb *)ctx, &map);
    if (ret < 0 || !cmap)
	return ret;
    for (int i = 0; i < 256; i++) {
	cmap->red[i] = map.cm_red[i];
	cmap->green[i] = map.cm_green[i];
	cmap->blue[i] = map.cm_blue[i];
    }
    return ret;
}

static int
legacy_wmap(void *ctx, const struct fbserv_colormap *cmap)
{
    if (!cmap)
	return fb_wmap((struct fb *)ctx, COLORMAP_NULL);

    ColorMap map;
    for (int i = 0; i < 256; i++) {
	map.cm_red[i] = cmap->red[i];
	map.cm_green[i] = cmap->green[i];
	map.cm_blue[i] = cmap->blue[i];
    }
    return fb_wmap((struct fb *)ctx, &map);
}

static int
legacy_flush(void *ctx)
{
    return fb_flush((struct fb *)ctx);
}

static int
legacy_poll(void *ctx)
{
    return fb_poll((struct fb *)ctx);
}

static int
legacy_help(void *ctx)
{
    return fb_help((struct fb *)ctx);
}

static const struct fbserv_fb_ops legacy_ops = {
    .info = legacy_info,
    .clear = legacy_clear,
    .read = legacy_read,
    .write = legacy_write,
    .readrect = legacy_readrect,
    .writerect = legacy_writerect,
    .bwreadrect = legacy_bwreadrect,
    .bwwriterect = legacy_bwwriterect,
    .cursor = legacy_cursor,
    .getcursor = legacy_getcursor,
    .setcursor = legacy_setcursor,
    .scursor = legacy_scursor,
    .window = legacy_window,
    .zoom = legacy_zoom,
    .view = legacy_view,
    .getview = legacy_getview,
    .rmap = legacy_rmap,
    .wmap = legacy_wmap,
    .flush = legacy_flush,
    .poll = legacy_poll,
    .help = legacy_help
};

int
dm_fbserv_set_framebuffer(struct fbserv_obj *fbsp, struct fb *fbp)
{
    if (!fbsp)
	return -1;

    if (!fbp) {
	if (fbsp->fbs_fb_ops == &legacy_ops)
	    fbserv_clear_backend(fbsp);
	fbserv_set_legacy_framebuffer(fbsp, NULL);
	return 0;
    }

    fbserv_set_legacy_framebuffer(fbsp, fbp);
    return fbserv_set_backend(fbsp, &legacy_ops, fbp);
}
