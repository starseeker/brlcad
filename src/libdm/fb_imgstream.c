/*                  F B _ I M G S T R E A M . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libdm/fb_imgstream.c
 *
 * Internal compatibility helpers for legacy struct fb implementations that can
 * delegate standalone display-device traffic to libimgstream display hosts.
 */

#include "common.h"

#include "bu/str.h"
#include "imgstream/fb_compat.h"

#include "./include/private.h"


int
fb_imgstream_compat_supported(const char *spec)
{
    return imgstream_fb_spec_supported(spec);
}


int
fb_imgstream_compat_active(const struct fb_imgstream_compat *compat)
{
    return compat && compat->fb;
}


static int
fb_imgstream_width(struct fb *ifp, int width)
{
    if (width > 0)
	return width;
    if (ifp && ifp->i && ifp->i->if_width > 0)
	return ifp->i->if_width;
    return 1;
}


static int
fb_imgstream_height(struct fb *ifp, int height)
{
    if (height > 0)
	return height;
    if (ifp && ifp->i && ifp->i->if_height > 0)
	return ifp->i->if_height;
    return 1;
}


static void
fb_imgstream_apply_state(struct fb *ifp, struct fb_imgstream_compat *compat)
{
    if (!ifp || !ifp->i || !compat || !compat->fb)
	return;

    ifp->i->if_width = ifp->i->if_max_width =
	(int)imgstream_fb_width(compat->fb);
    ifp->i->if_height = ifp->i->if_max_height =
	(int)imgstream_fb_height(compat->fb);
    ifp->i->if_xzoom = 1;
    ifp->i->if_yzoom = 1;
    ifp->i->if_xcenter = ifp->i->if_width / 2;
    ifp->i->if_ycenter = ifp->i->if_height / 2;
    ifp->i->if_cursmode = 0;
    ifp->i->if_xcurs = 0;
    ifp->i->if_ycurs = 0;
    ifp->i->if_poll_refresh_rate = imgstream_fb_poll_rate(compat->fb);
}


int
fb_imgstream_compat_open(struct fb *ifp,
	struct fb_imgstream_compat *compat,
	const char *spec,
	int width,
	int height)
{
    if (!ifp || !ifp->i || !compat || !spec || !spec[0])
	return -1;

    compat->fb = NULL;
    compat->spec = spec;
    if (!imgstream_fb_spec_supported(spec))
	return 1;

    int swidth = fb_imgstream_width(ifp, width);
    int sheight = fb_imgstream_height(ifp, height);
    imgstream_fb_t *fb = imgstream_fb_open(spec, (size_t)swidth,
	    (size_t)sheight);
    if (!fb)
	return 1;

    compat->fb = fb;
    fb_imgstream_apply_state(ifp, compat);
    return 0;
}


int
fb_imgstream_compat_close(struct fb_imgstream_compat *compat)
{
    if (!compat || !compat->fb)
	return 0;

    imgstream_fb_close(compat->fb);
    compat->fb = NULL;
    return 0;
}


int
fb_imgstream_compat_configure(struct fb *ifp,
	struct fb_imgstream_compat *compat,
	int width,
	int height)
{
    if (!ifp || !ifp->i || !compat || !compat->fb)
	return -1;

    int swidth = fb_imgstream_width(ifp, width);
    int sheight = fb_imgstream_height(ifp, height);
    if (ifp->i->if_active_clients > 0) {
	return imgstream_fb_viewport(compat->fb, 0, 0, swidth, sheight);
    }

    if (swidth == (int)imgstream_fb_width(compat->fb) &&
	    sheight == (int)imgstream_fb_height(compat->fb)) {
	return imgstream_fb_viewport(compat->fb, 0, 0, swidth, sheight);
    }

    imgstream_fb_t *new_fb = imgstream_fb_open(compat->spec, (size_t)swidth,
	    (size_t)sheight);
    if (!new_fb)
	return -1;

    imgstream_fb_close(compat->fb);
    compat->fb = new_fb;
    fb_imgstream_apply_state(ifp, compat);
    return 0;
}


int
fb_imgstream_compat_flush(struct fb_imgstream_compat *compat)
{
    return compat && compat->fb ? imgstream_fb_flush(compat->fb) : -1;
}


int
fb_imgstream_compat_poll(struct fb *ifp, struct fb_imgstream_compat *compat)
{
    if (!compat || !compat->fb)
	return -1;

    if (ifp && ifp->i)
	ifp->i->if_poll_refresh_rate = imgstream_fb_poll_rate(compat->fb);
    return imgstream_fb_poll(compat->fb);
}


int
fb_imgstream_compat_clear(struct fb_imgstream_compat *compat, unsigned char *rgb)
{
    return compat && compat->fb ? imgstream_fb_clear(compat->fb, rgb) : -1;
}


ssize_t
fb_imgstream_compat_read(const struct fb_imgstream_compat *compat,
	int x,
	int y,
	unsigned char *rgb,
	size_t count)
{
    return compat && compat->fb ?
	imgstream_fb_read(compat->fb, x, y, rgb, count) : -1;
}


ssize_t
fb_imgstream_compat_write(struct fb_imgstream_compat *compat,
	int x,
	int y,
	const unsigned char *rgb,
	size_t count)
{
    return compat && compat->fb ?
	imgstream_fb_write(compat->fb, x, y, rgb, count) : -1;
}


int
fb_imgstream_compat_readrect(const struct fb_imgstream_compat *compat,
	int xmin,
	int ymin,
	int width,
	int height,
	unsigned char *rgb)
{
    return compat && compat->fb ?
	imgstream_fb_readrect(compat->fb, xmin, ymin, width, height, rgb) : -1;
}


int
fb_imgstream_compat_writerect(struct fb_imgstream_compat *compat,
	int xmin,
	int ymin,
	int width,
	int height,
	const unsigned char *rgb)
{
    return compat && compat->fb ?
	imgstream_fb_writerect(compat->fb, xmin, ymin, width, height, rgb) : -1;
}


int
fb_imgstream_compat_bwreadrect(const struct fb_imgstream_compat *compat,
	int xmin,
	int ymin,
	int width,
	int height,
	unsigned char *bw)
{
    return compat && compat->fb ?
	imgstream_fb_bwreadrect(compat->fb, xmin, ymin, width, height, bw) : -1;
}


int
fb_imgstream_compat_bwwriterect(struct fb_imgstream_compat *compat,
	int xmin,
	int ymin,
	int width,
	int height,
	const unsigned char *bw)
{
    return compat && compat->fb ?
	imgstream_fb_bwwriterect(compat->fb, xmin, ymin, width, height, bw) : -1;
}


int
fb_imgstream_compat_view(struct fb *ifp,
	struct fb_imgstream_compat *compat,
	int xcenter,
	int ycenter,
	int xzoom,
	int yzoom)
{
    if (!ifp || !ifp->i || !compat || !compat->fb)
	return -1;

    if (xzoom < 1)
	xzoom = 1;
    if (yzoom < 1)
	yzoom = 1;
    if (xcenter < 0 || xcenter >= ifp->i->if_width)
	return -1;
    if (ycenter < 0 || ycenter >= ifp->i->if_height)
	return -1;
    if (xzoom >= ifp->i->if_width || yzoom >= ifp->i->if_height)
	return -1;

    int ret = imgstream_fb_view(compat->fb, xcenter, ycenter, xzoom, yzoom);
    if (ret == 0) {
	ifp->i->if_xcenter = xcenter;
	ifp->i->if_ycenter = ycenter;
	ifp->i->if_xzoom = xzoom;
	ifp->i->if_yzoom = yzoom;
    }
    return ret;
}


int
fb_imgstream_compat_getview(const struct fb_imgstream_compat *compat,
	int *xcenter,
	int *ycenter,
	int *xzoom,
	int *yzoom)
{
    return compat && compat->fb ?
	imgstream_fb_getview(compat->fb, xcenter, ycenter, xzoom, yzoom) : -1;
}


int
fb_imgstream_compat_cursor(struct fb *ifp,
	struct fb_imgstream_compat *compat,
	int mode,
	int x,
	int y)
{
    if (!ifp || !ifp->i || !compat || !compat->fb)
	return -1;

    int ret = imgstream_fb_cursor(compat->fb, mode, x, y);
    if (ret == 0) {
	ifp->i->if_cursmode = mode;
	ifp->i->if_xcurs = x;
	ifp->i->if_ycurs = y;
    }
    return ret;
}


int
fb_imgstream_compat_getcursor(const struct fb_imgstream_compat *compat,
	int *mode,
	int *x,
	int *y)
{
    return compat && compat->fb ?
	imgstream_fb_getcursor(compat->fb, mode, x, y) : -1;
}


int
fb_imgstream_compat_setcursor(struct fb_imgstream_compat *compat,
	const unsigned char *bits,
	int xbits,
	int ybits,
	int xorig,
	int yorig)
{
    return compat && compat->fb ?
	imgstream_fb_setcursor(compat->fb, bits, xbits, ybits, xorig, yorig) : -1;
}


int
fb_imgstream_compat_rmap(const struct fb_imgstream_compat *compat,
	ColorMap *cmap)
{
    if (!compat || !compat->fb || !cmap)
	return -1;

    struct imgstream_fb_colormap imap;
    if (imgstream_fb_rmap(compat->fb, &imap) != 0)
	return -1;

    for (int i = 0; i < 256; i++) {
	cmap->cm_red[i] = imap.red[i];
	cmap->cm_green[i] = imap.green[i];
	cmap->cm_blue[i] = imap.blue[i];
    }
    return 0;
}


int
fb_imgstream_compat_wmap(struct fb_imgstream_compat *compat,
	const ColorMap *cmap)
{
    if (!compat || !compat->fb)
	return -1;

    if (!cmap)
	return imgstream_fb_wmap(compat->fb, NULL);

    struct imgstream_fb_colormap imap;
    for (int i = 0; i < 256; i++) {
	imap.red[i] = cmap->cm_red[i];
	imap.green[i] = cmap->cm_green[i];
	imap.blue[i] = cmap->cm_blue[i];
    }
    return imgstream_fb_wmap(compat->fb, &imap);
}
