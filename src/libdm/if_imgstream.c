/*                  I F _ I M G S T R E A M . C
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

#include "common.h"

#include "imgstream/fb_compat.h"
#include "./include/private.h"
#include "dm.h"


#define IMGSTREAM_FB(_ifp) ((imgstream_fb_t *)(_ifp)->i->pp)


static int
if_imgstream_open(struct fb *ifp, const char *file, int width, int height)
{
    if (!ifp || !file || width < 0 || height < 0)
	return -1;

    imgstream_fb_t *fb = imgstream_fb_open(file, (size_t)width, (size_t)height);
    if (!fb)
	return -1;

    size_t actual_width = imgstream_fb_width(fb);
    size_t actual_height = imgstream_fb_height(fb);
    if (actual_width > INT_MAX || actual_height > INT_MAX) {
	imgstream_fb_close(fb);
	return -1;
    }

    ifp->i->pp = fb;
    ifp->i->if_width = (int)actual_width;
    ifp->i->if_height = (int)actual_height;
    ifp->i->if_xzoom = 1;
    ifp->i->if_yzoom = 1;
    ifp->i->if_xcenter = ifp->i->if_width / 2;
    ifp->i->if_ycenter = ifp->i->if_height / 2;
    ifp->i->if_cursmode = 0;
    ifp->i->if_xcurs = 0;
    ifp->i->if_ycurs = 0;
    return 0;
}


static int
if_imgstream_open_existing(struct fb *UNUSED(ifp), int UNUSED(width),
	int UNUSED(height), struct fb_platform_specific *UNUSED(fb_p))
{
    return -1;
}


static int
if_imgstream_close_existing(struct fb *UNUSED(ifp))
{
    return -1;
}


static struct fb_platform_specific *
if_imgstream_existing_get(uint32_t UNUSED(magic))
{
    return NULL;
}


static void
if_imgstream_existing_put(struct fb_platform_specific *UNUSED(fb_p))
{
}


static int
if_imgstream_close(struct fb *ifp)
{
    if (ifp->i->dmp && ifp->i->dmp->i->fbp == ifp)
	ifp->i->dmp->i->fbp = FB_NULL;
    imgstream_fb_close(IMGSTREAM_FB(ifp));
    ifp->i->pp = NULL;
    return 0;
}


static int
if_imgstream_clear(struct fb *ifp, unsigned char *rgb)
{
    return imgstream_fb_clear(IMGSTREAM_FB(ifp), rgb);
}


static ssize_t
if_imgstream_read(struct fb *ifp, int x, int y, unsigned char *rgb,
	size_t count)
{
    return imgstream_fb_read(IMGSTREAM_FB(ifp), x, y, rgb, count);
}


static ssize_t
if_imgstream_write(struct fb *ifp, int x, int y, const unsigned char *rgb,
	size_t count)
{
    return imgstream_fb_write(IMGSTREAM_FB(ifp), x, y, rgb, count);
}


static int
if_imgstream_rmap(struct fb *ifp, ColorMap *cmap)
{
    struct imgstream_fb_colormap imgstream_cmap;
    int ret = imgstream_fb_rmap(IMGSTREAM_FB(ifp), &imgstream_cmap);
    if (ret != 0)
	return ret;

    for (size_t i = 0; i < 256; i++) {
	cmap->cm_red[i] = imgstream_cmap.red[i];
	cmap->cm_green[i] = imgstream_cmap.green[i];
	cmap->cm_blue[i] = imgstream_cmap.blue[i];
    }
    return 0;
}


static int
if_imgstream_wmap(struct fb *ifp, const ColorMap *cmap)
{
    if (!cmap)
	return imgstream_fb_wmap(IMGSTREAM_FB(ifp), NULL);

    struct imgstream_fb_colormap imgstream_cmap;
    for (size_t i = 0; i < 256; i++) {
	imgstream_cmap.red[i] = cmap->cm_red[i];
	imgstream_cmap.green[i] = cmap->cm_green[i];
	imgstream_cmap.blue[i] = cmap->cm_blue[i];
    }
    return imgstream_fb_wmap(IMGSTREAM_FB(ifp), &imgstream_cmap);
}

static int
if_imgstream_view(struct fb *ifp, int xcenter, int ycenter, int xzoom,
	int yzoom)
{
    int ret = imgstream_fb_view(IMGSTREAM_FB(ifp), xcenter, ycenter, xzoom,
	yzoom);
    if (ret == 0) {
	ifp->i->if_xcenter = xcenter;
	ifp->i->if_ycenter = ycenter;
	ifp->i->if_xzoom = xzoom;
	ifp->i->if_yzoom = yzoom;
    }
    return ret;
}


static int
if_imgstream_getview(struct fb *ifp, int *xcenter, int *ycenter, int *xzoom,
	int *yzoom)
{
    return imgstream_fb_getview(IMGSTREAM_FB(ifp), xcenter, ycenter, xzoom,
	yzoom);
}


static int
if_imgstream_setcursor(struct fb *ifp, const unsigned char *bits, int xbits,
	int ybits, int xorig, int yorig)
{
    return imgstream_fb_setcursor(IMGSTREAM_FB(ifp), bits, xbits, ybits,
	xorig, yorig);
}


static int
if_imgstream_cursor(struct fb *ifp, int mode, int x, int y)
{
    int ret = imgstream_fb_cursor(IMGSTREAM_FB(ifp), mode, x, y);
    if (ret == 0) {
	ifp->i->if_cursmode = mode;
	ifp->i->if_xcurs = x;
	ifp->i->if_ycurs = y;
    }
    return ret;
}


static int
if_imgstream_getcursor(struct fb *ifp, int *mode, int *x, int *y)
{
    return imgstream_fb_getcursor(IMGSTREAM_FB(ifp), mode, x, y);
}


static int
if_imgstream_configure_window(struct fb *UNUSED(ifp), int UNUSED(width),
	int UNUSED(height))
{
    return 0;
}


static int
if_imgstream_refresh(struct fb *UNUSED(ifp), int UNUSED(x), int UNUSED(y),
	int UNUSED(width), int UNUSED(height))
{
    return 0;
}


static int
if_imgstream_poll(struct fb *ifp)
{
    return imgstream_fb_poll(IMGSTREAM_FB(ifp));
}


static int
if_imgstream_flush(struct fb *ifp)
{
    return imgstream_fb_flush(IMGSTREAM_FB(ifp));
}


static int
if_imgstream_help(struct fb *ifp)
{
    fb_log("Description: %s\n", ifp->i->if_type);
    fb_log("Device: %s\n", ifp->i->if_name);
    fb_log("Max width/height: %d %d\n", ifp->i->if_max_width,
	ifp->i->if_max_height);
    fb_log("Default width/height: %d %d\n", ifp->i->if_width,
	ifp->i->if_height);
    return 0;
}

struct fb_impl imgstream_memory_interface_impl = {
    0,
    FB_MEMORY_MAGIC,
    if_imgstream_open,
    if_imgstream_open_existing,
    if_imgstream_close_existing,
    if_imgstream_existing_get,
    if_imgstream_existing_put,
    if_imgstream_close,
    if_imgstream_clear,
    if_imgstream_read,
    if_imgstream_write,
    if_imgstream_rmap,
    if_imgstream_wmap,
    if_imgstream_view,
    if_imgstream_getview,
    if_imgstream_setcursor,
    if_imgstream_cursor,
    if_imgstream_getcursor,
    fb_sim_readrect,
    fb_sim_writerect,
    fb_sim_bwreadrect,
    fb_sim_bwwriterect,
    if_imgstream_configure_window,
    if_imgstream_refresh,
    if_imgstream_poll,
    if_imgstream_flush,
    if_imgstream_close,
    if_imgstream_help,
    "Memory Buffer",
    FB_XMAXSCREEN,
    FB_YMAXSCREEN,
    "/dev/mem",
    512,
    512,
    -1,
    -1,
    1, 1,
    256, 256,
    0, 0, 0,
    PIXEL_NULL,
    PIXEL_NULL,
    PIXEL_NULL,
    -1,
    0,
    0L,
    0L,
    0,
    0,
    NULL,
    NULL,
    0,
    NULL,
    {0},
    {0},
    {0},
    {0},
    {0},
    {0},
    0,
    {{FB_EVENT_NONE, 0, 0, 0, 0, 0}},
    0,
    0
};


struct fb imgstream_memory_interface = { &imgstream_memory_interface_impl };

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
