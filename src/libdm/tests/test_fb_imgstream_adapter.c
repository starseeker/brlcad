/*              T E S T _ F B _ I M G S T R E A M _ A D A P T E R . C
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

#include <stdio.h>
#include <string.h>

#include "bu/app.h"
#include "dm.h"
#include "../include/private.h"


#define CHECK(_condition, _message) do { \
    if (!(_condition)) { \
	fprintf(stderr, "%s\n", (_message)); \
	return 1; \
    } \
} while (0)


int
main(int argc, char **argv)
{
    static const unsigned char pixels[12] = {
	1, 2, 3, 4, 5, 6,
	7, 8, 9, 10, 11, 12
    };
    static const unsigned char clipped_pixels[6] = {13, 14, 15, 16, 17, 18};
    static const unsigned char corner_pixels[12] = {
	21, 22, 23, 24, 25, 26,
	27, 28, 29, 30, 31, 32
    };
    unsigned char readback[12] = {0};
    ColorMap cmap;
    int mode = 0;
    int x = 0;
    int y = 0;
    int xcenter = 0;
    int ycenter = 0;
    int xzoom = 0;
    int yzoom = 0;

    bu_setprogname(argv[0]);
    (void)argc;

    struct fb *fbp = fb_open("/dev/mem", 2, 2);
    CHECK(fbp != FB_NULL, "failed to open memory framebuffer");
    CHECK(fbp->i->pp != NULL,
	"memory framebuffer did not select the image stream adapter");
    CHECK(fb_getwidth(fbp) == 2 && fb_getheight(fbp) == 2,
	"memory framebuffer dimensions changed");
    CHECK(fb_writerect(fbp, 0, 0, 2, 2, pixels) == 4,
	"rectangle write failed");
    CHECK(fb_readrect(fbp, 0, 0, 2, 2, readback) == 4 &&
	memcmp(pixels, readback, sizeof(pixels)) == 0,
	"rectangle readback changed");
    CHECK(fb_writerect(fbp, 1, 0, 2, 1, clipped_pixels) == 1,
	"right-edge rectangle clipping changed");
    CHECK(fb_writerect(fbp, -1, -1, 2, 2, corner_pixels) == 1 &&
	fb_read(fbp, 0, 0, readback, 1) == 1 &&
	memcmp(readback, corner_pixels + 9, 3) == 0,
	"negative-coordinate rectangle clipping changed");
    CHECK(fb_view(fbp, 1, 1, 2, 3) == 0 &&
	fb_getview(fbp, &xcenter, &ycenter, &xzoom, &yzoom) == 0 &&
	xcenter == 1 && ycenter == 1 && xzoom == 2 && yzoom == 3,
	"view state changed");
    CHECK(fb_cursor(fbp, 1, 1, 0) == 0 &&
	fb_getcursor(fbp, &mode, &x, &y) == 0 &&
	mode == 1 && x == 1 && y == 0,
	"cursor state changed");
    fb_make_linear_cmap(&cmap);
    cmap.cm_red[1] = 1234;
    CHECK(fb_wmap(fbp, &cmap) == 0, "colormap write failed");
    memset(&cmap, 0, sizeof(cmap));
    CHECK(fb_rmap(fbp, &cmap) == 0 && cmap.cm_red[1] == 1234,
	"colormap readback changed");
    CHECK(fb_close(fbp) == 0, "memory framebuffer close failed");

    fbp = fb_open("/dev/memw /dev/null", 1, 1);
    CHECK(fbp != FB_NULL, "failed to open attached legacy memory framebuffer");
    CHECK(fbp->i->pp == NULL,
	"attached memory framebuffer did not retain the legacy backend");
    CHECK(fb_close(fbp) == 0, "attached memory framebuffer close failed");
    return 0;
}
