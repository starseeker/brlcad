/*                         D D I S P . C
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file ddisp.c
 *
 * Data Display - shows doubles on a framebuffer in various ways.
 *
 */

#include "common.h"

#include <string.h>
#include <stdlib.h>

#include "bio.h"

#include "bu/app.h"
#include "bu/malloc.h"
#include "bu/color.h"
#include "bu/getopt.h"
#include "bu/str.h"
#include "bu/exit.h"
#include "bu/snooze.h"
#include "imgstream/fb_compat.h"


#define VERT 1
#define BARS 2


static void
lineout(imgstream_fb_t *fbp, double *dat, int n)
{
    static int y = 0;
    int i, value;
    unsigned char lbuf[1024*4][3];

    if (n > (int)imgstream_fb_width(fbp)) n = (int)imgstream_fb_width(fbp);

    for (i = 0; i < n; i++) {
	/* Magnitude version */
	value = dat[i] * 255.9;
	if (value < 0) value = 0;
	else if (value > 255) value = 255;
	lbuf[i][RED] = lbuf[i][GRN] = lbuf[i][BLU] = value;
    }
    imgstream_fb_write(fbp, 0, y, (unsigned char *)lbuf, (size_t)n);

    /* Next screen position */
    y = (y + 1) % (int)imgstream_fb_height(fbp);
}


/*
 * Display doubles.
 * +/- 1.0 in, becomes +/- 128 from center Y.
 */
static void
disp_inten(imgstream_fb_t *fbp, double *buf, int size)
{
    int x, y;
    unsigned char color[3];

/* color.red = color.green = color.blue = 255;*/

    if (size > (int)imgstream_fb_width(fbp)) size = (int)imgstream_fb_width(fbp);

    for (x = 0; x < size; x++) {
	y = buf[x] * 128;
#ifdef OVERLAY
	imgstream_fb_read(fbp, x, y+255, color, 1);
#else
	color[RED] = color[BLU] = 0;
#endif
	color[GRN] = 255;
	imgstream_fb_write(fbp, x, y+255, color, 1);
    }
}


/*
 * Display doubles.
 * +/- 1.0 in, becomes +/- 128 from center Y.
 */
static void
disp_bars(imgstream_fb_t *fbp, double *buf, int size)
{
    int x, y;
    unsigned char color[3];

/* color.red = color.green = color.blue = 255;*/

    if (size > (int)imgstream_fb_width(fbp)) size = (int)imgstream_fb_width(fbp);

    for (x = 0; x < size; x++) {
	if (buf[x] > 1.0) {
	    y = 128;
	} else if (buf[x] < -1.0) {
	    y = -128;
	} else {
	    y = buf[x] * 128;
	}
#ifdef OVERLAY
	imgstream_fb_read(fbp, x, y+255, color, 1);
#else
	color[RED] = color[BLU] = 0;
#endif
	color[GRN] = 255;
	if (y > 0) {
	    while (y >= 0) {
		imgstream_fb_write(fbp, x, y+255, color, 1);
		y--;
	    }
	} else {
	    while (y <= 0) {
		imgstream_fb_write(fbp, x, y+255, color, 1);
		y++;
	    }
	}
    }
}


int
main(int argc, char **argv)
{
    static const char usage[] = "Usage: ddisp [-F framebuffer] [-v -b -p -c] [width] < inputfile\n";

    imgstream_fb_t *fbp = NULL;
    double buf[BU_PAGE_SIZE];

    int n, L;
    int Clear = 0;
    int pause_time = 0;
    int mode = 0;
    int fbsize = 512;
    const char *framebuffer = NULL;
    int c;

    bu_setprogname(argv[0]);

    if (isatty(fileno(stdin))) {
	bu_exit(1, "%s", usage);
    }


    while ((c = bu_getopt(argc, argv, "F:vbpcHh?")) != -1) {
	if (c == 'F') {
	    framebuffer = bu_optarg;
	} else if (c == 'v') {
	    mode = VERT;
	    pause_time = 0;
	    Clear = 0;
	} else if (c == 'b') {
	    mode = BARS;
	} else if (c == 'p') {
	    pause_time = 3;
	} else if (c == 'c') {
	    Clear++;
	} else {
	    bu_exit(1, "%s", usage);
	}
    }

    if (bu_optind < argc) {
	L = atoi(argv[bu_optind++]);
	if (L <= 0 || L > BU_PAGE_SIZE)
	    bu_exit(1, "ddisp: width must be between 1 and %d\n", BU_PAGE_SIZE);
    } else {
	L = BU_PAGE_SIZE;
    }
    if (bu_optind < argc)
	bu_exit(1, "%s", usage);

    if ((fbp = imgstream_fb_open(framebuffer, (size_t)fbsize,
				 (size_t)fbsize)) == NULL) {
	bu_exit(2, "Unable to open framebuffer\n");
    }

    while ((n = fread(buf, sizeof(*buf), L, stdin)) > 0) {
	/* XXX - width hack */
	if (n > (int)imgstream_fb_width(fbp))
	    n = (int)imgstream_fb_width(fbp);

	if (Clear)
	    imgstream_fb_clear(fbp, NULL);
	if (mode == VERT)
	    disp_inten(fbp, buf, n);
	else if (mode == BARS)
	    disp_bars(fbp, buf, n);
	else
	    lineout(fbp, buf, n);
	if (pause_time)
	    bu_snooze(BU_SEC2USEC(pause_time));
    }
    imgstream_fb_close(fbp);

    return 0;
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
