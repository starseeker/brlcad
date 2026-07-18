/*                        F B G R I D . C
 * BRL-CAD
 *
 * Copyright (c) 1986-2026 United States Government as represented by
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
 *
 */
/** @file fbgrid.c
 *
 */

#include "common.h"

#include <errno.h>
#include <limits.h>
#include <stdlib.h>

#include "bio.h"

#include "bu/app.h"
#include "bu/color.h"
#include "bu/getopt.h"
#include "bu/exit.h"
#include "imgstream/fb_compat.h"


static unsigned char *white_line, *grey_line, *dark_line;
static imgstream_fb_t *fbp;
static char *framebuffer = NULL;

#define OLD 0
#define BINARY 1
#define DECIMAL 2

static int fbwidth = 0;
static int fbheight = 0;
static int flavor = DECIMAL;
static int clear = 0;

void grid(imgstream_fb_t *fbiop, unsigned char *line, int spacing), oldflavor(void);

static char usage[] = "\
Usage: fbgrid [-c] [-b | -d | -o] [-F framebuffer]\n\
	[-S squaresize] [-W width] [-N height]\n";

static int
parse_positive_int_arg(const char *label, const char *value, int *result)
{
    long parsed;
    char *endptr = NULL;

    errno = 0;
    parsed = strtol(value, &endptr, 10);
    if (errno != 0 || endptr == value || *endptr != '\0' || parsed <= 0 || parsed > INT_MAX) {
	fprintf(stderr, "fbgrid: invalid %s '%s'\n", label, value);
	return 0;
    }

    *result = (int)parsed;
    return 1;
}

int
get_args(int argc, char **argv)
{
    int c;
    int remaining = 0;

    while ((c = bu_getopt(argc, argv, "cbdoF:s:w:n:S:W:N:h?")) != -1) {
	switch (c) {
	    case 'c':
		clear = 1;
		break;
	    case 'b':
		flavor = BINARY;
		break;
	    case 'd':
		flavor = DECIMAL;
		break;
	    case 'o':
		flavor = OLD;
		break;
	    case 'F':
		framebuffer = bu_optarg;
		break;
	    case 'S':
	    case 's':
		/* square size */
		if (!parse_positive_int_arg("size", bu_optarg, &fbwidth))
		    return 0;
		fbheight = fbwidth;
		break;
	    case 'W':
	    case 'w':
		if (!parse_positive_int_arg("width", bu_optarg, &fbwidth))
		    return 0;
		break;
	    case 'N':
	    case 'n':
		if (!parse_positive_int_arg("height", bu_optarg, &fbheight))
		    return 0;
		break;

	    default:		/* '?' */
		return 0;
	}
    }

    if (argc == 1 && isatty(fileno(stdin)) && isatty(fileno(stdout)))
	return 0;

    remaining = argc - bu_optind;
    if (remaining != 0) {
	fprintf(stderr, "fbgrid: excess argument(s) not supported\n");
	return 0;
    }

    return 1;		/* OK */
}


int
main(int argc, char **argv)
{
    int i;

    bu_setprogname(argv[0]);

    if (!get_args(argc, argv)) {
	(void)fputs(usage, stderr);
	bu_exit(1, NULL);
    }

    if (flavor == OLD)
	oldflavor();	/* exits */

    if ((fbp = imgstream_fb_open(framebuffer, (size_t)fbwidth,
				 (size_t)fbheight)) == NULL)
	bu_exit(2, NULL);

    fbwidth = (int)imgstream_fb_width(fbp);
    fbheight = (int)imgstream_fb_height(fbp);

    /* Initialize the color lines */
    white_line = (unsigned char *)malloc((size_t)fbwidth * 3);
    grey_line  = (unsigned char *)malloc((size_t)fbwidth * 3);
    dark_line  = (unsigned char *)malloc((size_t)fbwidth * 3);
    for (i = 0; i < fbwidth; i++) {
	white_line[3*i+RED] = white_line[3*i+GRN] = white_line[3*i+BLU] = 255;
	grey_line[3*i+RED] = grey_line[3*i+GRN] = grey_line[3*i+BLU] = 128;
	dark_line[3*i+RED] = dark_line[3*i+GRN] = dark_line[3*i+BLU] = 64;
    }

    if (clear)
	imgstream_fb_clear(fbp, NULL);

    if (flavor == BINARY) {
	/* Dark lines every 8 */
	grid(fbp, dark_line, 8);
	/* Grey lines every 64 */
	grid(fbp, grey_line, 64);
	/* White line every 128 */
	grid(fbp, white_line, 128);
    } else {
	/* DECIMAL */
	/* Dark lines every 10 */
	grid(fbp, dark_line, 10);
	/* Grey lines every 50 */
	grid(fbp, grey_line, 50);
	/* White line every 100 */
	grid(fbp, white_line, 100);
    }

    imgstream_fb_close(fbp);
    return 0;
}


void
grid(imgstream_fb_t *fbiop, unsigned char *line, int spacing)
{
    int x, y;

    for (y = 0; y < fbheight; y += spacing)
	imgstream_fb_write(fbiop, 0, y, line, (size_t)fbwidth);
    for (x = 0; x < fbwidth; x += spacing) {
	imgstream_fb_writerect(fbiop, x, 0, 1, fbheight, line);
    }
}


void
oldflavor(void)
{
    imgstream_fb_t *fbiop;
    int x, y;
    int middle;
    int mask;
    int fb_sz;
    unsigned char *line;
    static unsigned char black[3] = {0, 0, 0};
    static unsigned char white[3] = {255, 255, 255};
    static unsigned char red[3] = {255, 0, 0};

    fbiop = imgstream_fb_open(framebuffer, (size_t)fbwidth,
			     (size_t)fbheight);
    if (fbiop == NULL) {
	bu_exit(1, NULL);
    }

    fb_sz = (int)imgstream_fb_width(fbiop);
    middle = fb_sz/2;
    line = (unsigned char *)malloc((size_t)fb_sz * 3);
    if (line == NULL)
	bu_exit(1, "fbgrid: unable to allocate scanline\n");
    if (fb_sz <= 512)
	mask = 0x7;
    else
	mask = 0xf;

    for (y = fb_sz-1; y >= 0; y--) {
	for (x = 0; x < fb_sz; x++) {
	    const unsigned char *pixel;
	    if (x == y || x == fb_sz - y) {
		pixel = white;
	    } else
		if (x == middle || y == middle) {
		    pixel = red;
		} else
		    if ((x & mask) && (y & mask)) {
			pixel = black;
		    } else {
			pixel = white;
		    }
	    line[3*x+RED] = pixel[RED];
	    line[3*x+GRN] = pixel[GRN];
	    line[3*x+BLU] = pixel[BLU];
	}
	if (imgstream_fb_write(fbiop, 0, y, line, (size_t)fb_sz) != fb_sz)
	    bu_exit(1, "fbgrid: framebuffer write failed\n");
    }
    free(line);
    imgstream_fb_close(fbiop);
    bu_exit(0, NULL);
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
