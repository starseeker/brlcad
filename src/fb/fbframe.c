/*                       F B F R A M E . C
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
/** @file fbframe.c
 *
 * Overwrite a frame (border) on the current framebuffer.  CCW from
 * the bottom: Red, Green, Blue, White
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
#include "bu/malloc.h"
#include "imgstream/fb_compat.h"

const char *Usage="[-F framebuffer] [-s|S squareframesize] [-w|W frame_width] [-n|N frame_height]\n";

#define USAGE_EXIT(p) { fprintf(stderr, "Usage: %s %s\n", (p), Usage); \
	bu_exit(-1, NULL); }

static int
parse_positive_int_arg(const char *arg, int *value, const char *label)
{
    char *end = NULL;
    long parsed = 0;

    errno = 0;
    parsed = strtol(arg, &end, 10);
    if (arg[0] == '\0' || end == arg || *end != '\0' || errno != 0) {
	fprintf(stderr, "%s: invalid %s '%s'\n", bu_getprogname(), label, arg);
	return 0;
    }
    if (parsed <= 0) {
	fprintf(stderr, "%s: %s must be greater than zero, got '%s'\n", bu_getprogname(), label, arg);
	return 0;
    }
    if (parsed > INT_MAX) {
	fprintf(stderr, "%s: %s out of range '%s'\n", bu_getprogname(), label, arg);
	return 0;
    }

    *value = (int)parsed;
    return 1;
}

int
main(int argc, char **argv)
{
    int c;
    int x;
    imgstream_fb_t *fbp;
    int xsize, ysize;
    int len;
    char *framebuffer = (char *)NULL;
    unsigned char *line;
    static unsigned char white[3] = { 255, 255, 255 };
    static unsigned char red[3] = { 255, 0, 0 };
    static unsigned char green[3] = { 0, 255, 0 };
    static unsigned char blue[3] = { 0, 0, 255 };

    bu_setprogname(argv[0]);

    xsize = ysize = 0;
    while ((c = bu_getopt(argc, argv, "F:s:w:n:S:W:N:h?")) != -1) {
	switch (c) {
	    case 'F':
		framebuffer = bu_optarg;
		break;
	    case 's':
	    case 'S':
		/* square file size */
		if (parse_positive_int_arg(bu_optarg, &len, "square frame size"))
		    xsize = ysize = len;
		else
		    USAGE_EXIT(*argv);
		break;
	    case 'w':
	    case 'W':
		if (parse_positive_int_arg(bu_optarg, &len, "frame width"))
		    xsize = len;
		else
		    USAGE_EXIT(*argv);
		break;
	    case 'n':
	    case 'N':
		if (parse_positive_int_arg(bu_optarg, &len, "frame height"))
		    ysize = len;
		else
		    USAGE_EXIT(*argv);
		break;
	    default:	/* '?' 'h' */
		USAGE_EXIT(*argv);
		break;
	}
    }

    if (argc == 1 && isatty(fileno(stdin)) && isatty(fileno(stdout)))
	USAGE_EXIT(*argv);
    if (argc > bu_optind) {
	fprintf(stderr, "%s: excess argument(s) not supported\n", bu_getprogname());
	USAGE_EXIT(*argv);
    }

    if ((fbp = imgstream_fb_open(framebuffer, (size_t)xsize,
	    (size_t)ysize)) == NULL)
	bu_exit(1, NULL);

    if (xsize <= 0)
	xsize = (int)imgstream_fb_width(fbp);
    if (ysize <= 0)
	ysize = (int)imgstream_fb_height(fbp);

    /* malloc buffer for pixel lines */
    len = (xsize > ysize) ? xsize : ysize;
    line = (unsigned char *)bu_calloc((size_t)len, 3, "line");
    if (!line) {
	fprintf(stderr, "fbframe:  malloc failure\n");
	return 1;
    }

#define FLOOD(col) { for (x=len-1; x >= 0; x--) { \
    line[3*x+RED] = (col)[RED]; line[3*x+GRN] = (col)[GRN]; \
    line[3*x+BLU] = (col)[BLU]; } }

    /*
     * Red:	(0->510,      0)
     * Green:	(511, 0->510)
     * Blue:	(511->1,    511)
     * White:	(0, 511->1)
     */
    FLOOD(red);
    imgstream_fb_writerect(fbp, 0, 0, xsize-1, 1, line);
    FLOOD(green);
    imgstream_fb_writerect(fbp, xsize-1, 0, 1, ysize-1, line);
    FLOOD(blue);
    imgstream_fb_writerect(fbp, 1, ysize-1, xsize-1, 1, line);
    FLOOD(white);
    imgstream_fb_writerect(fbp, 0, 1, 1, ysize-1, line);

    imgstream_fb_close(fbp);

    bu_free(line, "line");

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
