/*                        B W H I S T . C
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
 */
/** @file util/bwhist.c
 *
 * Display, and optionally dump to tty, a histogram of a
 * black and white file.  Black is top of screen, white bottom.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "bio.h"

#include "bu/app.h"
#include "bu/getopt.h"
#include "bu/str.h"
#include "bu/exit.h"
#include "imgstream/fb_compat.h"


long bin[256];
int verbose = 0;
imgstream_fb_t *fbp;

static const char *Usage = "Usage: bwhist [-F framebuffer] [-v] [file.bw]\n";


int
main(int argc, char **argv)
{
    size_t i;
    size_t n;
    long max;
    static double scale;
    unsigned char buf[512];
    unsigned char *bp;
    unsigned char white[3*512];
    FILE *fp;
    const char *framebuffer = NULL;
    int c;

    bu_setprogname(argv[0]);
    setmode(fileno(stdin), O_BINARY);
    setmode(fileno(stdout), O_BINARY);

    while ((c = bu_getopt(argc, argv, "F:vh?")) != -1) {
	switch (c) {
	    case 'F': framebuffer = bu_optarg; break;
	    case 'v': verbose++; break;
	    default: bu_exit(1, "%s", Usage);
	}
    }

    /* look for optional input file */
    if (bu_optind < argc) {
	if ((fp = fopen(argv[bu_optind], "rb")) == 0) {
	    bu_exit(1, "bwhist: can't open '%s'\n", argv[bu_optind]);
	}
	bu_optind++;
    } else {
	fp = stdin;
    }
    /* check usage */
    if (bu_optind < argc || isatty(fileno(fp)))
	bu_exit(1, "%s", Usage);

    for (i = 0; i < 3*512; i++)
	white[i] = 255;

    while ((n = fread(buf, sizeof(*buf), sizeof(buf), fp)) > 0) {
	bp = &buf[0];
	for (i = 0; i < n; i++)
	    bin[ *bp++ ]++;
    }

    /* find max */
    max = 1;
    for (i = 0; i < 256; i++)
	if (bin[i] > max) max = bin[i];
    scale = 511.0 / (double)max;

    /* Display the max? */
    printf("Full screen = %ld pixels\n", max);

    if ((fbp = imgstream_fb_open(framebuffer, 512, 512)) == NULL) {
	bu_exit(12, "fb_open failed\n");
    }

    /* Display them */
    for (i = 0; i < 256; i++) {
	int value;
	value = bin[i]*scale;
	if (value == 0 && bin[i] != 0) value = 1;
	imgstream_fb_write(fbp, 0, (int)(2*i), white, (size_t)value);
	imgstream_fb_write(fbp, 0, (int)(2*i+1), white, (size_t)value);
	if (verbose)
	    printf("%3lu: %10ld (%10f)\n", (long unsigned)i, bin[i], (float)bin[i]/(float)max);
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
