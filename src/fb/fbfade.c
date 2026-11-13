/*                        F B F A D E . C
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
 *
 */
/** @file fbfade.c
 *
 * fade in or out a frame buffer image
 *
 * This program displays a frame buffer image gradually, randomly
 * selecting the pixel display sequence.  (Suggested by Gary Moss.)
 * It requires fast single-pixel write support for best effect.
 *
 * Options:
 *
 * "-H" assumes 1024x1024 default input size instead of 512x512
 *
 * "-f in_fb_file" reads from the specified frame buffer file instead
 * of assuming constant black ("fade out") value
 *
 * "-s size" is the input size (width & height)
 *
 * "-w width" is the input width
 *
 * "-n height" is the input height
 *
 * "-F out_fb_file" writes to the specified frame buffer file instead of
 * the one specified by the FB_FILE environment variable (the default
 * frame buffer, if no FB_FILE)
 *
 * "-S size" is the output size (width & height)
 *
 * "-W width" is the output width
 *
 * "-N height" is the output height
 *
 * "out_fb_file" is the same as -F out_fb_file, for convenience
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "bio.h"

#include "bu/app.h"
#include "bu/getopt.h"
#include "bu/interrupt.h"
#include "vmath.h"
#include "imgstream/fb_compat.h"
#include "pkg.h"


#define USAGE1 "Usage: fbfade [ -s size ] [ -w width ] [ -n height ] [ -f in_fb_file ]\n\
[ -H ] [ -S size ] [ -W width ] [ -N height ] [ [ -F ] out_fb_file ]"
#define OPTSTR "f:F:Hn:N:s:S:w:W:h?"


typedef int bool_t;

static bool_t hires = 0;		/* set for 1024x1024; clear for 512x512 */
static char *in_fb_file = NULL;		/* input image name */
static char *out_fb_file = NULL;	/* output frame buffer name */
static imgstream_fb_t *fbp = NULL;	/* framebuffer input/output handle */
static int src_width = 0;		/* input image width */
static int src_height = 0;		/* input image height */
static int dst_width = 0;		/* output frame buffer size */
static int dst_height = 0;		/* output frame buffer size */
static unsigned char (*pix)[3];		/* input image */
static unsigned char bg[3] = { 0, 0, 0 };	/* background */

/* in ioutil.c */
extern void Message(const char *format, ...);
extern void Fatal(imgstream_fb_t *fbiop, const char *format, ...);


static void
Sig_Catcher(int sig)
{
    (void)signal(sig, SIG_DFL);

    /* The following is not guaranteed to work, but it's worth a try. */
    Fatal(fbp, "Interrupted by signal %d", sig);
}


int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);

    /* Plant signal catcher. */
    {
	/* signals to catch */
	static int getsigs[] = {
#ifdef SIGHUP
	    SIGHUP,			/* hangup */
#endif
#ifdef SIGINT
	    SIGINT,			/* interrupt */
#endif
#ifdef SIGQUIT
	    SIGQUIT,		/* quit */
#endif
#ifdef SIGPIPE
	    SIGPIPE,		/* write on a broken pipe */
#endif
#ifdef SIGTERM
	    SIGTERM,		/* software termination signal */
#endif
	    0
	};
	int i;

	for (i = 0; getsigs[i] != 0; ++i)
	    if (signal(getsigs[i], SIG_IGN) != SIG_IGN)
		(void)signal(getsigs[i], Sig_Catcher);
    }

    /* Process arguments. */
    {
	int c;
	bool_t errors = 0;

	while ((c = bu_getopt(argc, argv, OPTSTR)) != -1)
	    switch (c) {
		default:	/* '?': invalid option */
		    errors = 1;
		    break;

		case 'f':	/* -f in_fb_file */
		    in_fb_file = bu_optarg;
		    break;

		case 'F':	/* -F out_fb_file */
		    out_fb_file = bu_optarg;
		    break;

		case 'H':	/* -H */
		    hires = 1;
		    break;

		case 'n':	/* -n height */
		    if ((src_height = atoi(bu_optarg)) <= 0)
			errors = 1;

		    break;

		case 'N':	/* -N height */
		    if ((dst_height = atoi(bu_optarg)) <= 0)
			errors = 1;

		    break;

		case 's':	/* -s size */
		    if ((src_height = src_width = atoi(bu_optarg)) <= 0)
			errors = 1;

		    break;

		case 'S':	/* -S size */
		    if ((dst_height = dst_width = atoi(bu_optarg)) <= 0)
			errors = 1;

		    break;

		case 'w':	/* -w width */
		    if ((src_width = atoi(bu_optarg)) <= 0)
			errors = 1;

		    break;

		case 'W':	/* -W width */
		    if ((dst_width = atoi(bu_optarg)) <= 0)
			errors = 1;

		    break;
	    }

	if (argc == 1 && isatty(fileno(stdin)) && isatty(fileno(stdout)))
    	    errors = 1;
	if (errors)
	    Fatal(fbp, USAGE1);
    }

    if (bu_optind < argc) {
	/* out_fb_file */
	if (bu_optind < argc - 1 || out_fb_file != NULL) {
	    Message(USAGE1);
	    Fatal(fbp, "Can't handle multiple output frame buffers!");
	}

	out_fb_file = argv[bu_optind];
    }

    /* Open frame buffer for unbuffered input. */

    if (src_width == 0)
	src_width = hires ? 1024 : 512;		/* starting default */

    if (src_height == 0)
	src_height = hires ? 1024 : 512;	/* starting default */

    if (in_fb_file != NULL) {

	if ((fbp = imgstream_fb_open(in_fb_file, (size_t)src_width,
				      (size_t)src_height)) == NULL)
	    Fatal(fbp, "Couldn't open input frame buffer");
	else {
	    int y;
	    int wt = (int)imgstream_fb_width(fbp);
	    int ht = (int)imgstream_fb_height(fbp);

	    /* Use smaller actual input size instead of request. */

	    V_MIN(src_width, wt);
	    V_MIN(src_height, ht);

	    if ((pix = (unsigned char (*)[3])malloc((size_t)src_width * (size_t)src_height * sizeof(*pix))) == NULL)
		Fatal(fbp, "Not enough memory for pixel array");

	    for (y = 0; y < src_height; ++y) {
		if (imgstream_fb_read(fbp, 0, y, pix[y * src_width],
				      (size_t)src_width) < 0)
		    Fatal(fbp, "Error reading raster");
	    }

	    imgstream_fb_close(fbp);
	    fbp = NULL;
	}
    }

    /* Open frame buffer for unbuffered output. */

    if (dst_width == 0)
	dst_width = src_width;		/* default */

    if (dst_height == 0)
	dst_height = src_height;	/* default */

    if ((fbp = imgstream_fb_open(out_fb_file, (size_t)dst_width,
				 (size_t)dst_height)) == NULL)
	Fatal(fbp, "Couldn't open output frame buffer");
    else {
	int wt = (int)imgstream_fb_width(fbp);
	int ht = (int)imgstream_fb_height(fbp);

	/* Use smaller actual frame buffer size for output. */

	V_MIN(dst_width, wt);
	V_MIN(dst_height, ht);

	/* Avoid selecting pixels outside the input image. */

	V_MIN(dst_width, src_width);
	V_MIN(dst_height, src_height);
    }

    /* The following is probably an optimally fast shuffling
     * algorithm; unfortunately, it requires a huge auxiliary array.
     * The way it works is to start with an array of all pixel
     * indices, then repeat: select an entry at random from the array,
     * output that index, replace that entry with the last array
     * entry, then reduce the array size.
     */
    {
	long *loc;		/* keeps track of pixel shuffling */
	long wxh = (long)dst_width * (long)dst_height;
	/* down-counter */

	if ((loc = (long *)malloc((size_t)wxh * sizeof(long))) == NULL)
	    Fatal(fbp, "Not enough memory for location array");

	/* Initialize pixel location array to sequential order. */

	while (--wxh >= 0L)
	    loc[wxh] = wxh;

	/* Select a pixel at random, paint it, and adjust the location
	 * array.
	 */

	for (wxh = (long)dst_width * (long)dst_height; --wxh >= 0L;) {
	    long r = (long)((double)wxh * drand48());
	    long x = loc[r] % dst_width;
	    long y = loc[r] / dst_width;

	    if (imgstream_fb_write(fbp, (int)x, (int)y,
				   in_fb_file == NULL ? bg : pix[x + y * src_width],
				   1) < 0) {
		Fatal(fbp, "Error writing pixel");
	    }

	    loc[r] = loc[wxh];	/* track the shuffle */
	}
    }

    /* Close the frame buffer. */

    imgstream_fb_close(fbp);
    fbp = NULL;

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
