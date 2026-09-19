/*                        F B H E L P . C
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
/** @file fbhelp.c
 *
 * Print out info about the selected image-stream target.
 *
 */

#include "common.h"

#include <stdlib.h>

#include "bio.h"
#include "bu/app.h"
#include "bu/getopt.h"
#include "imgstream/fb_compat.h"


static char *framebuffer = NULL;

static char usage[] = "\
Usage: fbhelp [-F framebuffer]\n";

int
main(int argc, char **argv)
{
    int c;
    imgstream_fb_t *fbp;

    bu_setprogname(argv[0]);

    while ((c = bu_getopt(argc, argv, "F:h?")) != -1) {
	switch (c) {
	    case 'F':
		framebuffer = bu_optarg;
		break;
	    default:		/* '?' */
		(void)fputs(usage, stderr);
		return 1;
	}
    }
    if (argc > bu_optind) {
	fprintf(stderr, "fbhelp: excess argument(s) not supported\n");
	(void)fputs(usage, stderr);
	return 1;
    }
    fprintf(stdout, "Image-stream targets are selected with -F.\n"
	"Memory, file, remote, and diagnostic streams are GUI-independent.\n"
	"Toolkit display targets require an application-owned host.\n");

    fprintf(stdout, "=============== Available Targets ================\n");
    fprintf(stdout, "memory: /dev/mem\nfile: path or /dev/disk:path\n");
    fprintf(stdout, "remote: host:port, tcp:host:port, or ipc:address\n");
    fprintf(stdout, "diagnostic: /dev/null, /dev/debug, or /dev/txt\n");

    fprintf(stdout, "=============== Current Selection ================\n");
    if ((fbp = imgstream_fb_open(framebuffer, 0, 0)) == NULL) {
	fprintf(stderr, "fbhelp: Can't open frame buffer\n");
	return 1;
    }
    fprintf(stdout, "target: %s\nkind: %s\nsize: %lux%lu\n",
	imgstream_fb_name(fbp),
	imgstream_fb_spec_kind_name(imgstream_fb_spec_kind(imgstream_fb_name(fbp))),
	(unsigned long)imgstream_fb_width(fbp),
	(unsigned long)imgstream_fb_height(fbp));
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
