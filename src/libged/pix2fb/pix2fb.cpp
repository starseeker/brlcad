/*                        P I X 2 F B . C
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
/** @file libged/pix2fb.cpp
 *
 * Program to take bottom-up pixel files and send them to a framebuffer.
 *
 */

#include "common.h"

#include <stdlib.h>
#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif
#ifdef HAVE_WINSOCK_H
#  include <winsock.h>
#endif

#include "bio.h"

#include "bu/getopt.h"
#include "bu/snooze.h"
#include "BObol/BFramebuffer.h"
#include "imgstream/fb_compat.h"

#include "pkg.h"
#include "ged.h"
#include "../ged_bobol_private.hpp"

struct pix2fb_state {
    size_t file_width;	/* default input width */
    size_t file_height;	/* default input height */
    const char *file_name;
    int multiple_lines;	/* Streamlined operation */
    int infd;
    int fileinput;		/* file or pipe on input? */
    int autosize;		/* !0 to autosize input */
    int scr_width;		/* screen tracks file if not given */
    int scr_height;
    int file_xoff;
    int file_yoff;
    int scr_xoff;
    int scr_yoff;
    int clear;
    int zoom;
    int inverse;		/* Draw upside-down */
    int one_line_only;	/* insist on 1-line writes */
    int pause_sec; 	/* Pause that many seconds before
			   closing the FB and exiting */
};

#define PIX2FB_STATE_INIT_ZERO {512, 512, NULL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

static char pix2fb_usage[] = "\
Usage: pix-fb [-a -h -i -c -z -1] [-m #lines] [-F framebuffer]\n\
	[-s squarefilesize] [-w file_width] [-n file_height]\n\
	[-x file_xoff] [-y file_yoff] [-X scr_xoff] [-Y scr_yoff]\n\
	[-S squarescrsize] [-W scr_width] [-N scr_height] [-p seconds]\n\
	[file.pix]\n";

static int
pix2fb_get_args(struct pix2fb_state *s, int argc, char **argv)
{
    int c;

    if (!s)
	return 0;

    bu_optind = 1;
    while ((c = bu_getopt(argc, argv, "1m:aiczF:p:s:w:n:x:y:X:Y:S:W:N:h?")) != -1) {
	switch (c) {
	    case '1':
		s->one_line_only = 1;
		break;
	    case 'm':
		s->multiple_lines = atoi(bu_optarg);
		break;
	    case 'a':
		s->autosize = 1;
		break;
	    case 'i':
		s->inverse = 1;
		break;
	    case 'c':
		s->clear = 1;
		break;
	    case 'z':
		s->zoom = 1;
		break;
	    case 's':
		/* square file size */
		s->file_height = s->file_width = atoi(bu_optarg);
		s->autosize = 0;
		break;
	    case 'w':
		s->file_width = atoi(bu_optarg);
		s->autosize = 0;
		break;
	    case 'n':
		s->file_height = atoi(bu_optarg);
		s->autosize = 0;
		break;
	    case 'x':
		s->file_xoff = atoi(bu_optarg);
		break;
	    case 'y':
		s->file_yoff = atoi(bu_optarg);
		break;
	    case 'X':
		s->scr_xoff = atoi(bu_optarg);
		break;
	    case 'Y':
		s->scr_yoff = atoi(bu_optarg);
		break;
	    case 'S':
		s->scr_height = s->scr_width = atoi(bu_optarg);
		break;
	    case 'W':
		s->scr_width = atoi(bu_optarg);
		break;
	    case 'N':
		s->scr_height = atoi(bu_optarg);
		break;
	    case 'p':
		s->pause_sec = atoi(bu_optarg);
		break;

	    default:		/* 'h' '?' */
		return 0;
	}
    }

    if (bu_optind >= argc) {
	if (isatty(fileno(stdin)))
	    return 0;
	s->file_name = "-";
	s->infd = 0;
	setmode(s->infd, O_BINARY);
    } else {
	s->file_name = argv[bu_optind];
	s->infd = open(s->file_name, 0);
	if (s->infd < 0) {
	    perror(s->file_name);
	    fprintf(stderr,
			  "pix-fb: cannot open \"%s\" for reading\n",
			  s->file_name);
	    bu_exit(1, NULL);
	}
	setmode(s->infd, O_BINARY);
	s->fileinput++;
    }

    if (argc > ++bu_optind)
	fprintf(stderr, "pix-fb: excess argument(s) ignored\n");

    return 1;		/* OK */
}


static int
pix2fb_apply(BObolFramebufferStream &stream, void *userdata)
{
    struct pix2fb_state *state = (struct pix2fb_state *)userdata;
    if (!state)
	return BRLCAD_ERROR;
    struct imgstream_fb_import_options options =
	IMGSTREAM_FB_IMPORT_OPTIONS_INIT;
    options.file_xoff = state->file_xoff;
    options.file_yoff = state->file_yoff;
    options.screen_xoff = state->scr_xoff;
    options.screen_yoff = state->scr_yoff;
    options.clear = state->clear;
    options.zoom = state->zoom;
    options.inverse = state->inverse;
    return imgstream_fb_import_pix_fd(stream.framebuffer(), state->infd, state->file_name,
	state->file_width, state->file_height, state->autosize,
	&options) == 0 ? BRLCAD_OK : BRLCAD_ERROR;
}

int
ged_pix2fb_core(struct ged *gedp, int argc, const char *argv[])
{
    int ret;

    struct pix2fb_state p2fbs = PIX2FB_STATE_INIT_ZERO;

    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, ": no current view set\n");
	return BRLCAD_ERROR;
    }

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], pix2fb_usage);
	return GED_HELP;
    }

    if (!pix2fb_get_args(&p2fbs, argc, (char **)argv)) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], pix2fb_usage);
	return GED_HELP;
    }

    ret = ged_bobol_framebuffer_apply(gedp, view_ctx,
	pix2fb_apply, &p2fbs, true);

    if (p2fbs.infd != 0)
	close(p2fbs.infd);

    bu_snooze(BU_SEC2USEC(p2fbs.pause_sec));


    if (ret == BRLCAD_OK)
	return BRLCAD_OK;

    bu_vls_printf(gedp->ged_result_str,
	"unable to import PIX data into the active Obol framebuffer");

    return BRLCAD_ERROR;
}


#include "../include/plugin.h"

#define GED_PIX2FB_COMMANDS(X, XID) \
    X(pix2fb, ged_pix2fb_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_PIX2FB_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_pix2fb", 1, GED_PIX2FB_COMMANDS)

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
