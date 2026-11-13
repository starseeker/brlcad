/*                        P N G - F B . C
 * BRL-CAD
 *
 * Copyright (c) 1998-2026 United States Government as represented by
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
/** @file png2fb.cpp
 *
 * Program to take PNG (Portable Network Graphics) files and send them to a framebuffer.
 *
 */

#include "common.h"

#include <stdlib.h>

#include "bio.h"
#include "bu/getopt.h"
#include "icv.h"
#include "ged.h"
#include "imgstream/fb_compat.h"
#include "ged/display.h"

struct png2fb_state {
    double def_screen_gamma;	/* Don't add more gamma, default = 1.0*/
    /* Particularly because on SGI, the system provides gamma correction,
     * so programs like this one don't have to.
     */
    const char *file_name;
    icv_image_t *image;
    int multiple_lines;	/* Streamlined operation */
    int file_xoff;
    int file_yoff;
    int scr_xoff;
    int scr_yoff;
    int clear;
    int zoom;
    int inverse;	/* Draw upside-down */
    int one_line_only;	/* insist on 1-line writes */
    int verbose;
    int header_only;
};

#define PNG2FB_STATE_INIT_ZERO {1.0, NULL, NULL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

static char png2fb_usage[] = "\
Usage: png2fb [-H -i -c -v -z -1] [-m #lines]\n\
	[-g screen_gamma]\n\
	[-x file_xoff] [-y file_yoff] [-X scr_xoff] [-Y scr_yoff]\n\
	[-S squarescrsize] [file.png]\n";

static int
png2fb_get_args(struct png2fb_state *s, int argc, char **argv)
{
    int c;

    if (!s)
	return 0;

    bu_optind = 1;
    while ((c = bu_getopt(argc, argv, "1m:g:HicvzF:x:y:X:Y:S:W:N:h?")) != -1) {
	switch (c) {
	    case '1':
		s->one_line_only = 1;
		break;
	    case 'm':
		s->multiple_lines = atoi(bu_optarg);
		break;
	    case 'g':
		s->def_screen_gamma = atof(bu_optarg);
		break;
	    case 'H':
		s->header_only = 1;
		break;
	    case 'i':
		s->inverse = 1;
		break;
	    case 'c':
		s->clear = 1;
		break;
	    case 'v':
		s->verbose = 1;
		break;
	    case 'z':
		s->zoom = 1;
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
	    default:		/* '?''h' */
		return 0;
	}
    }

    if (bu_optind >= argc) {
	if (isatty(fileno(stdin)))
	    return 0;
	s->file_name = "-";
	setmode(fileno(stdin), O_BINARY);
    } else {
	s->file_name = argv[bu_optind];
	FILE *probe = fopen(s->file_name, "rb");
	if (!probe) {
	    perror(s->file_name);
	    fprintf(stderr,
		    "png-fb: cannot open \"%s\" for reading\n",
		    s->file_name);
	    return 0;
	}
	fclose(probe);
    }

    if (argc > ++bu_optind)
	fprintf(stderr, "png-fb: excess argument(s) ignored\n");

    return 1;		/* OK */
}


static int
png2fb_apply(struct imgstream_fb *fb, void *userdata)
{
    struct png2fb_state *state = (struct png2fb_state *)userdata;
    if (!state || !state->image)
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
    return imgstream_fb_import_icv(fb, state->image, &options) == 0 ?
	BRLCAD_OK : BRLCAD_ERROR;
}

int
ged_png2fb_core(struct ged *gedp, int argc, const char *argv[])
{
    int ret;

    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    struct png2fb_state p2fbs = PNG2FB_STATE_INIT_ZERO;

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "no current view set\n");
	return BRLCAD_ERROR;
    }

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], png2fb_usage);
	return GED_HELP;
    }

    if (!png2fb_get_args(&p2fbs, argc, (char **)argv)) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], png2fb_usage);
	return GED_HELP;
    }

    p2fbs.image = icv_read(BU_STR_EQUAL(p2fbs.file_name, "-") ? NULL :
	p2fbs.file_name, BU_MIME_IMAGE_PNG, 0, 0);
    if (!p2fbs.image) {
	bu_vls_printf(gedp->ged_result_str, "unable to decode PNG input");
	return BRLCAD_ERROR;
    }
    if (p2fbs.verbose)
	bu_vls_printf(gedp->ged_result_str, "Image size: %zu X %zu\n",
	    p2fbs.image->width, p2fbs.image->height);
    if (p2fbs.header_only) {
	bu_vls_printf(gedp->ged_result_str, "WIDTH=%zu HEIGHT=%zu\n",
	    p2fbs.image->width, p2fbs.image->height);
	icv_destroy(p2fbs.image);
	return BRLCAD_OK;
    }
    /* libpng historically used 0.5 when an input did not provide gAMA.
     * Preserve png2fb's display conversion while libicv supplies decoding. */
    p2fbs.image->gamma_corr = (float)(0.5 * p2fbs.def_screen_gamma);
    ret = ged_view_framebuffer_apply(gedp, view_ctx,
	png2fb_apply, &p2fbs, 1);
    icv_destroy(p2fbs.image);

    if (ret == BRLCAD_OK)
	return BRLCAD_OK;

    bu_vls_printf(gedp->ged_result_str,
	"unable to import PNG data into the active Obol framebuffer");
    return BRLCAD_ERROR;
}


#include "../include/plugin.h"

#define GED_PNG2FB_COMMANDS(X, XID) \
    X(png2fb, ged_png2fb_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_PNG2FB_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_png2fb", 1, GED_PNG2FB_COMMANDS)

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
