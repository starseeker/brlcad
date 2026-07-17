/*                         S C R E E N G R A B . C
 * BRL-CAD
 *
 * Copyright (c) 2008-2026 United States Government as represented by
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
/** @file libged/screengrab.c
 *
 * The screengrab command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "BObol/BDisplayEndpoint.h"
#include "icv.h"
#include "ged/draw_obol.h"

#include "../ged_private.h"

static int
screengrab_image_mime(struct bu_vls *msg, size_t argc, const char **argv, void *set_mime)
{
    int type_int;
    bu_mime_image_t type = BU_MIME_IMAGE_UNKNOWN;
    bu_mime_image_t *set_type = (bu_mime_image_t *)set_mime;

    BU_OPT_CHECK_ARGV0(msg, argc, argv, "mime format");

    type_int = bu_file_mime(argv[0], BU_MIME_IMAGE);
    type = (type_int < 0) ? BU_MIME_IMAGE_UNKNOWN : (bu_mime_image_t)type_int;
    if (type == BU_MIME_IMAGE_UNKNOWN) {
	if (msg) {
	    bu_vls_sprintf(msg, "Error - unknown geometry file type: %s \n", argv[0]);
	}
	return -1;
    }
    if (set_type) {
	(*set_type) = type;
    }
    return 1;
}

int
ged_screen_grab_core(struct ged *gedp, int argc, const char *argv[])
{
    int print_help = 0;
    int grab_fb = 0;
    unsigned char *idata = NULL;
    struct icv_image *bif = NULL;	/**< icv image container for saving images */
    struct bu_vls view_name = BU_VLS_INIT_ZERO;
    bu_mime_image_t type = BU_MIME_IMAGE_AUTO;
    static char usage[] = "Usage: screengrab [-h] [-F] [-V view] [--format fmt] [file.img]\n";

    struct bu_opt_desc d[5];
    BU_OPT(d[0], "h", "help",           "",     NULL,             &print_help,       "Print help and exit");
    BU_OPT(d[1], "F", "fb",             "",     NULL,             &grab_fb,          "screengrab framebuffer instead of scene display");
    BU_OPT(d[2], "V", "view",           "name", &bu_opt_vls,      &view_name,        "view endpoint to capture");
    BU_OPT(d[3], "",  "format",         "fmt",  &screengrab_image_mime,      &type,             "output image file format");
    BU_OPT_NULL(d[4]);

    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    GED_CHECK_DRAWABLE(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    argc-=(argc>0); argv+=(argc>0); /* done with command name argv[0] */

    int opt_ret = bu_opt_parse(NULL, argc, argv, d);

    if (print_help) {
	_ged_cmd_help(gedp, usage, d);
	return GED_HELP;
    }

    argc = opt_ret;

    void *view_ctx = bu_vls_strlen(&view_name) ?
	ged_view_find_ctx(gedp, bu_vls_cstr(&view_name)) :
	ged_view_active_ctx(gedp);
    bu_vls_free(&view_name);

    /* must be wanting help */
    if (!argc) {
	_ged_cmd_help(gedp, usage, d);
	return GED_HELP;
    }

    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "view endpoint not found\n");
	return BRLCAD_ERROR;
    }
    bobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(view_ctx);
    if (!endpoint) {
	bu_vls_printf(gedp->ged_result_str,
		"view '%s' has no Obol display endpoint\n",
		bv_name_get(bv_context_view_const(
		    (const struct bv_context *)view_ctx)));
	return BRLCAD_ERROR;
    }

    (void)ged_draw_obol_framebuffer_present(gedp);
    if (!grab_fb && !bobol_display_endpoint_view_sync(endpoint, view_ctx)) {
	bu_vls_printf(gedp->ged_result_str,
		"%s: could not synchronize the endpoint camera.", argv[0]);
	return BRLCAD_ERROR;
    }
    size_t image_size = 0;
    unsigned int width = 0;
    unsigned int height = 0;
    unsigned int components = 0;
    enum bobol_capture_plane capture_plane = grab_fb ?
	BOBOL_CAPTURE_FRAMEBUFFER : BOBOL_CAPTURE_COMPOSITE;
    if (!bobol_display_endpoint_capture_plane(endpoint, capture_plane,
	    &idata, &image_size, &width, &height, &components) || !idata ||
	    !width || !height ||
	    (components != 3 && components != 4) ||
	    image_size < (size_t)width * height * components) {
	bu_vls_printf(gedp->ged_result_str,
		"%s: Obol display endpoint did not return %s image data.",
		argv[0], grab_fb ? "framebuffer" : "composite");
	if (idata)
	    bu_free(idata, "endpoint image data");
	return BRLCAD_ERROR;
    }

    bif = icv_create((int)width, (int)height, ICV_COLOR_SPACE_RGB);
    if (!bif) {
	bu_free(idata, "endpoint image data");
	bu_vls_printf(gedp->ged_result_str,
		": could not create icv_image write structure.");
	return BRLCAD_ERROR;
    }
    unsigned char *row = (unsigned char *)bu_malloc((size_t)width * 3,
	"endpoint image row");
    for (unsigned int y = 0; y < height; y++) {
	for (unsigned int x = 0; x < width; x++) {
	    const unsigned char *source = idata +
		((size_t)y * width + x) * components;
	    row[(size_t)x * 3] = source[0];
	    row[(size_t)x * 3 + 1] = source[1];
	    row[(size_t)x * 3 + 2] = source[2];
	}
	icv_writeline(bif, (int)y, row, ICV_DATA_UCHAR);
    }
    bu_free(row, "endpoint image row");
    bu_free(idata, "endpoint image data");

    icv_write(bif, argv[0], type);
    icv_destroy(bif);

    return BRLCAD_OK;
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
