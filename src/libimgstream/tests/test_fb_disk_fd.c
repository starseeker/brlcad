/*                   T E S T _ F B _ D I S K _ F D . C
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

#include "bio.h"
#include "bu/app.h"
#include "bu/file.h"
#include "imgstream/fb_compat.h"


int
main(int argc, char **argv)
{
    char path[MAXPATHLEN] = {0};
    FILE *fp;
    imgstream_fb_t *fbp;

    bu_setprogname(argv[0]);
    (void)argc;

    fp = bu_temp_file(path, sizeof(path));
    if (fp == NULL || fclose(fp) != 0 || !bu_file_delete(path)) {
	fprintf(stderr, "failed to prepare disk framebuffer test path\n");
	return 1;
    }
    if (close(fileno(stdin)) != 0) {
	fprintf(stderr, "failed to release descriptor zero\n");
	return 1;
    }

    fbp = imgstream_fb_open(path, 1, 1);
    if (fbp == NULL) {
	fprintf(stderr, "disk framebuffer rejected descriptor zero\n");
	return 1;
    }
    imgstream_fb_close(fbp);
    if (!bu_file_delete(path)) {
	fprintf(stderr, "failed to close or remove descriptor-zero framebuffer\n");
	return 1;
    }

    return 0;
}
