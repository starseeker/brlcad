/*                         A S I Z E . C
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
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
/** @addtogroup libstruct fb */
/** @{ */
/** @file libdm/asize.c
 *
 * Legacy framebuffer image-sizing compatibility functions.
 *
 */
/** @} */

#include "common.h"

#include "dm.h"
#include "imgstream/fb_compat.h"


int
fb_common_file_size(size_t *width, size_t *height, const char *filename,
	int bytes_per_pixel)
{
    return imgstream_image_file_size(width, height, filename,
	bytes_per_pixel);
}


int
fb_common_name_size(size_t *width, size_t *height, const char *name)
{
    return imgstream_image_name_size(width, height, name);
}


int
fb_common_image_size(size_t *width, size_t *height, size_t pixel_count)
{
    return imgstream_image_size(width, height, pixel_count);
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
